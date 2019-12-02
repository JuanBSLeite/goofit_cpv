// ROOT stuff
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TComplex.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TH2Poly.h>

//Minuit

#include <Minuit2/MnStrategy.h>

// System stuff
#include <CLI/Timer.hpp>
#include <fstream>
#include <string>
#include <time.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/fitting/FitManagerMinuit2.h>
#include <goofit/fitting/FitManagerMinuit1.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/combine/EventWeightedAddPdf.h>
#include <goofit/PDFs/combine/CompositePdf.h>
#include <goofit/PDFs/physics/IncoherentSumPdf.h>

#include <thrust/transform_reduce.h>

//including input model
#include "input.cpp"

using namespace std;
using namespace GooFit;
using namespace ROOT;

UnbinnedDataSet *toyMC = nullptr;
DalitzPlotPdf* signalpdf = nullptr;
GooPdf* backgroundpdf = nullptr;
GooPdf* efficiency = nullptr;
GooPdf* Veto = nullptr;

// PWA INPUT FILE NAME

// Data File
const string bkghisto_file = "../../../dados/bkg_histo_16.root";
const string effhisto_file = "../../../dados/eff_16.root";
const string bkghisto_name = "h_eff";
const string effhisto_name = "h_eff";


GooPdf* makeEfficiencyPdf() {

    vector<Observable> lvars;
    lvars.push_back(s12);
    lvars.push_back(s13);
    BinnedDataSet *binEffData = new BinnedDataSet(lvars);
    
    TFile *f     = TFile::Open(effhisto_file.c_str());
    auto effHistogram = (TH2F *)f->Get(effhisto_name.c_str()); 

    for(int i = 0; i < s12.getNumBins(); ++i) {
        s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / s12.getNumBins());
        for(int j = 0; j < s13.getNumBins(); ++j) {
            s13.setValue(s13.getLowerLimit() + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / s13.getNumBins());
            if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
            double weight = effHistogram->GetBinContent(effHistogram->FindBin(s12.getValue(), s13.getValue()));
            binEffData->addWeightedEvent(weight);

        }
    }


    // Smooth
    Variable *effSmoothing = new Variable("effSmoothing_eff", 0.0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, *effSmoothing);
    return ret;
}

GooPdf* makeBackgroundPdf() {

   BinnedDataSet *binBkgData = new BinnedDataSet({s12, s13});
    
    TFile *f = new TFile(bkghisto_file.c_str());
    auto bkgHistogram = (TH2F*)f->Get(bkghisto_name.c_str());

    for(int i = 0; i < s12.getNumBins(); ++i) {
        s12.setValue(s12.getLowerLimit() + (s12.getUpperLimit() - s12.getLowerLimit()) * (i + 0.5) / s12.getNumBins());
        for(int j = 0; j < s13.getNumBins(); ++j) {
            s13.setValue(s13.getLowerLimit() + (s13.getUpperLimit() - s13.getLowerLimit()) * (j + 0.5) / s13.getNumBins());
            if(!inDalitz(s12.getValue(), s13.getValue(), D_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
            double weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(s12.getValue(), s13.getValue()));
            binBkgData->addWeightedEvent(weight);

        }
    }

    Variable *effSmoothing  = new Variable("effSmoothing_bkg",0.);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("background", binBkgData, *effSmoothing);

    return ret;
}

GooPdf *makeDstar_veto() {
   
    VetoInfo Dstar_veto12(Variable("Dstar_veto12_min", 2.85), Variable("Dstar_veto12_max", s12_max), PAIR_12);
    VetoInfo Dstar_veto13(Variable("Dstar_veto13_min", 2.85), Variable("Dstar_veto13_max", s12_max), PAIR_13);
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 0.), Variable("Fiducial_veto12_max", s12.getLowerLimit()), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("DFiducial_veto13_min", 0.), Variable("Fiducial_veto13_max", s13.getLowerLimit()), PAIR_13);
    
    vector<VetoInfo> vetos;
    vetos.push_back(Dstar_veto12);
    vetos.push_back(Dstar_veto13);
    vetos.push_back(Fiducial_veto12);
    vetos.push_back(Fiducial_veto13);

    DalitzVetoPdf* Dstar_veto = new DalitzVetoPdf("Dstar_veto", s12, s13, 
        Variable("Mother_Mass",D_MASS), Variable("Daughter1_Mass",d1_MASS), 
        Variable("Daughter2_Mass",d2_MASS), Variable("Daughter3_Mass",d3_MASS), 
        vetos);

    return Dstar_veto;
}


DalitzPlotPdf* NR_DP(GooPdf* eff = 0){

	  DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

        Variable nonr_real("nonr_real",1., 0.01,0,0);
      Variable nonr_imag("nonr_imag",0., 0.01,0,0);


       ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_real, nonr_imag);
       dtoppp.resonances.push_back(nonr);
     
       if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;
        Variable constantOne("constantOne", 1);
        Variable constantZero("constantZero", 0);
        observables.push_back(s12);
        observables.push_back(s13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); //No efficiency
    }

       return new DalitzPlotPdf("NRPDF", s12, s13, eventNumber, dtoppp, eff);

}



std::vector<std::vector<fptype>> fractions(DalitzPlotPdf* signalpdf){

    cout << "Fit Fractions" << endl;
    auto NRPDF = NR_DP(0);
    ProdPdf Prod_NRPDF{"prodNRPDF",{NRPDF}};
    DalitzPlotter dp_NR(&Prod_NRPDF,NRPDF);
    auto flatMC = new UnbinnedDataSet({s12,s13,eventNumber});
    dp_NR.fillDataSetMC(*flatMC,100000);

    ProdPdf overallsignal{"overallsignal",{signalpdf}};
    overallsignal.setData(flatMC);
    signalpdf->setDataSize(flatMC->getNumEvents());


    s12.setNumBins(400);
    s13.setNumBins(400);
    
    signalpdf->setParameterConstantness(true); 
    FitManagerMinuit2 fitter(&overallsignal);
    fitter.setVerbosity(0);
    fitter.fit();
    
    return signalpdf->fit_fractions();
   
}

void to_root(UnbinnedDataSet* toyMC , std::string name ){


    double _s12, _s13,_s23;
    TFile *f = new TFile(name.c_str(),"recreate");
    TTree *t = new TTree("DecayTree","toyMC");
    TBranch *b_s12 = t->Branch("s12",&_s12,"s12/D");
    TBranch *b_s13 = t->Branch("s13",&_s13,"s13/D");
    TBranch *b_s23 = t->Branch("s23",&_s23,"s23/D");
    for(int i = 0; i < toyMC->getNumEvents(); i++){
		toyMC->loadEvent(i);
		t->GetEntry(i);
		_s12 = s12.getValue();
		_s13 = s13.getValue();
		_s23 = POW2(D_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS) - s12.getValue() - s13.getValue() ;
		t->Fill();
    }
	t->Write();
	f->Write();
	f->Close();

}

void gentoyMC(std::string name, size_t nevents){
    
    Veto = makeDstar_veto();
    efficiency = makeEfficiencyPdf();
    ProdPdf *effWithVeto = new ProdPdf("effWithVeto", {Veto,efficiency});
    signalpdf = makesignalpdf(effWithVeto); 
 
    backgroundpdf = makeBackgroundPdf();
    backgroundpdf->setParameterConstantness(true);
    ProdPdf bkgWithVeto{"bkgWithVeto",{backgroundpdf,Veto}}; 
    
    Variable frac("Signal_Purity",0.925);
    vector<Variable> weights;
    weights.push_back(frac);

    vector<PdfBase*> comps = {signalpdf,&bkgWithVeto};
    AddPdf* overallPdf = new AddPdf("overallPdf",weights,comps);

    DalitzPlotter dp(overallPdf,signalpdf);
    toyMC = new UnbinnedDataSet({s12,s13,eventNumber});
    dp.fillDataSetMC(*toyMC,nevents);

    to_root(toyMC,name);
    
    fractions(signalpdf); 

}

int main(int argc, char **argv){

    GooFit::Application app{"D2hhh_cpv",argc,argv};

    std::string name = "MC/toyMC.root";
    size_t events = 1000000;
    auto gen = app.add_subcommand("gen","fit toy data");
    gen->add_option("-n,--name",name,"output file name");
    gen->add_option("-e,--events",events,"number of events");

    GOOFIT_PARSE(app);

    /// Make the mc directory if it does not exist
    std::string command = "mkdir -p MC";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `MC` directory failed");

    
    if(*gen){
        CLI::AutoTimer timer("Gen");
        gentoyMC(name,events);
    }
}