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
DalitzPlotPdf* NRPDF = nullptr;
GooPdf* backgroundpdf = nullptr;
GooPdf* efficiency = nullptr;
UnbinnedDataSet* Data = nullptr;

const string effhisto_file = eff_file;
const string effhisto_name = eff_name;
vector<fptype> HH_bin_limits;
vector<Variable> pwa_coefs_amp;
vector<Variable> pwa_coefs_phs;


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
            if(!inDalitz(s12.getValue(), s13.getValue(), Mother_MASS, d1_MASS, d2_MASS, d3_MASS)){continue;}
            double weight = effHistogram->GetBinContent(effHistogram->FindBin(s12.getValue(), s13.getValue()));
            binEffData->addWeightedEvent(weight);

        }
    }


    // Smooth
    Variable *effSmoothing = new Variable("effSmoothing_eff", 0.0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, *effSmoothing);
    return ret;
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
		_s23 = POW2(Mother_MASS) + POW2(d1_MASS) + POW2(d2_MASS) + POW2(d3_MASS) - s12.getValue() - s13.getValue() ;
		t->Fill();
    }
	t->Write();
	f->Write();
	f->Close();

}

void gentoyMC(std::string name, size_t nevents){
    
    s12.setNumBins(1500);
    s13.setNumBins(1500);
    ProdPdf *eff = nullptr;
 
    if(effOn){  
    	efficiency = makeEfficiencyPdf();
  	eff = new ProdPdf("eff", {efficiency});
    	signalpdf = makesignalpdf(eff); 
    }else{ 
    	signalpdf = makesignalpdf(0); 
    }

    
    ProdPdf *ProdSignal = new ProdPdf("ProdSignal",{signalpdf});
    signalpdf->setParameterConstantness(true);
    
    DalitzPlotter dp(ProdSignal,signalpdf);
    toyMC = new UnbinnedDataSet({s12,s13,eventNumber});
    dp.fillDataSetMC(*toyMC,nevents);

    TCanvas foo;
    auto h2 = dp.make2D();
    h2->Rebin2D(10,10);
    h2->Draw("colz"); foo.SaveAs("MC/DP.png");
    auto h1 = (TH1F*)h2->ProjectionX("s12");
    h1->Draw("HIST"); foo.SaveAs("MC/s12.png");


    ProdSignal->setData(toyMC);
    signalpdf->setDataSize(toyMC->getNumEvents());
    
    FitManager fitter(ProdSignal);
    fitter.setVerbosity(0);
    fitter.fit();
    
    signalpdf->fit_fractions();

    to_root(toyMC,name);
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
