// input.cpp
#include <vector>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>
#include "input.h"


namespace GooFit{

GooPdf *makeDstar_veto() {
   
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 0.), Variable("Fiducial_veto12_max", s12.getLowerLimit()), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("DFiducial_veto13_min", 0.), Variable("Fiducial_veto13_max", s13.getLowerLimit()), PAIR_13);
    
    std::vector<VetoInfo> vetos;
    //vetos.push_back(Fiducial_veto12);
    //vetos.push_back(Fiducial_veto13);

    DalitzVetoPdf* vetoPdf = new DalitzVetoPdf("Dstar_veto", s12, s13, 
    						Variable("Mother_Mass",D_MASS), Variable("Daughter1_Mass",d1_MASS), 
    						Variable("Daughter2_Mass",d2_MASS), Variable("Daughter3_Mass",d3_MASS),vetos);

    return vetoPdf;
}



DalitzPlotPdf* makesignalpdf(GooPdf* eff = 0){

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    ResonancePdf* sigma		= new Resonances::POLE("sigma",
							Variable("v_sigma_real",1.),
							Variable("v_sigma_img",0.),
							Variable("v_sigma__pole_real",0.47),
							Variable("v_sigma_pole_img",0.22),0,PAIR_12,true);
    
	ResonancePdf* f0_980	= new Resonances::FLATTE("f0_980",
							Variable("f0_980_real",1.),
							Variable("f0_980_img",0.),
							Variable("f0_980_mass",0.990),
							Variable("f0_980_gpp",0.02),
							Variable("f0_980_gkk",0.09),PAIR_12,true);

	ResonancePdf* f0_980_BW	= new Resonances::RBW("f0_980",
							Variable("f0_980_real",1.),
							Variable("f0_980_img",0.),
							Variable("f0_980_mass",0.990),
							Variable("f0_980_width",0.04),int(0),PAIR_12,true);

	ResonancePdf* f0_1370	= new Resonances::RBW("f0_1370",
							Variable("f0_1370_real",1.),
							Variable("f0_1370_img",0.),
							Variable("f0_1370_mass",1.370),
							Variable("f0_1370_width",0.3),int(0),PAIR_12,true);
	
	ResonancePdf* f0_1500	= new Resonances::RBW("f0_1500",
							Variable("f0_1500_real",1.),
							Variable("f0_1500_img",0.),
							Variable("f0_1500_mass",1.505),
							Variable("f0_1500_width",0.109),int(0),PAIR_12,true);
	
	ResonancePdf* omega_782	= new Resonances::RBW("omega_782",
							Variable("omega_782_real",1.),
							Variable("omega_782_img",0.),
							Variable("omega_782_mass",0.782),
							Variable("omega_782_width",0.0085),int(1),PAIR_12,true);
	
    ResonancePdf* rho_770	= new Resonances::RBW("rho_770",
							Variable("rho_770_real",1.),
							Variable("rho_770_img",0.),
							Variable("rho_770_mass",.7755),
							Variable("rho_770_width",0.149),int(1),PAIR_12,true);
	
	ResonancePdf* rho_1450	= new Resonances::RBW("rho_1450",
							Variable("rho_1450_real",1.),
							Variable("rho_1450_img",0.),
							Variable("rho_1450_mass",1.465),
							Variable("rho_1450_width",0.4),int(1),PAIR_12,true);
	
	ResonancePdf* f2_1270	= new Resonances::RBW("f2_1270",
							Variable("f2_1270_real",1.),
							Variable("f2_1270_img",0.),
							Variable("f2_1270_mass",1.2751),
							Variable("f2_1270_width",0.185),int(2),PAIR_12,true);
	
	ResonancePdf* f2_1525	= new Resonances::RBW("f2_1525",
							Variable("f2_1525_real",1.),
							Variable("f2_1525_img",0.),
							Variable("f2_1525_mass",1.525),
							Variable("f2_1525_width",0.073),int(2),PAIR_12,true);
		
  
    ResonancePdf *nonr 		= new Resonances::NonRes("nonr",
							Variable("nonr_real",1.), 
							Variable("nonr_img",0));

 
    //Pushing Resonances 
    //dtoppp.resonances.push_back(sigma);
    //dtoppp.resonances.push_back(f0_980);
    dtoppp.resonances.push_back(f0_980_BW);
    //dtoppp.resonances.push_back(f0_1500);
    dtoppp.resonances.push_back(f0_1370);
    //dtoppp.resonances.push_back(omega_782); 
    dtoppp.resonances.push_back(rho_770); 
    dtoppp.resonances.push_back(rho_1450);
    dtoppp.resonances.push_back(f2_1270);
    //dtoppp.resonances.push_back(f2_1525);
    //dtoppp.resonances.push_back(nonr);


    if(!eff) {
        // By default create a constant efficiency.
        std::vector<Variable> offsets;
        std::vector<Observable> observables;
        std::vector<Variable> coefficients;
        Variable constantOne("constantOne", 1);
        Variable constantZero("constantZero", 0);
        observables.push_back(s12);
        observables.push_back(s13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0); //No efficiency
    }

    return new DalitzPlotPdf("signalPDF", s12, s13, eventNumber, dtoppp, eff);
}

}

