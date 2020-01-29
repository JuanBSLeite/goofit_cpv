// input.cpp
#include <vector>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>
#include "input.h"


namespace GooFit{


double Signal_Purity = 0.9;


//include veto
GooPdf *makeDstar_veto() {
   								//Veto range from a to b
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 2.85), Variable("Fiducial_veto12_max", s12.getUpperLimit()), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("DFiducial_veto13_min", 2.9), Variable("Fiducial_veto13_max", s13.getUpperLimit()), PAIR_13);
    
    std::vector<VetoInfo> vetos;
    vetos.push_back(Fiducial_veto12);
    //vetos.push_back(Fiducial_veto13);

    DalitzVetoPdf* vetoPdf = new DalitzVetoPdf("Dstar_veto", s12, s13, 
    						Variable("Mother_Mass",Mother_MASS), Variable("Daughter1_Mass",d1_MASS), 
    						Variable("Daughter2_Mass",d2_MASS), Variable("Daughter3_Mass",d3_MASS),vetos);

    return vetoPdf;
}



DalitzPlotPdf* makesignalpdf(GooPdf* eff = 0){

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = Mother_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    //Masses and Widths from PDG

    ResonancePdf* sigma		= new Resonances::POLE("sigma",
							Variable("sigma_real",1.),
							Variable("sigma_img",0.),
							Variable("sigma_pole_real",0.47,.001,0,0),
							Variable("sigma_pole_img",0.22,.001,0,0),0,PAIR_13,symdp);
    
	ResonancePdf* kappa         = new Resonances::POLE("kappa",
                                                        Variable("kappa_real",1.,.001,0,0),
                                                        Variable("kappa_img",0.,.001,0,0),
                                                        Variable("kappa_pole_real",0.71),
                                                        Variable("kappa_pole_img",-0.31),0,PAIR_13,symdp);

    ResonancePdf* kappa_BW	= new Resonances::RBW("kappa",
							Variable("kappa_real",1.,.001,0,0),
							Variable("kappa_img",0.,.001,0,0),
							Variable("kappa_mass",0.797),
							Variable("kappa_width",0.410),int(0),PAIR_13,symdp);

	ResonancePdf* f0_980	= new Resonances::FLATTE("f0_980",
							Variable("f0_980_real",1.,.001,0,0),
							Variable("f0_980_img",0.,.001,0,0),
							Variable("f0_980_mass",0.965+0.02),
							Variable("f0_980_gpp",0.165),
							Variable("f0_980_gkk",0.695),PAIR_12,symdp);

	ResonancePdf* f0_980_BW	= new Resonances::RBW("f0_980",
							Variable("f0_980_real",1.,.001,0,0),
							Variable("f0_980_img",0.,.001,0,0),
							Variable("f0_980_mass",0.965+0.02),
							Variable("f0_980_width",0.04),int(0),PAIR_13,symdp);

	ResonancePdf* f0_1370	= new Resonances::RBW("f0_1370",
							Variable("f0_1370_real",1.,.001,0,0),
							Variable("f0_1370_img",0.,.001,0,0),
							Variable("f0_1370_mass",1.435),
							Variable("f0_1370_width",0.135),int(0),PAIR_12,symdp);
    
    ResonancePdf* a0_1450	= new Resonances::RBW("a0_1450",
							Variable("a0_1450_real",1.,.001,0,0),
							Variable("a0_1450_img",0.,.001,0,0),
							Variable("a0_1450_mass",1.474),
							Variable("a0_1450_width",0.265),int(0),PAIR_12,symdp);
	
	ResonancePdf* f0_1500	= new Resonances::RBW("f0_1500",
							Variable("f0_1500_real",1.,.001,0,0),
							Variable("f0_1500_img",0.,.001,0,0),
							Variable("f0_1500_mass",1.505),
							Variable("f0_1500_width",0.109),int(0),PAIR_12,symdp);

	ResonancePdf* ks_1680	= new Resonances::RBW("ks_1680",
							Variable("ks_1680_real",1.,.001,0,0),
							Variable("ks_1680_img",0.,.001,0,0),
							Variable("ks_1680_mass",1.717),
							Variable("ks_1680_width",0.322),int(1),PAIR_13,symdp);

    ResonancePdf* f0_1710	= new Resonances::RBW("f0_1710",
							Variable("f0_1710_real",1.,.001,0,0),
							Variable("f0_1710_img",0.,.001,0,0),
							Variable("f0_1710_mass",1.723),
							Variable("f0_1710_width",0.139),int(0),PAIR_12,symdp);

    ResonancePdf* ks_1430      = new Resonances::RBW("ks_1430",
                                                        Variable("ks_1430_real",1.,.001,0,0),
                                                        Variable("ks_1430_img",0.,.001,0,0),
                                                        Variable("ks_1430_mass",1.425),
                                                        Variable("ks_1430_width",0.270),int(0),PAIR_13,symdp);


	ResonancePdf* omega_782	= new Resonances::RBW("omega_782",
							Variable("omega_782_real",1.,.001,0,0),
							Variable("omega_782_img",0.,.001,0,0),
							Variable("omega_782_mass",0.782),
							Variable("omega_782_width",0.0085),int(1),PAIR_12,symdp);
	
    ResonancePdf* rho_770	= new Resonances::RBW("rho_770",
							Variable("rho_770_real",1.,.001,0,0),
							Variable("rho_770_img",0.,.001,0,0),
							Variable("rho_770_mass",.7755),
							Variable("rho_770_width",0.149),int(1),PAIR_12,symdp);
	
	ResonancePdf* rho_1450	= new Resonances::RBW("rho_1450",
							Variable("rho_1450_real",1.,.001,0,0),
							Variable("rho_1450_img",0.,.001,0,0),
							Variable("rho_1450_mass",1.465),
							Variable("rho_1450_width",0.4),int(1),PAIR_12,symdp);

    ResonancePdf* rho_1700	= new Resonances::RBW("rho_1700",
							Variable("rho_1700_real",1.,.001,0,0),
							Variable("rho_1700_img",0.,.001,0,0),
							Variable("rho_1700_mass",1.740),
							Variable("rho_1700_width",0.1872),int(1),PAIR_12,symdp);

    ResonancePdf* ks_892	  = new Resonances::RBW("ks_892",
                                                        Variable("ks_892_real",1.),
                                                        Variable("ks_892_img",0.),
                                                        Variable("ks_892_mass",0.89581),
                                                        Variable("ks_892_width",0.0474),int(1),PAIR_13,symdp);


    ResonancePdf* ks_1410     = new Resonances::RBW("ks_1410",
                                                        Variable("ks_1410_real",1.,.001,0,0),
                                                        Variable("ks_1410_img",0.,.001,0,0),
                                                        Variable("ks_1410_mass",1.414),
                                                        Variable("ks_1410_width",0.232),int(1),PAIR_13,symdp);

    ResonancePdf* phi_1020     = new Resonances::RBW("phi_1020",
                                                        Variable("phi_1020_real",1.,.001,0,0),
                                                        Variable("phi_1020_img",0.,.001,0,0),
                                                        Variable("phi_1020_mass",1.019461),
                                                        Variable("phi_1020_width",0.004266),int(1),PAIR_12,symdp);


    ResonancePdf* phi_1680     = new Resonances::RBW("phi_1680",
                                                        Variable("phi_1680_real",1.,.001,0,0),
                                                        Variable("phi_1680_img",0.,.001,0,0),
                                                        Variable("phi_1680_mass",1.680),
                                                        Variable("phi_1680_width",0.150),int(1),PAIR_12,symdp);

	ResonancePdf* f2_1270	= new Resonances::RBW("f2_1270",
							Variable("f2_1270_real",1.,.001,0,0),
							Variable("f2_1270_img",0.,.001,0,0),
							Variable("f2_1270_mass",1.2751),
							Variable("f2_1270_width",0.185),int(2),PAIR_12,symdp);

    ResonancePdf* a2_1320	= new Resonances::RBW("a2_1320",
							Variable("a2_1320_real",1.,.001,0,0),
							Variable("a2_1320_img",0.,.001,0,0),
							Variable("a2_1320_mass",1.3183),
							Variable("a2_1320_width",0.107),int(2),PAIR_12,symdp);
	
	ResonancePdf* f2_1525	= new Resonances::RBW("f2_1525",
							Variable("f2_1525_real",1.,.001,0,0),
							Variable("f2_1525_img",0.,.001,0,0),
							Variable("f2_1525_mass",1.525),
							Variable("f2_1525_width",0.073),int(2),PAIR_12,symdp);
	
	ResonancePdf* k2s_1430   = new Resonances::RBW("k2s_1430",
                                                        Variable("k2s_1430_real",1.,.001,0,0),
                                                        Variable("k2s_1430_img",0.,.001,0,0),
                                                        Variable("k2s_1430_mass",1.4324),
                                                        Variable("k2s_1430_width",0.109),int(2),PAIR_13,symdp);	
  
    ResonancePdf *nonr 		= new Resonances::NonRes("nonr",
							Variable("nonr_real",1.,.001,0,0), 
							Variable("nonr_img",0,.001,0,0));


//Push Resonances of your model 
    //dtoppp.resonances.push_back(sigma);
    //dtoppp.resonances.push_back(kappa);
    dtoppp.resonances.push_back(kappa_BW);
    dtoppp.resonances.push_back(f0_980);
    //dtoppp.resonances.push_back(f0_980_BW);
    dtoppp.resonances.push_back(a0_1450);
    dtoppp.resonances.push_back(f0_1500);
    dtoppp.resonances.push_back(ks_1680);
    dtoppp.resonances.push_back(f0_1710);
    dtoppp.resonances.push_back(f0_1370);
    dtoppp.resonances.push_back(ks_1430);
    //dtoppp.resonances.push_back(omega_782); 
    //dtoppp.resonances.push_back(rho_770); 
    //dtoppp.resonances.push_back(rho_1450);
    //dtoppp.resonances.push_back(rho_1700);
    dtoppp.resonances.push_back(ks_892);
    dtoppp.resonances.push_back(phi_1020);
    dtoppp.resonances.push_back(ks_1410);
    dtoppp.resonances.push_back(phi_1680);
    dtoppp.resonances.push_back(f2_1270);
    dtoppp.resonances.push_back(a2_1320);
    dtoppp.resonances.push_back(f2_1525);
    dtoppp.resonances.push_back(k2s_1430);
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

