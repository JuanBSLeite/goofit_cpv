// input.cpp
#include <vector>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>
#include "input.h"


namespace GooFit{


double Signal_Purity = 0.8;

double Meio= 0.5;

double VarAmp_sigma=0.;
double VarPha_sigma=0.;
double VarAmp_f0_1370=0.;
double VarPha_f0_1370=0.;
double VarAmp_rho1450=0.;
double VarPha_rho1450=0.;

//Valores de amplitude e fase (graus)			Valores originais do modelo isobarico da Fernanda:
double sigma_amp = ((100.0 + VarAmp_sigma )/100.)*2.88974237246;		//	2.88974237246
double sigma_pha = ( VarPha_sigma *1.0)+21.3324464125;				//	21.3324464125

double f0_980_amp =  2.07391338443;		//	2.07391338443
double f0_980_pha = 17.3221794525;		//	17.3221794525

double f0_1500_amp =  0.942826506775;		//	0.942826506775
double f0_1500_pha = -11.904153174;		//	-11.904153174

double f0_1370_amp = ((100.0 + VarAmp_f0_1370 )/100.)*2.27256585778;		//	2.27256585778
double f0_1370_pha = ( VarPha_f0_1370 *1.0)+1.48931393812;			//	1.48931393812

double romix_amp = 1.0;				//	1.0
double romix_pha = 0.0;				//	0.0

double omega_amp = 0.0140304708354;		//	0.0140304708354
double omega_pha = -75.0347878841;		//	-75.0347878841

double f2_1270_amp = 1.8950683295;		//	1.8950683295
double f2_1270_pha = -97.1867270077;		//	-97.1867270077

double rho1450_amp = ((100.0 + VarAmp_rho1450 )/100.)*0.830919761702;		//	0.830919761702
double rho1450_pha = ( VarPha_rho1450 *1.0)+101.068997745;			//	101.068997745



//include veto
GooPdf *makeDstar_veto() {
   								//Veto range from a to b
    VetoInfo Fiducial_veto12(Variable("Fiducial_veto12_min", 0.), Variable("Fiducial_veto12_max", s12.getLowerLimit()), PAIR_12);
    VetoInfo Fiducial_veto13(Variable("DFiducial_veto13_min", 0.), Variable("Fiducial_veto13_max", s13.getLowerLimit()), PAIR_13);
    
    std::vector<VetoInfo> vetos;
    //vetos.push_back(Fiducial_veto12);
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
							Variable("sigma_real",sigma_amp * cos(2.*3.14159* sigma_pha * (1./360.))),
							Variable("sigma_img",sigma_amp * sin(2.*3.14159* sigma_pha * (1./360.))),
							Variable("sigma_pole_real",0.47,.001,0,0),
							Variable("sigma_pole_img",0.22,.001,0,0),0,PAIR_12,symdp);
    
	ResonancePdf* ks_800         = new Resonances::POLE("ks_800",
                                                        Variable("ks_800_real",1.),
                                                        Variable("ks_800_img",0.),
                                                        Variable("ks_800_pole_real",0.680,.001,0,0),
                                                        Variable("ks_800_pole_img",0.300,.001,0,0),0,PAIR_12,symdp);

	ResonancePdf* f0_980	= new Resonances::FLATTE("f0_980",
							Variable("f0_980_real",f0_980_amp * cos(2.*3.14159* f0_980_pha * (1./360.))),
							Variable("f0_980_img",f0_980_amp * sin(2.*3.14159* f0_980_pha * (1./360.))),
							Variable("f0_980_mass",0.965),
							Variable("f0_980_gpp",0.165),
							Variable("f0_980_gkk",0.694),PAIR_12,symdp);

	ResonancePdf* f0_980_BW	= new Resonances::RBW("f0_980",
							Variable("f0_980_real",1.,.001,0,0),
							Variable("f0_980_img",0.,.001,0,0),
							Variable("f0_980_mass",0.965),
							Variable("f0_980_width",0.04),int(0),PAIR_12,symdp);

	ResonancePdf* f0_1370	= new Resonances::RBW("f0_1370",
							Variable("f0_1370_real",f0_1370_amp * cos(2.*3.14159* f0_1370_pha * (1./360.)),.001,0,0),
							Variable("f0_1370_img",f0_1370_amp * sin(2.*3.14159* f0_1370_pha * (1./360.)),.001,0,0),
							Variable("f0_1370_mass",1.4154),
							Variable("f0_1370_width",0.304),int(0),PAIR_12,symdp);
	
	ResonancePdf* f0_1500	= new Resonances::RBW("f0_1500",
							Variable("f0_1500_real",f0_1500_amp * cos(2.*3.14159* f0_1500_pha * (1./360.)),.001,0,0),
							Variable("f0_1500_img",f0_1500_amp * sin(2.*3.14159* f0_1500_pha * (1./360.)),.001,0,0),
							Variable("f0_1500_mass",1.505),
							Variable("f0_1500_width",0.109),int(0),PAIR_12,symdp);
	
	ResonancePdf* omega_782	= new Resonances::RBW("omega_782",
							Variable("omega_782_real",omega_amp * cos(2.*3.14159* omega_pha * (1./360.)),.001,0,0),
							Variable("omega_782_img",omega_amp * sin(2.*3.14159* omega_pha * (1./360.)),.001,0,0),
							Variable("omega_782_mass",0.782),
							Variable("omega_782_width",0.0085),int(1),PAIR_12,symdp);
	
    	ResonancePdf* rho_770	= new Resonances::RBW("rho_770",
							Variable("rho_770_real",romix_amp * cos(2.*3.14159* romix_pha * (1./360.)),.001,0,0),
							Variable("rho_770_img",romix_amp * sin(2.*3.14159* romix_pha * (1./360.)),.001,0,0),
							Variable("rho_770_mass",.7755),
							Variable("rho_770_width",0.149),int(1),PAIR_12,symdp);
	
	ResonancePdf* rho_1450	= new Resonances::RBW("rho_1450",
							Variable("rho_1450_real",rho1450_amp * cos(2.*3.14159* rho1450_pha * (1./360.)),.001,0,0),
							Variable("rho_1450_img",rho1450_amp * sin(2.*3.14159* rho1450_pha * (1./360.)),.001,0,0),
							Variable("rho_1450_mass",1.465),
							Variable("rho_1450_width",0.4),int(1),PAIR_12,symdp);

        ResonancePdf* ks_892	  = new Resonances::RBW("ks_892",
                                                        Variable("ks_892_real",1.,.001,0,0),
                                                        Variable("ks_892_img",0.,.001,0,0),
                                                        Variable("ks_892_mass",0.8917),
                                                        Variable("ks_892_width",0.05),int(1),PAIR_12,symdp);

        ResonancePdf* ks_1430      = new Resonances::RBW("ks_1430",
                                                        Variable("ks_1430_real",1.,.001,0,0),
                                                        Variable("ks_1430_img",0.,.001,0,0),
                                                        Variable("ks_1430_mass",1.425),
                                                        Variable("ks_1430_width",0.270),int(1),PAIR_12,symdp);

        ResonancePdf* phi_1020     = new Resonances::RBW("phi_1020",
                                                        Variable("phi_1020_real",1.,.001,0,0),
                                                        Variable("phi_1020_img",0.,.001,0,0),
                                                        Variable("phi_1020_mass",1.019461),
                                                        Variable("phi_1020_width",0.0048),int(1),PAIR_12,symdp);

	ResonancePdf* f2_1270	= new Resonances::RBW("f2_1270",
							Variable("f2_1270_real",f2_1270_amp * cos(2.*3.14159* f2_1270_pha * (1./360.)),.001,0,0),
							Variable("f2_1270_img",f2_1270_amp * sin(2.*3.14159* f2_1270_pha * (1./360.)),.001,0,0),
							Variable("f2_1270_mass",1.2751),
							Variable("f2_1270_width",0.185),int(2),PAIR_12,symdp);
	
	ResonancePdf* f2_1525	= new Resonances::RBW("f2_1525",
							Variable("f2_1525_real",1.,.001,0,0),
							Variable("f2_1525_img",0.,.001,0,0),
							Variable("f2_1525_mass",1.525),
							Variable("f2_1525_width",0.073),int(2),PAIR_12,symdp);
	
	ResonancePdf* k2s_1430   = new Resonances::RBW("k2s_1430",
                                                        Variable("k2s_1430_real",1.,.001,0,0),
                                                        Variable("k2s_1430_img",0.,.001,0,0),
                                                        Variable("k2s_1430_mass",1.4256),
                                                        Variable("k2s_1430_width",0.0985),int(2),PAIR_12,symdp);	
  
    	ResonancePdf *nonr 		= new Resonances::NonRes("nonr",
							Variable("nonr_real",1.,.001,0,0), 
							Variable("nonr_img",0,.001,0,0));


//Push Resonances of your model 
    dtoppp.resonances.push_back(sigma);
    //dtoppp.resonances.push_back(ks_800);
    dtoppp.resonances.push_back(f0_980);
    //dtoppp.resonances.push_back(f0_980_BW);
    dtoppp.resonances.push_back(f0_1500);
    dtoppp.resonances.push_back(f0_1370);
    dtoppp.resonances.push_back(omega_782); 
    dtoppp.resonances.push_back(rho_770); 
    dtoppp.resonances.push_back(rho_1450);
    //dtoppp.resonances.push_back(ks_892);
    //dtoppp.resonances.push_back(ks_1430);
    //dtoppp.resonances.push_back(phi_1020);
    dtoppp.resonances.push_back(f2_1270);
    //dtoppp.resonances.push_back(f2_1525);
    //dtoppp.resonances.push_back(k2s_1430);
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

