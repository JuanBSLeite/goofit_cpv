// input.cpp
#include <vector>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>

#include "input.h"


namespace GooFit{

DalitzPlotPdf* makesignalpdf(GooPdf* eff = 0){

    DecayInfo3 dtoppp;
    dtoppp.motherMass   = D_MASS;
    dtoppp.daug1Mass    = d1_MASS;
    dtoppp.daug2Mass    = d2_MASS;
    dtoppp.daug3Mass    = d3_MASS;
    dtoppp.meson_radius = 1.5;

    //Mass and width

    double f0_980_MASS    = 0.990;
    double f0_980_GPP     = 0.02;
    double f0_980_GKK     = 4.5*0.02;
    double f0_980_WIDTH   = 0.04;
    double f0_980_amp     = 1.0;
    double f0_980_img     = 0.0;
  
    
    double f0_1370_MASS  = 1.370;
    double f0_1370_WIDTH = .3;
    double f0_1370_amp   = 0.75*cos(torad(198));
    double f0_1370_img = 0.75*sin(torad(198));

    double f0_1500_MASS  = 1.505;
    double f0_1500_WIDTH = .109;
    double f0_1500_amp   = 1.;
    double f0_1500_img = 0.;

    double omega_MASS   = 0.78265;
    double omega_WIDTH  = 0.00849;
    double omega_amp    = 1.;
    double omega_img  =   0.;

    double rho770_MASS   = .77549;
    double rho770_WIDTH  = .1491;
    //E791
    //double rho770_amp    = 0.32*cos(torad(109));
    //double rho770_img  = 0.32*sin(torad(109));
    //Babar
    double rho770_amp    = 0.19*cos(1.1);
    double rho770_img  = 0.19*sin(1.1);
  


    double rho1450_MASS   = 1.465 + 3*0.025;
    double rho1450_WIDTH  = 0.54;//0.4;
    //E791
    //double rho1450_amp    = 0.28*cos(torad(162));
    //double rho1450_img  = 0.28*sin(torad(162));
    //Babar
    double rho1450_amp    = 1.2*cos(4.1);
    double rho1450_img  = 1.2*sin(4.1);
 

    
    
    double f2_1270_MASS     = 1.2751;
    double f2_1270_WIDTH    = 0.1851;
    //E791
    //double f2_1270_amp      = 0.59*cos(torad(133));
    //double f2_1270_img    = 0.59*sin(torad(133));
    //Babar
    double f2_1270_amp      = 1.;
    double f2_1270_img    = 0.;

    double f2_1525_MASS     = 1.525;
    double f2_1525_WIDTH    = 0.073;
    double f2_1525_amp      = 1.;
    double f2_1525_img    = 0.;
    

    //omega(782)
    Variable v_omega_Mass("omega_MASS",omega_MASS,0.01,0,0);
    Variable v_omega_Width("omega_WIDTH",omega_WIDTH,0.01,0,0);
    Variable v_omega_real("omega_real",omega_amp, 0.01,0,0);
    Variable v_omega_img("omega_img",omega_img, 0.01,0,0);

    //rho(770)
    Variable v_rho770_Mass("rho770_MASS",rho770_MASS,0.01,0,0);
    Variable v_rho770_Width("rho770_WIDTH",rho770_WIDTH,0.01,0,0);
    Variable v_rho770_real("rho770_real",rho770_amp, 0.01,0,0);
    Variable v_rho770_img("rho770_img",rho770_img, 0.01,0,0);
    
    //rho(1450)
    Variable v_rho1450_Mass("rho1450_MASS",rho1450_MASS,0.01,0,0);
    Variable v_rho1450_Width("rho1450_WIDTH",rho1450_WIDTH,0.01,0,0);
    Variable v_rho1450_real("rho1450_real",rho1450_amp, 0.01,0,0);
    Variable v_rho1450_img("rho1450_img",rho1450_img, 0.01,0,0);

	//f2(1270)
    Variable v_f2_1270_Mass("f2_1270_MASS",f2_1270_MASS,0.01,0,0);
    Variable v_f2_1270_Width("f2_1270_WIDTH",f2_1270_WIDTH,0.01,0,0);
    Variable v_f2_1270_real("f2_1270_real",f2_1270_amp, 0.01,0,0);
    Variable v_f2_1270_img("f2_1270_img",f2_1270_img, 0.01,0,0);

    //f2(1525)
    Variable v_f2_1525_Mass("f2_1525_MASS",f2_1525_MASS,0.01,0,0);
    Variable v_f2_1525_Width("f2_1525_WIDTH",f2_1525_WIDTH,0.01,0,0);
    Variable v_f2_1525_real("f2_1525_real",f2_1525_amp, 0.01,0,0);
    Variable v_f2_1525_img("f2_1525_img",f2_1525_img, 0.01,0,0);

    //f0(980)
    Variable v_f0_980_Mass("f0_980_MASS",f0_980_MASS,.02,0,0);
    Variable v_f0_980_GPP("f0_980_GPP",f0_980_GPP,0.01,0.,0);
    Variable v_f0_980_GKK("f0_980_GKK",f0_980_GKK,0.01,0.,0.);
    Variable v_f0_980_Width("f0_980_WIDTH",f0_980_WIDTH,0.01,0,0);
    Variable v_f0_980_real("f0_980_real",f0_980_amp, 0.01,0,0);
    Variable v_f0_980_img("f0_980_img",f0_980_img, 0.01,0,0);

    v_f0_980_real.setFixed(true);
    v_f0_980_img.setFixed(true);
    //v_f2_1270_real.setFixed(true);
    //v_f2_1270_img.setFixed(true);

    
    //f0(1370)
    Variable v_f0_1370_Mass("f0_1370_MASS",f0_1370_MASS,0.025,0,0);
    Variable v_f0_1370_Width("f0_1370_Width",f0_1370_WIDTH,0.035,0,0);
    Variable v_f0_1370_real("f0_1370_real",f0_1370_amp, 0.01,0,0);
    Variable v_f0_1370_img("f0_1370_img",f0_1370_img, 0.01, 0,0);

    
    //f0(1500)
    Variable v_f0_1500_Mass("f0_1500_MASS",f0_1500_MASS,0.06,0,0);
    Variable v_f0_1500_Width("f0_1500_Width",f0_1500_WIDTH,0.07,0,0);
    Variable v_f0_1500_real("f0_1500_real",f0_1500_amp, 0.01,0,0);
    Variable v_f0_1500_img("f0_1500_img",f0_1500_img, 0.01, 0,0);

    //NR
    Variable nonr_real("nonr_real",0.09*cos(torad(181)), 0.01,0,0);
    Variable nonr_imag("nonr_imag",0.09*sin(torad(181)), 0.01,0,0);
    //Variable nonr_real("nonr_real",1., 0.01,0,0);
     // Variable nonr_imag("nonr_imag",0., 0.01,0,0);


    Variable be_real("be_real",1., 0.01,0,0);
    Variable be_imag("be_imag", 0., 0.01,0,0);
    Variable be_coef("be_coef", 1.9,0.01,0,0);
    //Masses and Widths fixed

    v_omega_Mass.setFixed(true);
    v_omega_Width.setFixed(true);
   
    v_rho770_Mass.setFixed(true);
    v_rho770_Width.setFixed(true);
    
    v_rho1450_Mass.setFixed(true);
    v_rho1450_Width.setFixed(true);
    
    v_f2_1270_Mass.setFixed(true);
    v_f2_1270_Width.setFixed(true);

    v_f2_1525_Mass.setFixed(true);
    v_f2_1525_Width.setFixed(true);

    //v_f0_980_Mass.setFixed(true);
    //v_f0_980_GKK.setFixed(true);
    //v_f0_980_GPP.setFixed(true);
    v_f0_980_Width.setFixed(true);

    v_f0_1370_Mass.setFixed(true);
    v_f0_1370_Width.setFixed(true);

    v_f0_1500_Mass.setFixed(true);
    v_f0_1500_Width.setFixed(true);

    be_coef.setFixed(true);
        
   //setting resonances
    ResonancePdf* sigma_12 = new Resonances::POLE("sigma",Variable("v_sigma_real",1.),Variable("v_sigma_img",0.),Variable("v_sigma__pole_real",0.47),Variable("v_sigma_pole_img",0.22),0,PAIR_12,true);
    
    ResonancePdf* omega_12 = new Resonances::RBW("omega",v_omega_real,v_omega_img,v_omega_Mass,v_omega_Width,1,PAIR_12,true);
    
    ResonancePdf* rho770_12 = new Resonances::RBW("rho770",v_rho770_real,v_rho770_img,v_rho770_Mass,v_rho770_Width,1,PAIR_12,true);
    
    ResonancePdf* rho1450_12 = new Resonances::RBW("rho1450",v_rho1450_real,v_rho1450_img,v_rho1450_Mass,v_rho1450_Width,1,PAIR_12,true);
    
    ResonancePdf* f2_1270_12 = new Resonances::RBW("f2",v_f2_1270_real,v_f2_1270_img,v_f2_1270_Mass,v_f2_1270_Width,2,PAIR_12,true);

    ResonancePdf* f2_1525_12 = new Resonances::RBW("f2",v_f2_1525_real,v_f2_1525_img,v_f2_1525_Mass,v_f2_1525_Width,2,PAIR_12,true);
  
    ResonancePdf* f0_980_12 = new Resonances::FLATTE("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_GPP,v_f0_980_GKK,PAIR_12,true);

    //ResonancePdf* f0_980_12 = new Resonances::RBW("f0_980",v_f0_980_real,v_f0_980_img,v_f0_980_Mass,v_f0_980_Width,(unsigned int)0,PAIR_12,true);
	
    ResonancePdf* f0_1370_12 = new Resonances::RBW("f0_1370_12",v_f0_1370_real,v_f0_1370_img,v_f0_1370_Mass,v_f0_1370_Width,(unsigned int)0,PAIR_12,true);

    ResonancePdf* f0_1500_12 = new Resonances::RBW("f0_1500_12",v_f0_1500_real,v_f0_1500_img,v_f0_1500_Mass,v_f0_1500_Width,(unsigned int)0,PAIR_12,true);  

    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_real, nonr_imag);

    ResonancePdf *be   = new Resonances::BoseEinstein("be",be_real,be_imag,be_coef);

   
    

  

    //Pushing Resonances 
    //dtoppp.resonances.push_back(sigma_12);
    dtoppp.resonances.push_back(omega_12); 
    dtoppp.resonances.push_back(rho770_12); 
    dtoppp.resonances.push_back(rho1450_12);
    //dtoppp.resonances.push_back(f2_1270_12);
    //dtoppp.resonances.push_back(f2_1525_12);
    //dtoppp.resonances.push_back(f0_980_12);
    //dtoppp.resonances.push_back(f0_1500_12);
    //dtoppp.resonances.push_back(f0_1370_12);
    //dtoppp.resonances.push_back(nonr);
    //dtoppp.resonances.push_back(be);
 



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

