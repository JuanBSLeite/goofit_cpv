// input.h
#ifndef INPUT_H
#define INPUT_H
 
namespace GooFit{

#define torad(x)(x*M_PI/180)
//Root Style
//Globals

const double pi_MASS  = 0.13957018; //GEV
const double D_MASS   = 1.96834; //GEV

const double d1_MASS  = pi_MASS;  //daughter 1 mass
const double d2_MASS  = pi_MASS;
const double d3_MASS  = pi_MASS;

fptype s12_min = POW2(d1_MASS  + d2_MASS);
fptype s12_max = POW2(D_MASS   - d2_MASS);
fptype s13_min = POW2(d1_MASS  + d3_MASS);
fptype s13_max = POW2(D_MASS   - d3_MASS);
fptype s23_min = POW2(d2_MASS  + d3_MASS);
fptype s23_max = POW2(D_MASS   - d3_MASS);

Observable s12("s12",s12_min,s12_max); //s12^{2}
Observable s13("s13",s13_min,s13_max);
Observable s23("s23",s23_min,s23_max);
EventNumber eventNumber("eventNumber");

DalitzPlotPdf* makesignalpdf(GooPdf* eff);

}

#endif
