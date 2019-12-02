// input.h
#ifndef INPUT_H
#define INPUT_H
 
namespace GooFit{

#define torad(x)(x*M_PI/180)
//Root Style
//Globals

const double pi_MASS  = 0.13957018; //GEV
const double k_MASS   = 0.493677;
const double DS_MASS  = 1.96834; 
const double D_MASS   = 1.86965;

const double d1_MASS  = pi_MASS;  //daughter 1 mass
const double d2_MASS  = pi_MASS;  //daughter 2 mass	
const double d3_MASS  = pi_MASS;  //daughter 3 mass

fptype s12_min = POW2(d1_MASS  + d2_MASS);
fptype s12_max = POW2(D_MASS   - d2_MASS);
fptype s13_min = POW2(d1_MASS  + d3_MASS);
fptype s13_max = POW2(D_MASS   - d3_MASS);
fptype s23_min = POW2(d2_MASS  + d3_MASS);
fptype s23_max = POW2(D_MASS   - d3_MASS);

Observable s12("s12",s12_min,s12_max); 
Observable s13("s13",s13_min,s13_max);
Observable s23("s23",s23_min,s23_max);
EventNumber eventNumber("eventNumber");

GooPdf *makeDstar_veto();
DalitzPlotPdf* makesignalpdf(GooPdf* eff);

}

#endif
