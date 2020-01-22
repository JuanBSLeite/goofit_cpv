// input.h
#ifndef INPUT_H
#define INPUT_H
 
namespace GooFit{

#define torad(x)(x*M_PI/180)
#define real(x,y)(x*cos(y) )
#define img(x,y)(x*sin(y))
const double pi_MASS  = 0.13957018; //GEV
const double k_MASS   = 0.493677;
const double DS_MASS  = 1.96849; 
const double D_MASS   = 1.86962;

//Set the final state
const double Mother_MASS = D_MASS;
const double d1_MASS  = k_MASS;  //daughter 1 mass
const double d2_MASS  = k_MASS;  //daughter 2 mass	
const double d3_MASS  = pi_MASS;  //daughter 3 mass

bool symdp	= false;

//if true include efficiency and(or) background
bool effOn      = false;
bool bkgOn	= true;

//data sample for fitting
std::string DataFile = "/data1000/lhcb_charm/data/2016/DKKP/DKKP_Up_2016_21_27.root";
std::string TreeName = "DecayTree";
std::string s12Name  = "s12_KK_DTF";//branch name
std::string s13Name  = "s13_Kpi_DTF";//branch name
//background and eff paths and histo names
std::string bkg_file = "/home/juan/juan/work/kkpi/BKG_wl.root";
std::string eff_file = "../../dados/eff_16.root";
std::string bkg_name = "h0";
std::string eff_name = "h";

fptype s12_min = POW2(d1_MASS  + d2_MASS);
fptype s12_max = POW2(Mother_MASS   - d3_MASS);
fptype s13_min = POW2(d1_MASS  + d3_MASS);
fptype s13_max = POW2(Mother_MASS    - d2_MASS);
fptype s23_min = POW2(d2_MASS  + d3_MASS);
fptype s23_max = POW2(Mother_MASS   - d1_MASS);

Observable s12("s12",s12_min,s12_max); 
Observable s13("s13",s13_min,s13_max);
Observable s23("s23",s23_min,s23_max);
EventNumber eventNumber("eventNumber");

GooPdf *makeDstar_veto();
DalitzPlotPdf* makesignalpdf(GooPdf* eff);

}

#endif
