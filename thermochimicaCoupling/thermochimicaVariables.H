std::vector<double> alphaSolid;

//Liquid
std::vector<double> SolnPhaseMol_SALTsoln;
std::vector<double> SolnPhaseMol_K_liquid;

//Solid
std::vector<double> SolnPhaseMol_solid1;
std::vector<double> SolnPhaseMol_LiF_NaF;
std::vector<double> SolnPhaseMol_LiF_KF;
std::vector<double> SolnPhaseMol_NaF_KF;
std::vector<double> SolnPhaseMol_NaF_solids;

//Pure Condensed Vectors
std::vector<double> SolnPhaseMol_LaK3F6;
std::vector<double> SolnPhaseMol_CeF3;
std::vector<double> SolnPhaseMol_CsF2Li_solids;
std::vector<double> SolnPhaseMol_FLi_solids;
std::vector<double> SolnPhaseMol_Cs_liquid;
std::vector<double> SolnPhaseMol_CsF_solids;

//Gases
std::vector<double> SolnPhaseMol_gas_ideal;

//Phase Components
std::vector<double> molFrac_KF;
std::vector<double> molFrac_K;
std::vector<double> molFrac_NaF;
std::vector<double> molFrac_LiF;
std::vector<double> molFrac_CsF;
std::vector<double> molFrac_Cs;
std::vector<double> molFrac_FLi;
std::vector<double> molFrac_F2;

char phaseName1[25] = "solid1";
char phaseName2[25] = "LiF-NaF";
//char phaseName3[25] = "";
char phaseName4[25] = "SALTsoln";
char phaseName5[25] = "gas_ideal";
char phaseName6[25] = "NaF-KF";
char phaseName7[25] = "LaK3F6_solid";
char phaseName8[25] = "CeF3_solid";
char phaseName9[25] = "CsF2Li_solids";
//char phaseName10[25] = "";
char phaseName11[25] = "FLi_solids";

char phaseName12[25] = "Cs_liquid";
char phaseName13[25] = "CsF_solids";

double molesPhase_solid1;
double molesPhase_LiF_NaF;
double molesPhase_NaF_KF;
double molesPhase_LaK3F6;
double molesPhase_CeF3;
double molesPhase_CsF2Li_solids;
double molesPhase_FLi_solids;
double molesPhase_SALTsoln;
//double molesPhase_K_liquid;
double molesPhase_CsF_solids;
double molesPhase_Cs_liquid;
double molesPhase_gas_ideal;

char speciesName1[8] = "LiF";
char speciesName2[8] = "KF";
char speciesName3[8] = "CsF";
char speciesName5[8] = "NaF";
char speciesName4[8] = "F2";

double moleFraction_NaF;
double moleFraction_LiF;
double moleFraction_CsF;
double moleFraction_KF;
double moleFraction_F2;
