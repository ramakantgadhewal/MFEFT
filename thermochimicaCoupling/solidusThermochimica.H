//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//solid1 (LiF-NaF-KF-CsF)
GetSolnPhaseMol(phaseName1, &molesPhase_solid1, &info);
SolnPhaseMol_solid1.push_back(molesPhase_solid1);
//printf("%s fraction as a stable solution phase: %e\n", phaseName1, molesPhase_solid1);
solid1[i] = molesPhase_solid1;

//LiF
GetOutputMolSpeciesPhase(phaseName1, speciesName1, &moleFraction_LiF, &info);
molFrac_LiF.push_back(moleFraction_LiF);
LiF[i] = moleFraction_LiF*solid1[i];

//NaF
GetOutputMolSpeciesPhase(phaseName1, speciesName5, &moleFraction_NaF, &info);
molFrac_NaF.push_back(moleFraction_NaF);
NaF[i] = moleFraction_NaF*solid1[i];

//KF
GetOutputMolSpeciesPhase(phaseName1, speciesName2, &moleFraction_KF, &info);
molFrac_KF.push_back(moleFraction_KF);
KF[i] = moleFraction_KF*solid1[i];

//CsF
GetOutputMolSpeciesPhase(phaseName1, speciesName3, &moleFraction_CsF, &info);
molFrac_CsF.push_back(moleFraction_CsF);
CsF[i] = moleFraction_CsF*solid1[i];

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//NaF-KF
GetSolnPhaseMol(phaseName6, &molesPhase_NaF_KF, &info);
SolnPhaseMol_NaF_KF.push_back(molesPhase_NaF_KF);
NaF_KF[i] = molesPhase_NaF_KF;

//NaF
GetOutputMolSpeciesPhase(phaseName6, speciesName5, &moleFraction_NaF, &info);
molFrac_NaF.push_back(moleFraction_NaF);
NaF[i] = NaF[i] + moleFraction_NaF*NaF_KF[i];

//KF
GetOutputMolSpeciesPhase(phaseName6, speciesName2, &moleFraction_KF, &info);
molFrac_KF.push_back(moleFraction_KF);
KF[i] = KF[i] + moleFraction_KF*NaF_KF[i];

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//LiF-NaF
GetSolnPhaseMol(phaseName2, &molesPhase_LiF_NaF, &info);
SolnPhaseMol_LiF_NaF.push_back(molesPhase_LiF_NaF);
LiF_NaF[i] = molesPhase_LiF_NaF;

//NaF
GetOutputMolSpeciesPhase(phaseName2, speciesName5, &moleFraction_NaF, &info);
molFrac_NaF.push_back(moleFraction_NaF);
NaF[i] = NaF[i] + moleFraction_NaF*LiF_NaF[i];

//LiF
GetOutputMolSpeciesPhase(phaseName2, speciesName1, &moleFraction_LiF, &info);
molFrac_LiF.push_back(moleFraction_LiF);
LiF[i] = LiF[i] + moleFraction_LiF*LiF_NaF[i];

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Pure Condensed Phases

//LaK3F6_solid
GetPureConPhaseMol(phaseName7, &molesPhase_LaK3F6, &info);
SolnPhaseMol_LaK3F6.push_back(molesPhase_LaK3F6);
LaK3F6[i] = molesPhase_LaK3F6;

//CeF3_solid
GetPureConPhaseMol(phaseName8, &molesPhase_CeF3, &info);
SolnPhaseMol_CeF3.push_back(molesPhase_CeF3);
CeF3[i] = molesPhase_CeF3;

//CsF2Li_solids
GetPureConPhaseMol(phaseName9, &molesPhase_CsF2Li_solids, &info);
SolnPhaseMol_CsF2Li_solids.push_back(molesPhase_CsF2Li_solids);
CsF2Li_solids[i] = molesPhase_CsF2Li_solids;

//FLi_solids
GetPureConPhaseMol(phaseName11, &molesPhase_FLi_solids, &info);
SolnPhaseMol_FLi_solids.push_back(molesPhase_FLi_solids);
FLi_solids[i] = molesPhase_FLi_solids;

//CsF_solid
GetPureConPhaseMol(phaseName13, &molesPhase_CsF_solids, &info);
SolnPhaseMol_CsF_solids.push_back(molesPhase_CsF_solids);
CsF_Solids[i] = (molesPhase_CsF_solids*0.151903/4120)/cv[i];


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Total Number of Solid molesPhase_solid1

FLiNaK_Solid[i] = solid1[i] + NaF_KF[i] + LiF_NaF[i];

//Info<< "FLiNaK Solid Fraction: " << FLiNaK_Solid << nl;

PureCondPhase[i] = CsF2Li_solids[i] + FLi_solids[i] + CeF3[i] + LaK3F6[i];
