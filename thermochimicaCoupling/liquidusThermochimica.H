//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//SALTsoln
GetSolnPhaseMol(phaseName4, &molesPhase_SALTsoln, &info);
SolnPhaseMol_SALTsoln.push_back(molesPhase_SALTsoln);
//printf("%s SALTsoln Phase Fraction: %e\n", phaseName4, molesPhase_SALTsoln);
FLiNaK_Liquid[i] = molesPhase_SALTsoln/3;
//Info << "Total Liquid FLiNaK= " << FLiNaK_Liquid[i] << " [mol]" << nl;

//Cs_liquid
GetPureConPhaseMol(phaseName12, &molesPhase_Cs_liquid, &info);
SolnPhaseMol_Cs_liquid.push_back(molesPhase_Cs_liquid);
Cs_Liquid[i] = (molesPhase_Cs_liquid*0.132905/1873)/cv[i];;
