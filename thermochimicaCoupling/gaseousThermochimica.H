
//Ideal Gas Phase
GetSolnPhaseMol(phaseName5, &molesPhase_gas_ideal, &info);
SolnPhaseMol_gas_ideal.push_back(molesPhase_gas_ideal);
idealGasPhase[i] = molesPhase_gas_ideal;
//Info << "Ideal Gas Phase Fraction " << idealGasPhase[i] << "[mol]" << nl;

//F2 in Ideal Gas
GetOutputMolSpeciesPhase(phaseName5, speciesName4, &moleFraction_F2, &info);
molFrac_F2.push_back(moleFraction_F2);
F2_gas[i] = (idealGasPhase[i]*moleFraction_F2*8.31446261815324*phases[0].thermo().T()[i]/phases[0].thermoRef().p()[i])/cv[i];

//F2_gas[i] = molesPhase_gas_ideal*moleFraction_F2;
