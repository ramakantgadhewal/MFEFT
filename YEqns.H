{
    autoPtr<phaseSystem::specieTransferTable>
        specieTransferPtr(fluid.specieTransfer());

    phaseSystem::specieTransferTable&
        specieTransfer(specieTransferPtr());

    fluid.correctReactions(); //if no chemical reaction, it jumps

    forAll(fluid.multiComponentPhases(), multiComponentPhasei) //multicomponent because 1 phase can be formed by several species
    {
        phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i) //if gas.AIR gas.CO gas.CO2 gas.H i=4 for multicomponentphasei=1(gas)
        {
            fvScalarMatrix YiEqn //So I suppose that the mass fraction of a species in a cell, is the value of Y times the value of the volume fraction in the cell (the value of alpha.phase). Is that the case?
            (
                phase.YiEqn(Y[i])  //all i species from current multiComponentPhase[multiComponentPhasei] multiComponentPhasei=1,2 and so on
             ==
               *specieTransfer[Y[i].name()]
              + fvModels.source(alpha, rho, Y[i])
            );

            YiEqn.relax();
            YiEqn.solve("Yi");
        }
    }

    fluid.correctSpecies(); //end all species fraction Y, updating CO CO2 H AIR for phase gas
}


// If a system is Fluid + Gas, there is two phases.
//If a system is Fluid(species: a, b and c) + Gas(d and e). Here it is calculating (Y_a + Y_b + Y_c = 1(fluid)) (Y_d + Y_e = 1(Gas)
//where Fluid is multiComponentPhase[1] and Gas is multiComponentPhase[2]
