//Calculate Phase Fractions ( && = and  ; || = or)

  forAll(phases, phasei) //for all phases being named from 1 to i(phasei)
  {

    phaseModel& phase = phases[phasei];

      if(phasei == 2) //==2 because here there is just three phases (0,1 and 2) and using this condition I avoid 2 unnecessary calc loops
      {
        forAll(phase, i)
        {
              //Acceleration Subroutine
              //This condition also limits the minimum concentration for a cell be reactive, which is less than calculating every cells
              //This also produce a different amount of solid, which also produce a different flow dynamics and transient solid formation
              if( (phases[0][i] > 0.99999) ||                                                           //PURE GAS
                                              (phases[1][i] > 0.99999) ||                               //PURE LIQUID
                 ((phases[0][i] > 0.00001) && (phases[1][i] < 0.00001) && (phases[2][i] > 0.00001)) || //MIXTURE WITHOUT LIQUID
                 ((phases[0][i] < 0.00001) && (phases[1][i] > 0.00001) && (phases[2][i] > 0.00001)))   //MIXTURE WITHOUT GAS
              {
                Info << "Cell: " << i << " IS NOT REACTIVE" << nl;
                Info<< nl;
              }
              else
              {
                Info << "Cell: " << i << " IS REACTIVE, UPDATING COEFFICIENTS" << nl;
                Info << "#: " << i << "/" << nCells << nl;
                Info<< "Time Step = " << runTime.value() << " s" << nl;
                Info<< nl;

                #include "thermochimicaCoupling/thermochimicaCalc.H"

                Info << "Gas TC    : " << F2_gas[i] << nl;
                Info << "Liquid  TC: " << Cs_Liquid[i] << nl;
                Info << "Solid TC  : " << CsF_Solids[i] << nl;

                    //Sometimes the very fist TC calculation return stores 0 values for all vectors,
                    //random numbers are allocated to avoid calc error
                    //Can it be related to the random value a vector stores before get it used?
                    //Do I need to use like F2_gas[i]{} or F2_gas[i]{0} to avoid id store junk value?
                    if((F2_gas[i] == Cs_Liquid[i]) && (Cs_Liquid[i] == CsF_Solids[i]) && (CsF_Solids[i] == 0))
                    {
                        F2_gas[i] = 0;
                        Cs_Liquid[i] = 0.994955;
                        CsF_Solids[i] = 6.19164e-06;

                        if((F2_gas[i] == Cs_Liquid[i]) && (Cs_Liquid[i] == CsF_Solids[i]) && (CsF_Solids[i] == 0))
                        {
                          return 0 ;
                        }
                    }

                    if(F2_gas[i]<0)
                    {
                      F2_gas[i] = double(0);
                    }
                    if(Cs_Liquid[i]<0)
                    {
                      Cs_Liquid[i] = double(0);
                    }
                    if(CsF_Solids[i]<0)
                    {
                      CsF_Solids[i] = double(0);
                    }

                Info << nl;
                Info << "Alpha TC: " << F2_gas[i] + Cs_Liquid[i] + CsF_Solids[i] << nl;
                Info << nl;

                    if(F2_gas[i] == CsF_Solids[i] && CsF_Solids[i] == Cs_Liquid[i])
                    {
                      Info << "FATAL ERROR: THERMOCHIMICA RETURNED 0" << nl;
                      Info << "CALCULATION ENDED" << nl;
                      return 0;
                    }

                //Old Phase Fractions before Thermochimica Calculation
                Info << "Gas OF    : " << phases[0][i] << nl;
                Info << "Liquid  OF: " << phases[1][i] << nl;
                Info << "Solid OF  : " << phases[2][i] << nl;

                Info<< nl;
                Info << "Alpha OF: " << phases[0][i]+phases[1][i]+phases[2][i] << nl;

                //Phase Fractions Update from TC Results

                phases[1][i] = Cs_Liquid[i];
                phases[2][i] = CsF_Solids[i] + phases[2][i];
                phases[0][i] = double(1) - phases[1][i] - phases[2][i]; //F2_gas[i];

                Info<< nl;
                Info << "Gas t    : " << phases[0][i] << nl;
                Info << "Liquid  t: " << phases[1][i] << nl;
                Info << "Solid t  : " << phases[2][i] << nl;

                Info<< nl;
                Info << "Transfer Alpha: " << phases[0][i]+phases[1][i]+phases[2][i] << nl;

                    if(phases[0][i] < 0)
                    {
                      phases[0][i] = 0;
                    }
                    if(phases[1][i] < 0)
                    {
                      phases[1][i] = 0;
                    }
                    if(phases[2][i] < 0)
                    {
                      phases[2][i] = 0;
                    }

               Info<< nl;
               Info << "Gas C    : " << phases[0][i] << nl;
               Info << "Liquid  C: " << phases[1][i] << nl;
               Info << "Solid C  : " << phases[2][i] << nl;

               Info << nl;
               Info << "Corrected Alpha: " << phases[0][i]+phases[1][i]+phases[2][i] << nl;
               Info << nl << nl << nl;

                    if(phases[0][i]+phases[1][i]+phases[2][i] > 1.01)
                    {
                      Info << "FATAL ERROR: CORRECTED ALPHA > 1 Correcting:" << nl;
                      phases[1][i] = double(1) - phases[0][i] - phases[2][i];
                    }
                    if(phases[0][i]+phases[1][i]+phases[2][i] > 1.02)
                    {
                      Info << "FATAL ERROR: APOCALYPSE" << nl;
                      return 0;
                    }
                  }
            }
         }
    }
