const char filename[120] = "/home/nikolas/Documents/JRCnoSUBM.dat";

SetThermoFilename(filename);
FORTRAN_CALL(setstandardunits)();

int iElement = 0;
double dMass = 0.0;

FORTRAN_CALL(setelementmass)(&iElement, &dMass);
iElement = 9; //F2
dMass = 2*(phases[0].thermoRef().p()[i]*phases[0][i]*cv[i])/(phases[0].thermo().T()[i]*8.31446261815324); //gas F2
Info<< "dMass F2: " << 2*(phases[0].thermoRef().p()[i]*phases[0][i]*cv[i])/(phases[0].thermo().T()[i]*8.31446261815324) << nl;
FORTRAN_CALL(setelementmass)(&iElement, &dMass);
iElement = 55; //Cs
dMass = (phases[1][i]*cv[i]*1873)/0.132905; //liquid
Info<< "dMass Cs: " << (phases[1][i]*cv[i]*1873)/0.132905 << nl;
FORTRAN_CALL(setelementmass)(&iElement, &dMass);


double dTemperature = phases[1].thermo().T()[i];
double dPressure = phases[0].thermoRef().p()[i]/1e5;   //[atm]
FORTRAN_CALL(settemperaturepressure)(&dTemperature, &dPressure);
FORTRAN_CALL(ssparsecsdatafile)();
FORTRAN_CALL(thermochimica)();

int iflag;
FORTRAN_CALL(checkinfothermo)(&iflag);
//printf("Debug flag: %i\n",iflag);
int iPrintMode = 2;

FORTRAN_CALL(setprintresultsmode)(&iPrintMode);
//FORTRAN_CALL(printresults)();

int info;
FORTRAN_CALL(checkinfothermo)(&info);
//printf("Information: %i\n",info);


  #include "thermochimicaCoupling/gaseousThermochimica.H"   //phasei = 0
  #include "thermochimicaCoupling/liquidusThermochimica.H"  //phasei = 1
  #include "thermochimicaCoupling/solidusThermochimica.H"   //phasei = 2

  //Info << nl << nl;

FORTRAN_CALL(resetthermoall)();
