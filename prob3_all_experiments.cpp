// looking at four cases: JUNO, T2K, NOvA, IceCube.
// we consider three models:
// - vacuum : GetVacuumProb
// - linear : propagateLinear (constant density)
// - DefinePath + propagate (earth layers)
// output: one row per (experiment, model, nubar, hierarchy, baseline, E)

#include <iostream>
#include <fstream>
#include <iomanip>
#include "BargerPropagator.h"
#include <cmath>
#include <string>

// CONSTANTS: 

const bool kSquared = true; // these are sin^2(theta), mark as true
const double theta12 = 0.3036;
const double theta13 = 0.02174797;
// does not effect P_ee in vacuum
const double theta23 = 0.5;	
// check these
double dm21 = 7.388e-5; //eV^2
double dm31 = 2.509e-3; //eV^2
// for normal. if inverted, dm31 = -2.482e-3
double dm32 = dm31 - dm21;
// int kNuBar = 1; // 1 for neutrino, -1 for antineutrino
double deltaCP = 0.0;
static const double R_Earth = 6371.0; //km, radius of earth

const double prodHeight = 0.0; //km, for atmospheric neutrinos
// fixed baseline (km) - JUNO
// double L = 52.5;
	
// constant density (g/cm^3) for linear model - ~earth's crust density at JUNO
double density = 2.61;
	
// write one line of the csv file:

void writeRow(std::ofstream& out, std::string experiment, std::string model, int nubar, int hierarchy, 
              double cosZ, double L, double E, double Pee, double Pmumu, double Pmue, double Pmutau)

{
  out << experiment << "," << model << "," << nubar << "," << hierarchy << "," << cosZ  << "," << L << "," << E << "," << L/E << "," << Pee << "," << Pmumu << "," << Pmue << "," << Pmutau << "\n";
}

// scanning energy for one experiment, baseline and neutrino type
// all types for both mass orderings in one thing

void runExperiment(BargerPropagator& neutrino, std::ofstream& out, std::string name, double L, double Emin, double Emax, int nPoints, double density, int nubar) {
  std::cout << "Running " << name << " L = " << L << " km nubar = " << nubar << "\n";

// path through earth is a chord -> zenith angle follows from L
// cosZ = -1 for straight through centre

  double cosZ = -L / (2.0 * R_EARTH); //check

  // hierarchy: +1 NO, then -1 IO
  for (int hierarchy = 1; hierarchy >= -1; hierarchy = hierarchy - 2){
    double DM2 = hierarchy * dm32 //currently, but to be more accurate i will set the const. in here later
    for (int i=0; i < nPoints; i++){
      // energies spaced evenly on log scale from Emin to Emax
      double E = Emin * std::pow(Emax / Emin, double(i) / double(nPoints-1));

      // VACUUM:

      neutrino.SetMNS(theta12,theta13,theta23,dm21,DM2,deltaCP,E,kSquared,nubar);
      writeRow(out, name, "vacuum", nubar, hierarchy, cosZ, L, E,
                neutrino.GetVacuumProb(1,1,E,L),  //Pee
                neutrino.GetVacuumProb(2,2,E,L),  //Pmumu
                neutrino.GetVacuumProb(2,1,E,L),  //Pmue
                neutrino.GetVacuumProb(2,3,E,L)); //Pmutau

      // CONSTANT DENSITY:

      neutrino.SetMNS(theta12,theta13,theta23,dm21,DM2,deltaCP,E,kSquared,nubar);
      neutrino.propagateLinear(1 * nubar, L, density);
      writeRow(out, name, "linear", nubar, hierarchy, cosZ, L, E,
              neutrino.GetProb(1,1),  //Pee
              neutrino.GetProb(2,2),  //Pmumu
              neutrino.GetProb(2,1),  //Pmue
              neutrino.GetProb(2,3)); //Pmutau

      // SPHERE: 
      // i am least confident on this, so this will need to be checked later

      neutrino.SetMNS(theta12,theta13,theta23,dm21,DM2,deltaCP,E,kSquared,nubar);
      neutrino.DefinePath(cosZ,prodHeight); //set to 0 rn
      neutrino.propagate(1*nubar);
      writeRow(out, name, "sphere", nubar, hierarchy, cosZ, L, E, 
              neutrino.GetProb(1,1),  //Pee
              neutrino.GetProb(2,2),  //Pmumu
              neutrino.GetProb(2,1),  //Pmue
              neutrino.GetProb(2,3)); //Pmutau
    }
  }
}

int main(){
  BargerPropagator neutrino;
  neutrino.UseMassEigenstates(false); //as in probLinear

  std::ofstream out("prob3_all_experiments.csv");
  if (!out.is_open()){ 
    std::cerr << "Could not create output file.\n";
    return 1;
  }
  out << std::setprecision(8); // 8 digits
  out << "experiment,model,nubar,hierarchy,cosZ,L_km,E_GeV,L_over_E,"
      << "Pee,Pmumu,Pmue,Pmutau\n";

  // one per experiment, using baselines and the energy spectrums we want to look at.
  //                            name   L(km)  Emin   Emax    N    density  nubar

  // JUNO is reactor anti nu_e, looking at abt 1-10MeV, baseline of 52.5km

  // DENSITY: 2.5g/cm-3, from the paper:
  // "Terrestrial Matter Effects on Reactor Antineutrino Oscillations: Constant vs. Fluctuated Density Profiles"
  // where i have averaged the two densities found to get my value
  const double densityJUNO = (2.541 + 2.554)/2
  runExperiment(neutrino, out, "JUNO", 52.5, 0.001, 0.010, 5000, densityJUNO, -1);

  // T2K: 0.1-3 GeV, baseline of 295km, nubar = +-1
  // DENSITY: 2.6g/cm-3, from the paper:
  // "the earth matter effects in neutrino oscillation experiments from Tokai to Kamioka and Korea"

  runExperiment(neutrino, out, "T2K", 295.0, 0.1, 3.0, 5000, 2.6, 1);
  runExperiment(neutrino, out, "T2K", 295.0, 0.1, 3.0, 5000, 2.6, -1);

  // NOvA: 0.5-5GeV, baseline of 810km, nubar=+-1, 
  // DENSITY: 2.84 g/cm^2 from CRUST2.0 model for avg depth of NuMI beam

  runExperiment(neutrino, out, "NOvA", 810.0, 0.5, 5.0, 5000, 2.84, 1);
  runExperiment(neutrino, out, "NOvA", 810.0, 0.5, 5.0, 5000, 2.84, 1);

  // IceCube: one baseline per cos(zenith).
  // taking 5.5g/cm-3 to be mean earth
  double icCosZ[5] = {-1.0, -0.8, -0.6, -0.4, -0.2};
  for (int k=0; k<5; k++){
    double L = -2.0 * R_Earth * icCosZ[k]; //chord length, km

    runExperiment(neutrino, out, "IceCube", L, 2.0, 100.0, 5000, 5.5, 1);
    runExperiment(neutrino, out, "IceCube", L, 2.0, 100.0, 5000, 5.5, -1);
  }

  out.close();
  std::cout << "Done. Wrote prob3_all_experiments.csv!\n";
  return 0;
}









