// alter test file for constant density (earth ~2.6 g/cm^3) 
// look at and swipe the mechanisms for finding it from probLinear
// constant slab of stuff
//
// also! integrating the old stuff for the vacuum here

#include <iostream>
#include <fstream>
#include <iomanip>
#include "BargerPropagator.h"

int main(){
	BargerPropagator neutrino;
	neutrino.UseMassEigenstates(false); //check
	
	// osc. parameters - taken from sept. 2026 JUNO paper
	// mixing parameters: sin^2(theta_ij)
	
	bool kSquared; // these are sin^2(theta)
	double theta12 = 0.3036;
	double theta13 = 0.02174797;
	
	// does not effect P_ee in vacuum
	double theta23 = 0.5;
	
	// check these
	double dm21 = 7.388e-5;
	double dm31 = 2.509e-3;
	// for normal. if inverted, dm31 = -2.482e-3
	// use DM2=dm31 directly
	double dm32 = dm31 - dm21;
	int kNuBar = 1; // 1 for neutrino, -1 for antineutrino

	double deltaCP = 0.0;
	
	// fixed baseline (km) - JUNO
	double L = 52.5;
	
	// constant density (g/cm^3) - ~earth's crust density at JUNO
	double density = 2.61;
	
	std::ofstream output("prob3_matter_constant.csv");
	
	if (!output.is_open())
	{
		std::cerr << "Could not create output file. \n";
		return 1;
	}
	
	output << std::setprecision(12);
	
	output  << "Energy_Gev, L_over_E," // possibly change this to MeV
		// add matter (sphere)
		// add matter (varying)
		<< "Pee_matter, Pmumu_matter, Pmue_matter, Pmutau_matter\n"
		<< "Pee_vac, Pmumu_vac, Pmue_vac, Pmutau_vac\n";
		
	// scanning L/E from (check what values we'd be looking at at a 52.5km baseline) in increments of 1
	// LE = L/E
	
	for (double LE = 1; LE <= 4000.0; LE +=1){ // GeV
	
		// corresponding energy
		double E = L / LE;
		// double E_MeV = E_GeV * 1000;
		
		// oscillation params
		neutrino.SetMNS(theta12,theta13,theta23,dm21,dm32,deltaCP,E,kSquared,kNuBar);
		
		// matter (linear)
		// check the propagation tech here
		neutrino.propagateLinear(1 * kNuBar, L, density);
		double Pee_m    = neutrino.GetProb(1,1);
		double Pmue_m   = neutrino.GetProb(2,1);
		double Pmumu_m  = neutrino.GetProb(2,2);
		double Pmutau_m = neutrino.GetProb(2,3);
		
		// vacuum
		double Pee_v    = neutrino.GetVacuumProb(1,1,E,L);
		double Pmue_v   = neutrino.GetVacuumProb(2,1,E,L);
		double Pmumu_v  = neutrino.GetVacuumProb(2,2,E,L);
		double Pmutau_v = neutrino.GetVacuumProb(2,3,E,L);
		
		// double LE = L / E_MeV;
		
		output  << E << ","
			<< LE << ","
			<< Pee_m << "," << Pmue_m << "," << Pmumu_m << "," << Pmutau_m << ","
			<< Pee_v << "," << Pmue_v << "," << Pmumu_v << "," << Pmutau_v << "\n";
	}
	
	output.close();
	
	std::cout << "Done! Wrote prob3_matter_constant.csv\n";
	// std::cout << "dm32 = " << dm32 << "eV^2\n";
	return 0;			
		
}
