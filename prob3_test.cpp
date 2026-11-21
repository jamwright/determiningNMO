#include <iostream>
#include <fstream>
#include <iomanip>
#include "BargerPropagator.h"

int main(){
	BargerPropagator neutrino;
	
	// osc. parameters - taken from sept. 2026 JUNO paper
	// mixing parameters: sin^2(theta_ij)
	
	double theta12 = 0.3036;
	double theta13 = 0.02174797;
	
	// does not effect P_ee in vacuum
	double theta23 = 0.5;
	
	// check these
	double dm21 = 7.388e-5;
	double dm31 = 2.509e-3;
	// for normal. if inverted, dm31 = -2.482e-3
	double dm32 = dm31 - dm21;

	double deltaCP = 0.0;
	
	// fixed baseline (km) - JUNO
	double L = 52.5;
	
	std::ofstream output("prob3_test.csv");
	
	if (!output.is_open())
	{
		std::cerr << "Could not create output file. \n";
		return 1;
	}
	
	output << std::setprecision(12);
	
	output  << "Energy_Mev, L_over_E,"
		<< "P_ee, P_mumu, P_mue\n";
		
	// scanning L/E from (check what values we'd be looking at at a 52.5km baseline) in increments of .1
	// LE = L/E
	
	for (double LE = 1; LE <= 1000.0; LE +=0.1){
	
		// corresponding energy
		double E_GeV = L / LE;
		double E_MeV = E_GeV * 1000;
		
		//oscillation params
					neutrino.SetMNS(theta12,theta13,theta23,dm21,dm32,deltaCP,E_GeV,true,-1);
		
		//calculate probabilities
		double P_ee = neutrino.GetVacuumProb(1,1,E_GeV,L);
		double P_mumu = neutrino.GetVacuumProb(2,2,E_GeV,L);
		double P_mue = neutrino.GetVacuumProb(2,1,E_GeV,L);
		
		// double LE = L / E_MeV;
		
		// writing to csv
		
		output  << E_MeV << ","
			<< LE << ","
			<< P_ee << ","
			<< P_mumu << ","
			<< P_mue << "\n";
	}
	
	output.close();
	
	std::cout << "Prob3 L/E calculation complete. \n";
	std::cout << "dm32 = " << dm32 << "eV^2\n";
	return 0;	
				
		
}
