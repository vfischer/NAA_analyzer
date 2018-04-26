####################################################
#
# Compute efficiency for NAA
# Use it to compute the efficiency of an HPGe
# Author: Vincent Fischer
#
####################################################

import math
import optparse


### Efficiencies from data ###

calibration_energy = [88.00, 122.10, 165.90, 391.70, 661.70, 898.00, 1173.20, 1332.50, 1836.10] 

# Geometry: Vial Face
efficiency_HPGe08pct = [8.84, 7.98, 6.19, 2.59, 1.42, 0.94, 0.74, 0.64, 0.46]
efficiency_HPGe25pct = [15.90, 14.00, 10.60, 5.26, 3.24, 2.21, 1.74, 1.55, 1.15]
efficiency_HPGe50pct = [15.40, 17.30, 15.40, 8.78, 5.83, 3.91, 3.16, 2.85, 2.21]
efficiency_HPGe99pct = [10.50, 14.90, 15.70, 10.40, 7.59, 5.11, 4.29, 3.92, 3.18]


# A simple linear interpolation function
def interpolate(energy, En_low, En_high, Eff_low, Eff_high):
    efficiency = Eff_low + (energy-En_low)*(Eff_high - Eff_low)/(En_high - En_low)
    return efficiency


if __name__=="__main__":
    parser = optparse.OptionParser(description='Computes an HPGe efficiency at a certain energy, python compute_efficiency.py <options>')
    parser.add_option('-e', dest="energy", help='Energy of interest [keV]')
    parser.add_option('-t', dest="hpge", help='HPGe type (%): 08, 25, 50, 99')
    (options, args) = parser.parse_args()

    hpge_type = options.hpge
    energy_interest = float(options.energy)
    hpge_type_list = ["08","25","50","99"]
    
    # Sanity tests
    if hpge_type not in hpge_type_list:
        parser.print_help()
        raise Exception("HPGe type not recognized ! Only 08, 25, 50, 99 are allowed !")
    if energy_interest < calibration_energy[0] or energy_interest > calibration_energy[len(calibration_energy)-1]:
        raise Exception("Requested energy out of the calibration bounds ! Efficiency cannot be extrapolated !")
      
    # Loads the efficiency lists depending on the HPGe type  
    if hpge_type == "08":
        efficiency_list = efficiency_HPGe08pct
    if hpge_type == "25": 
        efficiency_list = efficiency_HPGe25pct
    if hpge_type == "50": 
        efficiency_list = efficiency_HPGe50pct
    if hpge_type == "99": 
        efficiency_list = efficiency_HPGe99pct 
    
    # Loops on the calibration energies to find the values above and below the energy of interest
    for i in range(len(calibration_energy)):
        if calibration_energy[i] < energy_interest and calibration_energy[i+1] > energy_interest:
	    En_low = calibration_energy[i]
	    En_high = calibration_energy[i+1]
	    Eff_low = efficiency_list[i]
            Eff_high = efficiency_list[i+1]
            
    # The desired efficiency value        
    final_efficiency = interpolate(energy_interest, En_low, En_high, Eff_low, Eff_high)   
    
    print "Efficiency of the %s HPGe at %.3f keV is %.3f percents)" % (hpge_type, energy_interest, final_efficiency)
    
############ END ################    