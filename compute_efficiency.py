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
efficiency_HPGe08pct_face = [8.84, 7.98, 6.19, 2.59, 1.42, 0.94, 0.74, 0.64, 0.46]
efficiency_HPGe25pct_face = [15.90, 14.00, 10.60, 5.26, 3.24, 2.21, 1.74, 1.55, 1.15]
efficiency_HPGe50pct_face = [15.40, 17.30, 15.40, 8.78, 5.83, 3.91, 3.16, 2.85, 2.21]
efficiency_HPGe99pct_face = [10.50, 14.90, 15.70, 10.40, 7.59, 5.11, 4.29, 3.92, 3.18]

# Geometry: Vial @ 10 cm
efficiency_HPGe08pct_10cm = [0.58, 0.533, 0.41, 0.15, 0.0843, 0.057, 0.0463, 0.0405, 0.0283]
efficiency_HPGe25pct_10cm = [0.966, 0.831, 0.709, 0.348, 0.224, 0.167, 0.139, 0.125, 0.096]
efficiency_HPGe50pct_10cm = [0.962, 1.040, 0.990, 0.567, 0.392, 0.316, 0.261, 0.238, 0.184]
efficiency_HPGe99pct_10cm = [1.180, 1.480, 1.480, 0.989, 0.730, 0.617, 0.525, 0.486, 0.391]


# A simple linear interpolation function
def interpolate(energy, En_low, En_high, Eff_low, Eff_high):
    efficiency = Eff_low + (energy-En_low)*(Eff_high - Eff_low)/(En_high - En_low)
    return efficiency


if __name__=="__main__":
    parser = optparse.OptionParser(description='Computes an HPGe efficiency at a certain energy, python compute_efficiency.py <options>')
    parser.add_option('-e', dest="energy", help='Energy of interest [keV]')
    parser.add_option('-g', dest="geom", help='Geometry in the HPGe: Face, 10cm')
    parser.add_option('-t', dest="hpge", help='HPGe type (%): 08, 25, 50, 99')
    (options, args) = parser.parse_args()

    hpge_type = options.hpge
    geometry = options.geom
    energy_interest = float(options.energy)
    hpge_type_list = ["08","25","50","99"]
    geom_type_list = ["Face","10cm"]
    
    # Sanity tests
    if hpge_type not in hpge_type_list:
        parser.print_help()
        raise Exception("HPGe type not recognized ! Only 08, 25, 50, 99 are allowed !")
    if geometry not in geom_type_list:
        parser.print_help()
        raise Exception("Geometry not recognized ! Only Face and 10cm are allowed !")  
    if energy_interest < calibration_energy[0] or energy_interest > calibration_energy[len(calibration_energy)-1]:
        raise Exception("Requested energy out of the calibration bounds ! Efficiency cannot be extrapolated !")
      
    # Loads the efficiency lists depending on the HPGe type  
    if hpge_type == "08":
        if geometry == "Face":
            efficiency_list = efficiency_HPGe08pct_face
        if geometry == "10cm":
            efficiency_list = efficiency_HPGe08pct_10cm
    if hpge_type == "25": 
        if geometry == "Face":
            efficiency_list = efficiency_HPGe25pct_face
        if geometry == "10cm":
            efficiency_list = efficiency_HPGe25pct_10cm
    if hpge_type == "50": 
        if geometry == "Face":
            efficiency_list = efficiency_HPGe50pct_face
        if geometry == "10cm":
            efficiency_list = efficiency_HPGe50pct_10cm
    if hpge_type == "99": 
        if geometry == "Face":
            efficiency_list = efficiency_HPGe90pct_face
        if geometry == "10cm":
            efficiency_list = efficiency_HPGe99pct_10cm
    
    # Loops on the calibration energies to find the values above and below the energy of interest
    for i in range(len(calibration_energy)):
        if calibration_energy[i] < energy_interest and calibration_energy[i+1] > energy_interest:
	    En_low = calibration_energy[i]
	    En_high = calibration_energy[i+1]
	    Eff_low = efficiency_list[i]
            Eff_high = efficiency_list[i+1]
            
    # The desired efficiency value        
    final_efficiency = interpolate(energy_interest, En_low, En_high, Eff_low, Eff_high)   
    
    print "Efficiency of the %s HPGe (%s geometry) at %.3f keV is %.3f percents)" % (hpge_type, geometry, energy_interest, final_efficiency)
    
############ END ################    