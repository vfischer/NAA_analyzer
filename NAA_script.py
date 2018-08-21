####################################################
#
# NAA script
# Use it to compute an expected number of events
# or an exected ppm concentration
# Author: Vincent Fischer
#
####################################################

import math
import optparse

#### GENERAL INPUTS ####
N_Avogadro = 6.022140e23 # Avogadro constant

# Reactor infos
n_flux = 2e11 # Thermal neutron flux in n.cm^-2.s^-1 # Value till 13-02-18 was 1.73e11 (aluminium holder), 2.1e11 now
n_flux_err = 2.1e10 # Error on thermal neutron flux in n.cm^-2.s^-1

# Sample infos
sample_size = 20e-3 # sample size in liters (works for a solution only)
sample_density = 0.8 # sample density (1 for water, 0.8 for DDA, 0.9 for LAB)

# Time infos
t_irrad = 3600 # Irradiation time in seconds
t_cooling = 48*3600 # Cooling time in seconds
t_counting = 3*24*3600 # Counting time in seconds, should be the total time since the livetime is used to correct the peak integral

# HPGe infos
Epsilon_c = 0.01*0.12 # Efficiency of the germanium at this energy in percents norm. to 1

# User-defined element infos
sigma_ng_b_user = 15.4 # Neutron capture cross section (n,g) in barns
#sigma_ng = sigma_ng_b * 1e-24 # Neutron capture cross section (n,g) in cm^-2
I_gamma_user = 0.0991 # Intensity (branching ratio) of the gamma line in percents norm. to 1
element_half_life_user = 27.7*24*3600 # Half-life in seconds
#element_lambda = math.log(2)/element_half_life # lambda in s^-1
molar_mass_user = 52 # molar mass in grams per mol
abundance_user = 0.0435 # abundance of the isotope of interest in percents norm. to 1


# ============================================================================================================================================================================ #
# ============================================================================================================================================================================ #

def load_isotope(isotope_name):

    # Defining global variables (I know it's bad but I ain't no programmer)
    global sigma_ng_b
    global I_gamma
    global element_half_life
    global molar_mass
    global abundance
    global sigma_ng
    global element_lambda

    if isotope_name == "Na23":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 0.517, 0.9999, 14.997*3600, 22.99, 1
        print "Na-23 infos: Looking at the 1368,6 keV peak with 99.99% intensity"
    if isotope_name == "K41":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 1.46, 0.1808, 12.36*3600, 39.1, 0.0673
        print "K-41 infos: Looking at the 1524.6 keV peak with 18.08% intensity"
    if isotope_name == "Mn55":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 13.36, 0.9885, 2.5789*3600, 54.94, 1
        print "Mn-55 infos: Looking at the 846.8 keV peak with 98.85% intensity"
    if isotope_name == "Fe58":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 1.32, 0.565, 44.495*24*3600, 55.845, 0.00282
        print "Fe-58 infos: Looking at the 1099.2 keV peak with 56.5% intensity"
    if isotope_name == "Co59":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 74, 0.9998, 1925.28*24*3600, 58.93, 1
        print "Co-59 infos: Looking at the 1332.5 keV peak with 99.98% intensity"
    if isotope_name == "Cu63":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 4.5, 0.00475, 12.701*3600, 63.55, 0.6915
        print "Cu-63 infos: Looking at the 1345.77 keV peak with 0.475% intensity"
    if isotope_name == "Ni64":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 1.64, 0.2359, 2.5175*3600, 58.69, 0.009255
        print "Ni-64 infos: Looking at the 1481.84 keV peak with 23.59% intensity"
    if isotope_name == "Zn64":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 0.79, 0.5004, 243.93*24*3600, 65.38, 0.4917
        print "Zn-64 infos: Looking at the 1115.54 keV peak with 50.04% intensity"
    if isotope_name == "Br81":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 2.36, 0.834, 35.282*3600, 79.9, 0.4931
        print "Br-81 infos: Looking at the 776.5 keV peak with 83.4% intensity"
    if isotope_name == "Mo98":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 0.13, 0.122, 65.924*3600, 95.94, 0.2439
        print "Mo-98 infos: Looking at the 739.5 keV peak with 12.2% intensity"
    if isotope_name == "Te130":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 0.195, 0.815, 8.0252*24*3600, 129.9, 0.3408
        print "Te-130 infos: Looking at the 364.5 keV peak with 81.5% intensity (of I-131)"
    if isotope_name == "U238":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 2.68, 0.1451, 2.356*24*3600, 238, 0.9927
        print "U-238 infos: Looking at the 277.6 keV peak with 14.51% intensity (of Np-239)"
    if isotope_name == "Th232":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = 7.35, 0.385, 26.975*24*3600, 232, 1
        print "Th-232 infos: Looking at the 311.9 keV peak with 38.5% intensity (of Pa-233)"
    if isotope_name == "user_defined":
        sigma_ng_b, I_gamma, element_half_life, molar_mass, abundance = sigma_ng_b_user, I_gamma_user, element_half_life_user, molar_mass_user, abundance_user
        print "No particular isotope specified: Use the user-defined parameters instead"

    sigma_ng = sigma_ng_b * 1e-24 # Neutron capture cross section (n,g) in cm^-2
    element_lambda = math.log(2)/element_half_life # lambda in s^-1

    print "================================================"

# ============================================================================================================================================================================ #

def get_concentration(peak_integral,peak_integral_err):
    nb_atom_start_count = peak_integral/(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))
    nb_atom_start_cooling = nb_atom_start_count*math.exp(element_lambda*t_cooling)
    nb_atom_start_irrad = nb_atom_start_cooling*element_lambda/(sigma_ng*n_flux*(1-math.exp(-element_lambda*t_irrad)))
    nb_mol_sample = nb_atom_start_irrad/N_Avogadro
    nb_g_sample = nb_mol_sample*molar_mass
    concentration_perL = nb_g_sample/sample_size

    if peak_integral_err != 0:
        nb_atom_start_irrad_max = (peak_integral+peak_integral_err)/(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))*math.exp(element_lambda*t_cooling)/(sigma_ng*(n_flux-n_flux_err)*t_irrad)
        nb_g_sample_max = nb_atom_start_irrad_max*molar_mass/N_Avogadro
        concentration_perL_max = nb_g_sample_max/sample_size
        nb_atom_start_irrad_min = (peak_integral-peak_integral_err)/(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))*math.exp(element_lambda*t_cooling)/(sigma_ng*(n_flux+n_flux_err)*t_irrad)
        nb_g_sample_min = nb_atom_start_irrad_min*molar_mass/N_Avogadro
        concentration_perL_min = nb_g_sample_min/sample_size

    print "                                                "
    print "================================================"
    print "=================== RESULTS ===================="
    print "================================================"
    print "Number of atoms =               %e" % (nb_atom_start_irrad)
    print "Mass [g] =                      %e" % (nb_g_sample)
    print "Concentration [g/L] =           %e           (Reminder: 1 ppm = 1 mg/L)" % (concentration_perL)
    print "Corrected concentration [g/g] = %e" % (concentration_perL*1e-3*(1./sample_density)*(1./abundance))
    print "------------------------------------------------"
    print "Correction includes:"
    print "- Density (d = %f)           ==> Multiply by %e to go from [g/L] to [g/g]" % (sample_density, 1e-3/(sample_density))
    print "- Abundance (abundance = %f) ==> Multiply by %f to obtain the concentration of the element, not the isotope" % (abundance,1./abundance)
    print "================================================"

    if peak_integral_err != 0:
         print "Error range on number of atoms =               [%e; %e]" % (nb_atom_start_irrad_min,nb_atom_start_irrad_max)
	 print "Error range on size [g] =                      [%e; %e]" % (nb_g_sample_min,nb_g_sample_max)
	 print "Error range on concentration [g/L] =           [%e; %e]           (Reminder: 1 ppm = 1 mg/L)" % (concentration_perL_min,concentration_perL_max)
	 print "Corrected error range on concentration [g/g] = [%e; %e]" % (concentration_perL_min*1e-3*(1./sample_density)*(1./abundance),concentration_perL_max*1e-3*(1./sample_density)*(1./abundance))
         print "================================================"

# ============================================================================================================================================================================ #

def get_peak_integral(concentration,concentration_err):
    nb_g_sample = concentration*sample_size
    nb_mol_sample = nb_g_sample/molar_mass
    nb_atom_start_irrad = nb_mol_sample*N_Avogadro
    nb_atom_start_cooling = nb_atom_start_irrad*(sigma_ng*n_flux)/(element_lambda/(1-math.exp(-element_lambda*t_irrad)))
    nb_atom_start_count = nb_atom_start_cooling/math.exp(element_lambda*t_cooling)
    peak_integral = nb_atom_start_count*(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))

    if concentration_err != 0:
        peak_integral_max = (concentration+concentration_err)*sample_size/molar_mass*N_Avogadro\
                *(sigma_ng*(n_flux+n_flux_err)*t_irrad)/math.exp(element_lambda*t_cooling)*(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))
        peak_integral_min = (concentration-concentration_err)*sample_size/molar_mass*N_Avogadro\
                *(sigma_ng*(n_flux-n_flux_err)*t_irrad)/math.exp(element_lambda*t_cooling)*(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))

    print "                                                "
    print "================================================"
    print "=================== RESULTS ===================="
    print "================================================"
    print "Nb of expected counts under photopeak (peak integral) = %e (efficiency is %f percents)" % (peak_integral,Epsilon_c*100)
    print "================================================"

    if concentration_err != 0:
        print "Range of nb of expected counts under photopeak (peak integral) = [%e; %e] (efficiency is %f percents)" % (peak_integral_min,peak_integral_max,Epsilon_c*100)
        print "================================================"

# ============================================================================================================================================================================ #

def get_neutron_flux(concentration, peak_integral,concentration_err,peak_integral_err):
    nb_g_sample = concentration*sample_size
    nb_mol_sample = nb_g_sample/molar_mass
    nb_atom_start_irrad = nb_mol_sample*N_Avogadro
    nb_atom_start_count = peak_integral/(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))
    nb_atom_start_cooling = nb_atom_start_count*math.exp(element_lambda*t_cooling)
    n_flux = nb_atom_start_cooling*element_lambda/(sigma_ng*nb_atom_start_irrad*(1-math.exp(-element_lambda*t_irrad)))

    if concentration_err != 0 or peak_integral_err != 0:
        n_flux_max = ((peak_integral+peak_integral_err)/(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))*math.exp(element_lambda*t_cooling))\
                /(sigma_ng*t_irrad*(concentration-concentration_err)*sample_size/molar_mass*N_Avogadro)
        n_flux_min = ((peak_integral-peak_integral_err)/(Epsilon_c*I_gamma*(1-math.exp(-element_lambda*t_counting)))*math.exp(element_lambda*t_cooling))\
                /(sigma_ng*t_irrad*(concentration+concentration_err)*sample_size/molar_mass*N_Avogadro)

    print "                                                "
    print "================================================"
    print "=================== RESULTS ===================="
    print "================================================"
    print "Expected neutron flux = %e n.cm^-2.s^-1" % n_flux
    print "================================================"

    if concentration_err != 0 or peak_integral_err != 0:
        print "Range of expected neutron flux = [%e; %e] n.cm^-2.s^-1" % (n_flux_min,n_flux_max)
        print "================================================"

# ============================================================================================================================================================================ #
# ============================================================================================================================================================================ #

if __name__=="__main__":
    parser = optparse.OptionParser(description='Computes the concentration/peak integral of an NAA measurement, python NAA_script.py <options>')
    parser.add_option('-q', dest="quantity", help='Quantity to compute: \"concentration\", \"integral\" or \"flux\"')
    parser.add_option('-i', dest="isotope", help='Use a known isotope')
    parser.add_option('--error', action='store_true', help='Enable error propagation')
    (options, args) = parser.parse_args()

    quantity_list = ["concentration","integral","flux"]
    quantity_name = options.quantity
    isotope_list = ["Na23","K41","Mn55","Fe58","Co59","Cu63","Ni64","Zn64","Br81","Mo98","Te130","U238","Th232","user_defined"]
    isotope_name = options.isotope
    error_bool = options.error
    if quantity_name not in quantity_list:
        parser.print_help()
        raise Exception("Quantity not recognized ! Only \"concentration\", \"integral\" or \"flux\" are allowed !")
    if isotope_name not in isotope_list:
        parser.print_help()
        raise Exception("Isotope not recognized ! Only \"Na23\",\"K41\",\"Mn55\",\"Fe58\",\"Co59\",\"Cu63\",\"Ni64\",\"Zn64\",\"Br81\",\"Mo98\",\"Te130\",\"U238\",\"Th232\", \"user_defined\" are allowed for now!")

    if quantity_name == "concentration":
        print "================================================"
        print "================ ISOTOPE INFOS ================="
        print "================================================"
        print ("Using %s as the isotope of interest" % isotope_name)
        load_isotope(isotope_name)
        print "                                                "
        print "================================================"
        print "================= USER INPUTS =================="
        print "================================================"
        peak_integral = input("Enter the peak integral: ")
        if error_bool:
            peak_integral_err = input("Enter the peak integral error: ")
        else:
            peak_integral_err = 0
        print ("Peak integral is %d +/- %d" % (peak_integral,peak_integral_err))
        print "================================================"
        get_concentration(peak_integral,peak_integral_err)

    if quantity_name == "integral":
        print "================================================"
        print "================ ISOTOPE INFOS ================="
        print "================================================"
        print ("Using %s as the isotope of interest" % isotope_name)
        load_isotope(isotope_name)
        print "                                                "
        print "================================================"
        print "================= USER INPUTS =================="
        print "================================================"
        concentration = input("Enter the concentration in g/L (reminder: 1 ppm = 1e-3 g/L): ")
        if error_bool:
            concentration_err = input("Enter the concentration error in g/L (reminder: 1 ppm = 1e-3 g/L): ")
        else:
            concentration_err = 0
        print ("Concentration is %e +/- %e g/L (or %e +/- %e ppm)" % (concentration,concentration_err,concentration*1e3,concentration_err*1e3))
        print "================================================"
        get_peak_integral(concentration,concentration_err)

    if quantity_name == "flux":
        print "================================================"
        print "================ ISOTOPE INFOS ================="
        print "================================================"
        print ("Using %s as the isotope of interest" % isotope_name)
        load_isotope(isotope_name)
        print "                                                "
        print "================================================"
        print "================= USER INPUTS =================="
        print "================================================"
        concentration = input("Enter the concentration in g/L (reminder: 1 ppm = 1e-3 g/L): ")
        peak_integral = input("Enter the peak integral: ")
        if error_bool:
            concentration_err = input("Enter the concentration error in g/L (reminder: 1 ppm = 1e-3 g/L): ")
            peak_integral_err = input("Enter the peak integral error: ")
        else:
            concentration_err = 0
            peak_integral_err = 0
        print ("Concentration is %e +/- %e g/L (or %e +/- %e ppm)" % (concentration,concentration_err,concentration*1e3,concentration_err*1e3))
        print ("Peak integral is %d +/- %d" % (peak_integral,peak_integral_err))
        print "================================================"
        get_neutron_flux(concentration,peak_integral,concentration_err,peak_integral_err)



