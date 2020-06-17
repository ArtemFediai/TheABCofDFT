import numpy as np
import math 
import sympy 
import scipy.interpolate
import scipy.constants as phys
import scipy.integrate as integrate
import decimal
import time 
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import modules.mo as mo
from pandas import DataFrame as df
from functools import partial
from mpmath import mp

class iv_characteristics:

    def __init__(self,config):
        self.Spin     = config["Spin Polarized system"]
        self.NE       = config["Number of energy points"]
        self.Ea       = config["Lower energy border"]
        self.Eb       = config["Upper energy border"]
        self.G_calc   = config["Differential conductance calculation"]
        self.NVsd     = config["Number of sd-voltage points"]
        self.NVgs     = config["Number of gs-voltage points"]
        self.Vsd_a    = config["Lower source-drain voltage"]
        self.Vsd_b    = config["Upper source-drain voltage"]
        self.Vg_a     = config["Lower gate-source voltage"]
        self.Vg_b     = config["Upper gate-source voltage"]
        self.temp     = config["Temperature in Kelvin"]
        self.path_out = config["Path of output"]

    def iv_curves(self):
        path = self.path_out
        if not os.path.exists(path):
            os.makedirs(path)
        
        espace = '   '

        VSD = np.linspace(self.Vsd_a,self.Vsd_b,self.NVsd,dtype=float)
        VG  = np.linspace(self.Vg_a,self.Vg_b,self.NVgs,dtype=float)
        
        #Units eV,s
        e_char = phys.elementary_charge
        hbar   = phys.physical_constants["Planck constant over 2 pi in eV s"][0]
        Planck = phys.physical_constants["Planck constant in eV s"][0]
        kB     = phys.physical_constants["Boltzmann constant in eV/K"][0]
        G0     = phys.physical_constants["conductance quantum"][0]
        print ('electron charge: ', e_char, '\n', 'hbar: ', hbar, '\n', \
               'Boltzman constant[eV/K]: ', kB, '\n', 'Planck constant [eV*s]: ', Planck, \
               '\n', 'Quantum conductance: ', G0 )
        
        #Electrochemical potential in the left electrode
        mu = 0
        #Temperature in Kelvin
        temp = self.temp 
        tempstr = str(temp)
        
        if self.Spin in ["No","no","N","n"]:

            f1 = open(path+"/IVsd"+"-"+tempstr+"K"+".dat","+w")
            
            f1.write("Current" + espace + "Vsd" + espace + espace + "Vg" + espace + "Conductance \n")

            #Read out the energies and transmission values
            f2 = mo.read_files(path+"/out.dat",0,2,0)
            f2.write("# Vsd"+2*espace+"IVsd"+"\n")
            f2.write("# Temperature[K]: "+tempstr+"\n")
            transmission = f2[2]
            E = f2[0]
            
            for voltage in VSD:
                
                landau  = np.zeros(self.NE)
                
                '''
                Integrate Landau equation (using trapezoidal formula)
                '''
                for iE,(energy,trans) in enumerate(zip(E,transmission)):
                    try:
                        fermi_source = mo.fermi_func(energy, mu ,temp)
                        fermi_drain  = mo.fermi_func(energy, mu-voltage,temp)
                        landau[iE]   = trans * (fermi_source - fermi_drain)
                    except OverflowError:
                        fermi_source = float("inf")
                        fermi_drain  = float("inf")
                
                I  = (e_char**2)/Planck * np.trapz(landau,E) 
                
                f1.write("%.8f %.8f\n" %(voltage,I))
                f1.flush()
        

        elif self.Spin in ["Yes","yes","Y","y"]:
   
            if self.G_calc in ["No","no","N","n"]:
                
                #For a+b current
                for spin in ["alpha","beta"]:
                    #Read out the energies and transmission values
                    f2 = mo.read_files(path+"/out-"+spin+".dat",0,2,0)
                    if spin == "alpha":
                        transmission_alpha = f2[2]
                    if spin == "beta":
                        transmission_beta  = f2[2]
                    
                E = f2[0]
                transmission_tot = transmission_alpha + transmission_beta 
                ener_trans = df({'Energy':E, 'Transmission':transmission_tot})
                N_iv         = [0,10,60,160]
                voltage      = []
                result_apb   = []
                result_apbG  = []
                
                f3 = open(path+"/IVsd_a+b.dat","+w")
                f3.write("# Vsd"+2*espace+"IVsd"+"\n")
                f3.write("# Temperature[K]: "+tempstr+"\n")
                
                
                for i,Nvg in enumerate(N_iv):
                    print (-i)
                    endE             = transmission_tot[len(E)-1] 
                    gated_trans      = ener_trans.Transmission.shift(-Nvg * i).fillna(endE)
                    gated_ener_trans = df({'Energy':E, 'Transmission':gated_trans})
                    voltage_g        = round(i * Nvg * (gated_ener_trans.at[1, 'Energy'] -
                                             gated_ener_trans.at[0, 'Energy']), 8)

                    print ("Vg(i): ", i, "\n" ,"vg(eV): ", voltage_g, "\n\n")
                
                    landau       = np.zeros(len(E))
                    der_landau   = np.zeros(len(E))
                    current_apb  = []
                    difcond_apb  = []
                        
                    for vsdi,voltage_sd in enumerate(VSD):
                        for row in gated_ener_trans.itertuples():
                            '''
                            Integrate Landau equation (using trapezoidal formula)
                            '''
                            der_I  = 0.
                            I      = 0.
                            energy = row[1]
                            trans  = row[2]

                            fermi_source = mo.fermi_func(energy,mu,temp)
                            fermi_drain  = mo.fermi_func(energy,mu-voltage_sd/2.,temp)
                            der_I  = mo.exponent_fermi(energy, mu-voltage_sd/2. ,temp)/(1. \
                                   + mo.exponent_fermi(energy, mu-voltage_sd/2. ,temp))**2  
                            
                            landau[row[0]]      = trans * (fermi_source - fermi_drain)
                            der_landau[row[0]]  = trans * der_I
                        
                        I  =  (e_char/Planck) * np.trapz(landau,E) 
                        G_cond = (1./(2. * kB * temp)) * np.trapz(der_landau,E)
                            
                        current_apb.append(I)
                        difcond_apb.append(G_cond)
                        if i == 0:
                            voltage.append(voltage_sd)
                        elif i != 0:
                            pass
                    if i == 0:
                        result_apb.append(voltage)
                    elif i != 0:
                        pass
                    result_apb.append(current_apb)
                    result_apbG.append(difcond_apb)
                
                results_apb  = df(result_apb).T
                results_apbG = df(result_apbG).T
                
                results_apb = pd.concat([results_apb,results_apbG],axis=1)
                #l1 = str(list(Vgate_IV))
                #np.savetxt(f3, results_apb, header='Gate voltages[ev]: '+ l1, delimiter='    ', fmt='% .13f')
                np.savetxt(f3, results_apb, delimiter='    ', fmt='% .13f')
                f3.flush()
                f3.close()


            ### Shifting the Transmission ###
            elif self.G_calc in ["Yes","yes","Y","y"]:
                for spin in ["alpha","beta"]:
                    #Read out the energies and transmission values
                    f2 = mo.read_files(path+"/out-"+spin+".dat",0,2,0)
                    if spin == "alpha":
                        transmission_alpha = f2[2]
                    if spin == "beta":
                        transmission_beta = f2[2]
                    
                    E = f2[0]
                transmission_tot = transmission_alpha + transmission_beta 
                landau        = np.zeros(len(E))
                der_landau    = np.zeros(len(E))
                    
                f1 = open(path+"/Conductance.dat","+w")
                f1.write("# Vg"+2*espace+"Vsd"+2*espace+"Conductance"+"\n")
                f1.write("# Temperature[K]: "+tempstr+"\n")

                ener_trans = df({'Energy':E, 'Transmission':transmission_tot})
                
                for i in range(-self.NVgs+300,self.NVgs):
                    if i >= 0: 
                        endE = transmission_tot[0] 
                    elif i < 0:
                        endE = transmission_tot[len(E)-1] 
                    gated_trans = ener_trans.Transmission.shift(periods=i).fillna(endE)
                    gated_ener_trans = df({'Energy':E, 'Transmission':gated_trans})
                    voltage_g = round(-i * (ener_trans.at[1, 'Energy'] -
                                            ener_trans.at[0, 'Energy']), 8)
                    inter_trans = mo.inter_trans(gated_ener_trans['Energy'],gated_ener_trans['Transmission'])
                    
                    print ("Vg(i): ", i, "\n" ,"vg(eV): ", voltage_g, "\n\n")
                    #Enew = np.linspace(-0.5,0.5,10000)
                    #ynew = inter_trans(Enew)
                    #plt.plot(Enew,ynew,linewidth=1.5)
                    #plt.xticks(np.arange(-0.5,0.5,0.1))
                    #plt.minorticks_on()
                    #plt.show()
                    #sys.exit()
                    for voltage_sd in VSD:
                        
                        ''' Using trapz'''
                        #for row in gated_ener_trans.itertuples():
                        #    der_I  = 0.
                        #    energy = row[1]
                        #    trans  = row[2]
                        #    
                        #    fermi_source = mo.fermi_func(energy,mu,temp)
                        #    fermi_drain  = mo.fermi_func(energy,mu-voltage_sd/2.,temp)
                        #    der_I  = 0.5 * 1./(np.cosh((energy-(mu-voltage_sd/2.))/(kB *  temp)) + 1.)
                        #    der_I  = 0.25 * (1./np.cosh(0.5 * (energy-(mu-voltage_sd/2.))/(kB *  temp)))**2

                        #    landau[row[0]]      = trans * (fermi_source - fermi_drain)
                        #    der_landau[row[0]]  = trans * der_I
                        
                        #der_landau = int_der_landau(gated_ener_trans)
                        #I       = abs((e_char/Planck) * np.trapz(landau,E))
                        #G_cond  = G0 * (1./(2. * kB * temp)) * np.trapz(der_landau,E)

                        #f1.write("%.10f %.10f %.12f\n" %(voltage_g,voltage_sd,I))

                        '''Using integrate.quad'''
                        #int_landau     = lambda x,a,b,c,d: inter_trans(x) * ( \
                        #                         (1./(np.exp((x - a)/(b * c)) + 1.)) \
                        #                       - (1./(np.exp((x - (a-d/2.))/(b * c)) + 1.)) ) 
                        int_landau     = lambda x,a,b,c,d: inter_trans(x) * ( \
                                                 (1./(np.exp(x-a)**(1./(b * c)) + 1.) ) \
                                               - (1./(np.exp(x-(a-d/2.))**(1./b * c) + 1.)) ) 
                        
                        #int_der_landau = lambda y,a,b,c,d: inter_trans(y) * 0.25 * \
                        #                       (\
                        #                       1. - (mp.tanh((y-(a-d/2.))/(2.*b*c)))**2
                        #                       ) 
                        #int_der_landau = lambda y,a,b,c,d: 1. * 0.25 * inter_trans(y) if \
                        #                                   math.isclose(abs(y),abs((a-d/2.)),rel_tol=1E-1) \
                        #                                   else 0
                        #int_der_landau = lambda y,a,b,c,d: inter_trans(y) * 0.5 * \
                        #                       (\
                        #                       1./(1. + mp.cosh((y-(a-d/2.))/(b*c)))
                        #                       )
                        #int_der_landau = lambda y,a,b,c,d: inter_trans(y) * (0.25/(b*c)) * \
                        #                       (\
                        #                       1./(mp.cosh((y-(a-d/2.))/(2.*b*c)))**2
                        #                       )
                        #int_der_landau = lambda y,a,b,c,d: inter_trans(y) * \
                        #                       (\
                        #                       (np.exp(y-(a-d/2.)/(b*c)))/(1. +  \
                        #                        np.exp((y-(a-d/2.))/(b*c)))**2
                        #                       )
                        int_der_landau = lambda y,a,b,c,d: inter_trans(y) * (1./(b*c)) * \
                                               (
                                               (1./(np.exp(y-(a-d/2.))**(1./b * c) + 1.)) * \
                                               ( 1. - 1./(np.exp(y-(a-d/2.))**(1./b * c) + 1.) ) \
                                               )
                        
                        I,err_I       = integrate.quad(int_landau,np.float64(-0.5),np.float64(0.5),\
                                        args=(mu,kB,temp,voltage_sd,),limit=60)
                        I             = (2. * e_char/Planck) * I  
                        G_cond, G_err = integrate.quad(int_der_landau, np.float64(-0.5),
                                        np.float64(0.5),args=(mu,kB,temp,voltage_sd,),limit=60)
                        G_cond        = (2. * e_char/Planck) * G_cond
                        
                        f1.write("%.10f %.10f %.14f %.14f %.14f \n" %(voltage_g,voltage_sd,abs(I),G_cond,I))
                    f1.write("\n")
                    f1.flush()
                f1.close()
