# -*- coding: utf-8 -*-

import numpy as np
import math 
import scipy.interpolate
import scipy.constants as phys
import mpmath as mp

kB     = phys.physical_constants["Boltzmann constant in eV/K"][0]

#simple recursion expression for periodic systems
##########################################################################################
def gs_rec(E,eta,H_00,H_01,niterations):
    try:  
        Gs = np.linalg.inv((E+1j*eta)*np.eye(len(H_00))-H_00) 
        for i in range(niterations):
            Gs = np.linalg.inv(E*np.eye(len(H_00))-H_00-np.matrix.getH(H_01)@Gs@H_01)
    except TypeError:
        Gs = 1/((E+1.j*eta)-H_00)
        for i in range(niterations):
            Gs = 1./(E+1.j*eta-H_00-np.conjugate(H_01)*Gs*H_01)
        sigma = np.conjugate(H_01)*Gs*H_01
    return sigma  
##########################################################################################       

#sancho method
##########################################################################################    
#@jit(nopython=True)
def sancho(energy,h,t0_matrix,sh,st,eps):
    es = energy*sh-h
    e = energy*sh-h
    a = energy*st-t0_matrix
    b = energy*np.matrix.getH(st) - np.matrix.getH(t0_matrix)
    
    while (np.linalg.norm(abs(a), ord='fro') > eps):
        g = np.linalg.inv(e)
        bga = b @ g @ a
        agb = a @ g @ b
        e = e - bga - agb
        es = es - agb

        #a = -a @ g @ a
        #b = -b @ g @ b
        a = a @ g @ a
        b = b @ g @ b

    G = np.linalg.inv(es)
    return G

##########################################################################################    
def c_sancho(energy,h,t0_matrix,sh,st,eps):
    es = np.array(energy*sh-h,dtype=complex)
    e = np.array(energy*sh-h,dtype=complex)
    a = np.array(energy*st-t0_matrix,dtype=complex)
    b = np.array(energy*np.matrix.getH(st) - np.matrix.getH(t0_matrix),dtype=complex)
    
    while (np.linalg.norm(abs(a), ord='fro') > eps):
        g = np.linalg.inv(e)
        bga = b @ g @ a
        agb = a @ g @ b
        e = e - bga - agb
        es = es - agb

        a = -a @ g @ a
        b = -b @ g @ b

    G = np.linalg.inv(es)
    return G

##########################################################################################
def sancho_scalar(energy,h,t,sh,st,eps):
    es =np.array([energy*sh-h],dtype=np.complex)
    e = np.array([energy*sh-h],dtype=np.complex)
    a = np.array([energy*st-t],dtype=np.complex)
    b = np.array([energy*np.conjugate(st)-np.conjugate(t)],dtype=np.complex)
    while(abs(a) > eps):
        g = 1./e
        bga = b*g*a
        agb = a*g*b
        e = e - bga -agb
        es = es - agb
        
        a = -a*g*a
        b = -b*g*b
    G = 1./es

    return G  

##########################################################################################
def fermi_func(energy,e1,temp):
    kB = phys.physical_constants["Boltzmann constant in eV/K"][0]
    fermi_func = 1./(np.exp((energy - e1)/(kB * temp)) + 1)

    return fermi_func

##########################################################################################
def exponent_fermi(energy,e1,temp):
    kB = phys.physical_constants["Boltzmann constant in eV/K"][0]
    exponent_f = np.exp((energy - e1)/(kB * temp))

    return exponent_f

##########################################################################################
def inter_trans(x,y):
    interpolated_T = scipy.interpolate.interp1d(x,y,kind='quadratic')

    return interpolated_T

##########################################################################################
#def int_der_landau(x,mu,kB,temp,voltage_sd):
#    if math.isclose(y,(a-d/2.),abs_tol=1E-3): 
#        inter_trans(y) * 0.25 * \
#        1.
#    else:
#        0
#
#    return der_fermi

##########################################################################################
def read_files(filename, ctr1, ctr2, lines_2_skip):
    cols        =  [[] for i in range(ctr1, ctr2 + 1)]
    col_2_read  =  np.arange(ctr1, ctr2 + 1)
    cols = np.loadtxt(filename, skiprows=lines_2_skip,
    usecols=col_2_read)

    return np.transpose(cols)

##########################################################################################
