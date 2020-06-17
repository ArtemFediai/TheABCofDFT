import numpy as np
import scipy.constants as phys
import sys
import re


class SystemDef:
    
    def __init__(self,cfg):
        self.Sname              = cfg["System name"]
        self.Spin               = cfg["Spin Polarized system"]
        self.iatomC             = cfg["Central system"][0]["first atom"]
        self.fatomC             = cfg["Central system"][1]["last atom"]
        self.iatomLleads        = cfg["Leads"][0]['first atom left Lead']
        self.fatomLleads        = cfg["Leads"][1]['last atom left Lead']
        self.iatomRleads        = cfg["Leads"][2]['first atom right Lead']
        self.fatomRleads        = cfg["Leads"][3]['last atom right Lead']
        self.Nperiods           = cfg["Number of periods on the Leads"]
        self.CP2K_ao_mat        = cfg["Name of the CP2K ao matrices"]
        self.CP2Koutfile        = cfg["Name of the .out CP2K file"]
        self.CP2Kpdosfile_alpha = cfg["Name of the .pdos alpha CP2K file"]
        self.CP2Kpdosfile_beta  = cfg["Name of the .pdos beta CP2K file"]
        self.path2files         = cfg["Path to CP2K ao matrices"]
        self.path_inp           = cfg["Path to the system 14-matrices"]
        self.pDOS               = cfg["Projected DOS"]
        self.iatom2ign          = cfg["first atom to ignore"]
        self.fatom2ign          = cfg["last atom to ignore"]
   

# Reads the KOHN-SHAM and OVERLAP matrices from the CP2K *.ao file
# It returns the KOHN-SHAM matrix in eV and the OVERLAP matrix.
    def readCP2K2file(self,CP2K_ao_mat):
        overlap_elements = []
        ks_elements      = []
        where_to_put     = overlap_elements
        lines2ignore     = [1,2,3,4]
        with open (CP2K_ao_mat, "r") as f1:
            for line in f1:
                if "OVERLAP MATRIX" in line:
                    continue
                elif "KOHN-SHAM MATRIX" in line:
                    where_to_put = ks_elements
                    continue
                Mat = line.split() 
                if len(Mat) == 0:
                    continue
                elif len(Mat) in lines2ignore:
                    where_to_put.append([])
                    continue
                where_to_put[-1].append(Mat[4:])
                
            Overlap = np.hstack(overlap_elements).astype(float)
            KS      = (phys.physical_constants["Hartree energy in eV"][0])*(np.hstack(ks_elements)).astype(float)
        return Overlap, KS

# Reads specific orbital interactions. For example, if one is only interested in pz-pz
# interaction, just will read those orbitals
    def readCP2K2file_porbit(self,CP2K_ao_mat):
        overlap_elements = []
        ks_elements      = []
        where_to_put     = overlap_elements
        lines2ignore     = [1,2,3,4]
        with open (CP2K_ao_mat, "r") as f1:
            for line in f1:
                if "OVERLAP MATRIX" in line:
                    continue
                elif "KOHN-SHAM MATRIX" in line:
                    where_to_put = ks_elements
                    continue
                Mat    = line.split() 
                if len(Mat) == 0:
                    continue
                elif len(Mat) in lines2ignore:
                    where_to_put.append([])
                    continue
                orbit   = re.search(r".*?(3pz) (.*\d)", line)
                orbit1  = re.search(r".*?(4pz) (.*\d)", line)
                if orbit:
                    #print (orbit)
                    where_to_put[-1].append(list(orbit.group(2).split()))
                #if orbit1:
                    #print (orbit1)
                #    where_to_put[-1].append(list(orbit1.group(2).split()))
                #else:
                #    continue
            Overlap = np.hstack(overlap_elements).astype(float)
            KS      = (phys.physical_constants["Hartree energy in eV"][0])*(np.hstack(ks_elements)).astype(float)
        return Overlap, KS
    
# Reads the OVERLAP and KOHN-SHAM MATRICES (spin ALPHA and BETA)
    def readCP2K2fileSpin(self,CP2K_ao_mat):
        overlap_elements  = []
        ks_elements_alpha = []
        ks_elements_beta  = []
        lines2ignore = [1,2,3,4]
        where_to_put = overlap_elements
        with open (CP2K_ao_mat, "r") as f1:
            for line in f1:
                if "OVERLAP MATRIX" in line:
                    continue
                elif "KOHN-SHAM MATRIX FOR ALPHA SPIN" in line:
                    where_to_put = ks_elements_alpha
                    continue
                elif "KOHN-SHAM MATRIX FOR BETA SPIN" in line:
                    where_to_put = ks_elements_beta
                    continue
                Mat = line.split() 
                if len(Mat) == 0:
                    continue
                elif len(Mat) in lines2ignore:
                    where_to_put.append([])
                    continue
                where_to_put[-1].append(Mat[4:])
                
            Overlap  = np.hstack(overlap_elements).astype(float)
            KS_alpha = (phys.physical_constants["Hartree energy in eV"][0])*(np.hstack(ks_elements_alpha)).astype(float)
            KS_beta  = (phys.physical_constants["Hartree energy in eV"][0])*(np.hstack(ks_elements_beta)).astype(float)
        return Overlap, KS_alpha, KS_beta


# Get the Fermi Energy of the system. It uses as Input file the *.out or *.pdos file from CP2K
    def getEf(self,CP2Koutfile):
        Ef = []
        try:
            with open(CP2Koutfile) as f1:
                for line in f1.readlines():
                    Efs  = re.search(r".*?fermi\s+energy\s?.*?(-?\d+\.\d+).*?", line, flags=re.IGNORECASE)
                    Efs1 = re.search(r"(step i = 0), (\D*) (.*\d)", line)
                    if Efs:
                        Ef.append(float(Efs.group(1)))
                        Ef = [(phys.physical_constants["Hartree energy in eV"][0])*i for i in Ef]
                    if Efs1:
                        Ef.append(float(Efs1.group(3)))
                        Ef = [(phys.physical_constants["Hartree energy in eV"][0])*i for i in Ef]
                if self.Spin in ["No","no","N","n"]:
                    #Ef = [(phys.physical_constants["Hartree energy in eV"][0])*i for i in Ef]
                    Ef = Ef
                elif self.Spin in ["Yes","yes","Y","y"]:
                    Ef = Ef
            return Ef
        except FileNotFoundError:
            print ("Wrong file or path to *.out file!")
            sys.exit()


# This function substracts two things:
# 1) atom type and number of orbitals in those atoms
# 2) geometry of the system
    def get_region_geometry(self,CP2Koutfile):
        atom_kind = []
        number_orbitals = []
        geometry_atom = []
        atom_number = []
        try:
            with open(CP2Koutfile) as f1:
                if self.pDOS in ["No","no","N","n"]:
                    for line in f1.readlines():
                        at_info = re.search("(\d*) (Atomic kind): (.*\D) (Number of atoms): (.*\d)", line)
                        num_orb = re.search(" (Number of spherical basis functions): (.*\d)",line)
                        geo_at  = re.search("(.*\d) (.*\d) ([A-Z]|[A-Z][a-z]) (.*\d) ",line)
                        if at_info:
                            atom_kind.append(str(at_info.group(3)))
                            atom_kind = [element.strip() for element in atom_kind]
                        elif num_orb:
                            number_orbitals.append(int(num_orb.group(2)))
                        elif geo_at:
                            atom_number.append(int(geo_at.group(1)))
                            geometry_atom.append(geo_at.group(3))
                        dic_elements = dict(zip(atom_kind,number_orbitals))
                        geometry     = list(zip(atom_number,geometry_atom))
                elif self.pDOS in ["Yes","yes","Y","y"]:
                    for line in f1.readlines():
                        at_info = re.search("(\d*) (Atomic kind): (.*\D) (Number of atoms): (.*\d)", line)
                        num_orb = re.search(" (Number of spherical basis functions): (.*\d)",line)
                        geo_at  = re.search("(.*\d) (.*\d) ([A-Z]|[A-Z][a-z]) (.*\d) ",line)
                        at2ig   = np.array(range(self.iatom2ign,self.fatom2ign+1),dtype=int)
                        if at_info:
                            atom_kind.append(str(at_info.group(3)))
                            atom_kind = [element.strip() for element in atom_kind]
                        elif num_orb:
                            number_orbitals.append(int(num_orb.group(2)))
                        elif geo_at:
                            if int(geo_at.group(1)) in at2ig:
                                continue
                            elif int(geo_at.group(1)) not in at2ig:
                                atom_number.append(int(geo_at.group(1)))
                                geometry_atom.append(geo_at.group(3))
                        atom_number = list(range(1,len(atom_number)+1))

                        dic_elements = dict(zip(atom_kind,number_orbitals))
                        geometry     = list(zip(atom_number,geometry_atom))
                return dic_elements, geometry

        except FileNotFoundError:
            print ("Wrong file or path to *.out file")
            sys.exit()
