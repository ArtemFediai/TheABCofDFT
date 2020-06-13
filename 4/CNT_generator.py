import numpy as np
import argparse 
import copy
import pandas as pd
from pandas import DataFrame as df

def input_uc(uc_file,nper,tube_ax):
    if tube_ax == None or tube_ax not in ['x','y','z']:
        tube_ax = 'y'
    else:
        pass
    with open('CNT_'+nper+'periods.xyz', '+w') as outf:
        A = pd.read_table(uc_file,skiprows=[0,1], sep='\s+',names=['atom','x','y','z'])
        cnt = A
        if tube_ax == 'x':
            for i in range(int(nper)-1):
                uc_shift = 4.26 + i*4.26
                B = copy.copy(A)
                B['x'] = B['x'] + uc_shift
                cnt = cnt.append(B, ignore_index=True)
        elif tube_ax == 'y':
            for i in range(int(nper)-1):
                uc_shift = 4.26 + i*4.26
                B = copy.copy(A)
                B['y'] = B['y'] + uc_shift
                cnt = cnt.append(B, ignore_index=True)
        elif tube_ax == 'z':
            for i in range(int(nper)-1):
                uc_shift = 4.26 + i*4.26
                B = copy.copy(A)
                B['z'] = B['z'] + uc_shift
                cnt = cnt.append(B, ignore_index=True)
        
        natoms = str(len(cnt))
        outf.write(natoms+'\n\n')
        cnt.to_csv(outf, header=None, index=None, sep=' ', mode='a',\
                                               float_format='%.6f')
    outf.close()
            

def main(cnt_unitcell):
    uc_file  = cnt_unitcell.CNT_unit_cell
    nper     = cnt_unitcell.periods
    tube_ax  = cnt_unitcell.axis
    
    input_uc(uc_file,nper,tube_ax)

   
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("CNT_unit_cell")
    parser.add_argument("-p", "--periods", help="Number of periods for the CNT")
    parser.add_argument("-a", "--axis", help="Tube axis",type=str)
    args = parser.parse_args()
    main(args)
