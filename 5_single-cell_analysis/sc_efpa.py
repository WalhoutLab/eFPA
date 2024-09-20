import sys
import pandas as pd

from myCobra import *
from myHdf import hdf
from myIO import *
from MetabolicDistance import *

# INPUT
#rxns=['MAR03905', 'MAR03907'] #[G.matchrxn('ADSL1')] #['MAR03905', 'MAR03907']
#cells=['CTGATCCTCCAGCACG_TSP7_LymphNodes_Inguinal_10X_1_1','TTGAACGAGAAGAGCA_TSP7_Blood_NA_10X_2_1'] #
rxns=sys.argv[1].split(',') 
cells=sys.argv[2].split(',')
fname=sys.argv[3]

# SETTINGS
# see Help in myCobra.mnm.eFPA() and myCobra.mnm.setFPA() for the definitions of these terms
report=False #Make this true to save the log of eFPA analysis inclusing a table of fluxes, distances, and weights
dist_boundary=6
hardBoundary=False
base=2
plus=1.
allowance=1.
correctLoops=False

# FILE PATHS
jsonModelAt='Input/human_gem_solid.json'
byproductsAt='Input/byproducts.pkl'
gem2tsAt='Input/gemGene2TSGene.pkl'
supercellAt='Input/supercell.pkl'
outPath='Output/'
expAt='Input/TS_random100.h5ad'
distAt='Input/rxn_distance_weighted.h5ad'

# LOADING
## Auxiliary files
byproducts=loadobj(byproductsAt)
gem2ts=loadobj(gem2tsAt)
supercell=loadobj(supercellAt)

## Gene expression
print('Loading gene expression')
genes=list(gem2ts)
H=hdf(expAt)
H.setObs('obs/cell_id')
H.setVar('var/ensemblid')
expo=H[cells,:]
vsuper=np.array([supercell[var] for var in H.var]).reshape(1,-1)
expo=np.concatenate([expo,vsuper], axis=0)

cols_list = []
for gene in genes:
    # Calculate the sum of columns associated with this gene
    if gem2ts[gene]:
        cols_sum = np.sum(expo[:, [H.var_index[var] for var in gem2ts[gene]]], axis=1).reshape(-1,1)
    else:
        cols_sum = np.full((expo.shape[0], 1), np.nan) 
    # Append the result to the list
    cols_list.append(cols_sum)
    
exp=np.concatenate(cols_list, axis=1)
exp=exp.transpose()
matexp=matrix(exp,genes,cells+['supercell'])

## Model
print('Loading model')
M=mnm(jsonModelAt,'json')

## Distances
print('Loading distances')
dd=dist_data(distAt)
mdist=dd.getRxns(rxns)
dd.close()

###Inserting missing distances as NaN
Smissing_rxns=set(M.irrevIDs).difference(mdist.rows)
Lmissing_rxns=list(Smissing_rxns)
Nmiss=len(Lmissing_rxns)
row_vals=np.full((Nmiss,mdist.colno),np.nan)
mdist=mdist.appendRows(row_vals,Lmissing_rxns)

# PREPARATION 
print('\nSetting FPA')
M.setFPA(matexp,byproducts=byproducts,inherit_array=mdist,plus=plus)
conds_wo_supercell=[x for x in matexp.cols if not x=='supercell']
    
# RUNNING eFPA
print('Running eFPA reaction by reaction')
dircon={'_f':[0.,1000.],'_r':[-1000.,0.]}
Dresults,Dfinal={},{}
k=0
for rxn in rxns:
    k=k+1
    print('\n\n~~~',k,rxn,'~~~\n')
    #determining directional forms as target reactions
    targets=[]
    for sfx in ('_f','_r'):
        rxnd=rxn+sfx
        if rxnd in M.irrevIDs:
            targets.append(rxnd)
    #going over each target
    for target in targets:
        print(target)
        #overwrite finite boundary of target
        if M.optimizer.Dvars[target].ub<1000.:
            con_spc={rxn:dircon[target[-2:]]}
            print('target will be constrained by '+str(con_spc));
        else:
            con_spc=None;        
        #reporting
        if report:
            reporter=(outPath,'log_'+fname.split('.')[0]+'_','.csv',',')
        else:
            reporter=None                  
        #eFPA
        Dresults[target]=M.eFPA(target,dist_boundary,hardBoundary,base,allowance=allowance,conditions=conds_wo_supercell,\
        specialBoundaries=con_spc,reportingInstructions=reporter,correctDistanceLoops=False)

# SAVE
print('Saving results as dictionary')
saveobj(Dresults,outPath+fname)











