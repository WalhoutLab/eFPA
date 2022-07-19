from MetabolicDistance import *
import os
import time

def fcluster(cmd,q='large',n=20,nhost=1,W='2:00',mem=2000,outFile='temp.txt'):
    '''Returns proper cluster command for running a job with the specified input.
    <cmd> is the command that describes the job.'''
    Ls=['bsub'];
    if n:
        Ls.append('-n '+str(n));
    if mem:
        Ls.append('-R rusage[mem='+str(mem)+']');
    if W:
        Ls.append('-W '+str(W));
    if nhost:
        Ls.append('-R span[hosts='+str(nhost)+']');
    Ls.append('-q '+q);    
    if outFile:
        Ls.append('-o '+'"'+outFile+'"');
    Ls.append('\''+cmd+'\'');
    return ' '.join(Ls)

def cluster_pathTo(targetFolder,outFolder='Nohup/',progName='exampleCluster_function.py',q='short',\
    ncore=1,W='1:00',mem=2000,testN=None,quiet=False,pause=None): 
    '''High throughput usage of pathTo() in computer cluster. Use <testN> as a list
    of reaction inputs for test, or as an integer to test for first <textN> reactions.'''
    if not testN or type(testN)==int:
        S,rxns,metabolites,lb,ub,byproducts=loadVarsFromFilenames('Input/Smatrix_regular.npy',\
        'Input/reactions_regular.txt','Input/metabolites_regular.txt','Input/LB_regular.txt',\
        'Input/UB_regular.txt','Input/byproducts_regular.txt');  
        Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts);
        Lrxn=Net.Drxn2node.keys();
        Lrxn.sort();
    if testN:
        if type(testN)==list:
            Lrxn=testN;
        else:
            Lrxn=Lrxn[:testN];
    k=0;
    for rxn in Lrxn:
        k=k+1;
        cmd='python '+progName+' '+rxn;
        cmdf=fcluster(cmd,q=q,n=ncore,W=W,mem=mem,outFile=targetFolder+outFolder+'out_'+rxn+'.txt');
        os.system(cmdf);
        if not quiet:
            print str(k)+': '+cmdf
        if pause:
            time.sleep(pause);
            
cluster_pathTo('Output/',testN=None,pause=0.0025);