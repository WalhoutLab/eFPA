from MetabolicDistance import *
import numpy as np

def op(fname):
    '''Opens pickled object at <fname>.pkl and returns this object.'''
    fo=open(fname, 'rb');
    obj=pickle.load(fo);
    fo.close(); 
    return obj
    
def writeMatrix2File(A,DI,DJ,txtFile,dlm='\t',reportEvery=500,rowHead='Rowname'):
    '''Writes a numpy matrix <A>, with row dictionary (row name -> index) and column dictionary
    (column name -> index), to a <dlm> delimited file named <textFile> (filename includes path).
    '''
    N=len(DI);
    M=len(DJ);
    Lrows=N*[0];
    Lcols=M*[0];
    for row in DI:
        Lrows[DI[row]]=row;
    for col in DJ:
        Lcols[DI[col]]=col;    
    
    fo=open(txtFile,'w');
    fo.write(dlm.join([rowHead]+Lcols)+'\n');
    
    for i in range(N):
        if not i%reportEvery:
            print str(i)+'\t'+Lrows[i];
        fo.write(Lrows[i]+dlm+dlm.join([str(x) for x in A[i]])+'\n');
    fo.close();

def makeRxnDistanceMatrix(Folder='Output/Reactions/',saveTo='Output/distanceMatrix.txt'):
    '''Uses rxn by rxn results in <Folder> to generate a from rxn to rxn distance matrix.
    The matrix is saved in text format to <saveTo>. 
    '''
    print 'Preparing'
    S,rxns,metabolites,lb,ub,byproducts=loadVarsFromFilenames('Input/Smatrix_regular.npy',\
    'Input/reactions_regular.txt','Input/metabolites_regular.txt','Input/LB_regular.txt',\
    'Input/UB_regular.txt','Input/byproducts_regular.txt');  
    Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts);
    Lrxn=Net.Drxn2node.keys();
    Lrxn.sort();
    N=len(Lrxn);

    DI={};
    for i in range(N):
        DI[Lrxn[i]]=i;

    A=np.inf*np.ones((N,N));
    Lnotfound=[];
    print 'Filling in from-to matrix';
    Lbad=[];
    k=0;
    for rxn in Lrxn:
        k=k+1;
        if not k%500:
            print str(k)+'\t'+rxn;
        try:
            t=op(Folder+rxn+'.pkl');
            Dpath=t[0];
            goOn=True;
        except:
            Lnotfound.append(rxn);
            print 'No tuple file found for reaction '+rxn
            goOn=False;
        if goOn:
            if t[2]:
                print 'Errors exist for '+rxn;
            Lbad.append(len(t[1]));
            for torxn in Dpath:
                Lpath=Dpath[torxn];
                if not Lpath is None:
                    A[DI[rxn],DI[torxn]]=len(Lpath)-1;
    print str(np.sum(Lbad))+' reactions were at the end of questionable paths. Check the returned list.';
    print 'Saving from-to matrix';
    writeMatrix2File(A,DI,DI,saveTo);

    return Lbad,Lnotfound
    
Lbad,Lnotfound=makeRxnDistanceMatrix();
    