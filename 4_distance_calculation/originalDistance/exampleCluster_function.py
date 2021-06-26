from MetabolicDistance import *
import cPickle as pickle;
import sys

def sv(obj,fname):
    '''Pickles object <obj> as <fname>.pkl'''
    fo=open(fname+'.pkl','wb');
    pickle.dump(obj, fo, pickle.HIGHEST_PROTOCOL)
    fo.close();


def pathTo(rxnid,saveTo='Output/Reactions/'):
    '''If <rxnid> is given (not None or False), then generates results of findPath()
    including a dictionary of path to every reachable reaction starting from <rxnid>. 
    Returns the rxnnetwork instance generated and path data from findPath() as a tuple.
    If not, then only rxnnetwork instance is generated and returned. This instance is based on
    model from <modelFile> constrained with constraints in <constraintsFile>.
    
    If <saveTo> and <rxnid>, saves the tuple from findPath() in the <saveTo> folder
    with filename <rxnid>.'''
    S,rxns,metabolites,lb,ub,byproducts=loadVarsFromFilenames('Input/Smatrix_regular.npy',\
    'Input/reactions_regular.txt','Input/metabolites_regular.txt','Input/LB_regular.txt',\
    'Input/UB_regular.txt','Input/byproducts_regular.txt');
    Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts);
    t=findPaths(rxnid,Net);
    sv(t,saveTo+rxnid);
    
pathTo(sys.argv[1]);