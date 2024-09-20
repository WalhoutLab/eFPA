import numpy as np
import pickle
import re

###############################################################################
#########      HELPER FUNCTIONS        ########################################
###############################################################################

def readList(fname,isNumerical=False):
    '''Reads a list from a text file (any extension but .pkl) or a pickled list file.'''
    if fname[-3:]=='pkl':
        fo=open(fname, 'rb');
        L=pickle.load(fo);
        fo.close(); 
    else:
        L=[];
        fo=open(fname,'r');
        for ln in fo:
            L.append(re.sub('\\n','',ln));
    if isNumerical:
        for i in range(len(L)):
            L[i]=float(L[i]);
    return L
    
def readMatrix(fname,delimiter=None):
    '''Reads matrix from file with path <fname>. If a delimiter is not provided, it is assumed
    that the file is an npy file created by numpy.'''
    if delimiter is None:
        S=np.load(fname);
    else:
        S=np.loadtxt(fname,delimiter=delimiter);  
    return S


def addVal2DictList(D,key,val):
    '''
    Subroutine of connectRxns2Metabolites(). Adds str or numerical values to list in dict D.
    '''
    if key in D:
        if type(D[key])!=list:
            D[key]==[D[key]];
        if not set(D[key]).issuperset(set([val])):
            D[key].append(val);
    else:
        D[key]=[val];
        
def fillInMissingItems(D1,D2):
    '''
    Subroutine of connectRxns2Metabolites(). Makes two given dictionaries the same size 
    by filling in missing keys as D[key]=[].'''
    U=set(D1).union(D2);        
    for x in U:
        if not x in D1:
            D1[x]=[];
        if not x in D2:
            D2[x]=[];
            

###############################################################################
##################      CLASSES        ########################################
###############################################################################


class simplenode():
    def __init__(self,rxnid,parent,children,path):
        '''
        Represents a reaction (self.id) with the parent (self.parent) and children (self.children). 
        Also inherits a path, to which self.id is added.
        
        USAGE
        
        n=simplenode(id,parent,children,path); 
        
        where, first entry is string, parent is node or string, children is a list of strings or nodes, 
        path is a list of strings.
        '''
        self.id=rxnid;
        self.path=[x for x in path];
        self.path.append(rxnid);
        self.parent=parent;
        self.children=[x for x in children];
        self.childno=len(self.children);


class rxnnode():
    def __init__(self,rxnid,fromThese,toThese):
        '''
        A node class that describes a reaction (self.id) based on reactions that produce its reactants
        (self.fromThese) and reactions that consume its products (self.toThese). 
        
        Methods exist for adding and removing parents (fromThese elements) and children (toThese elements)
        one at a time. 
        
        USAGE
        
        n=rxnnode(rxnid,fromThese,toThese);
        
        where, rxnid is a string and the other two entries are lists of strings or nodes.
        
        '''
        self.id=rxnid;
        self.fromThese=fromThese;
        self.toThese=toThese;
        self.fromNo=len(self.fromThese);  
        self.toNo=len(self.toThese);
        self.fromThese_ids=[x for x in self.fromThese];
        self.toThese_ids=[x for x in self.toThese];

    def addParent(self,rxn):
        if type(rxn)==str:
            rxnid=rxn;
        else:
            rxnid=rxn.id;
        self.fromThese.append(rxn);
        self.fromThese_ids.append(rxnid);
        self.fromNo=len(self.fromThese);

    def addChild(self,rxn):
        if type(rxn)==str:
            rxnid=rxn;
        else:
            rxnid=rxn.id;
        self.toThese.append(rxn);
        self.toThese_ids.append(rxnid);
        self.toNo=len(self.toThese);        
            
    def removeParent(self,rxnid):
        if self.fromThese_ids.count(rxnid):
            L=[x for x in self.fromThese];
            self.fromThese=[];
            self.fromThese_ids=[];
            for x in L:
                if not type(x)==str:
                    xid=x.id;
                else:
                    xid=x;
                if not rxnid==xid:
                    self.fromThese_ids.append(xid);
                    self.fromThese.append(x);
            self.fromNo=len(self.fromThese);
        else:
            print('Given rxn id is not found in parents. No action done.')
            
    def removeChild(self,rxnid):
        if self.toThese_ids.count(rxnid):
            L=[x for x in self.toThese];
            self.toThese=[];
            self.toThese_ids=[];
            for x in L:
                if not type(x)==str:
                    xid=x.id;
                else:
                    xid=x;
                if not rxnid==xid:
                    self.toThese_ids.append(xid);
                    self.toThese.append(x);
            self.toNo=len(self.toThese);
        else:
            print('Given rxn id is not found in children. No action done.')
        
    def getFroms(self,quiet=False):
        return self._fgetFromsOrTos(self.fromThese,quiet)        
    
    def getTos(self,quiet=False):
        return self._fgetFromsOrTos(self.toThese,quiet)
    
    def _fgetFromsOrTos(self,L,quiet):
        n=len(L);
        out=n*[''];
        for i in range(n):
            if type(L[i])==str:
                out[i]=(L[i],'not a node');
            else:
                out[i]=(L[i].id,'node');
        if not quiet:
            print(out)        
        return out


class rxnnetwork():
    toggle={'f':'r','r':'f'};
    maxk=100;
    def __init__(self,Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct):
        '''
        Establishes a reaction network with reactions converted to nodes defined by the 
        rxnnode class. Input (see below) should define the mapping of reactions to metabolites
        and back. Each reaction ID in the input must end with "f" or "r", indicating foreward
        and reverse directions, respectively. So, reversible reactions must be converted to 
        two unidirectional reactions a priori.
        
        INPUT
        
        Drxn2reactants:     A dictionary that maps reaction ID to a list of reactants.
        
        Drxn2products:      A dictionary that maps reaction ID to a list of products.
        
        Dcpd2asreactant:    A dictionary that maps each metabolite to a list of reactions where it is used as a reactant.
        
        Dcpd2asproduct:     A dictionary that maps each metabolite to a list of reactions where it is used as a product.
        
        OUTPUT
        
        An instance of rxnnetwork that has a dictionary mappint reaction IDs to reaction nodes. Various
        methods work with this map to traverse the network.
  
        USAGE
        
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        
        To traverse the network starting from a particular reaction:
        
            Net.pathForward(rxn ID);
        
        Then, Net.Dfirst_forward maps reactions that are reachable from the input reaction
        to the shortest path going in forward direction. This path may not be unique.
        '''
        self.Drxn2node={};
        for rxn in Drxn2products:
            self.Drxn2node[rxn]=rxnnode(rxn,\
            self.mergeLists([Dcpd2asproduct[cpd] for cpd in Drxn2reactants[rxn]],rxn),\
            self.mergeLists([Dcpd2asreactant[cpd] for cpd in Drxn2products[rxn]],rxn));
        for rxn in self.Drxn2node:
            self.Drxn2node[rxn].fromThese=[self.Drxn2node[self.Drxn2node[rxn].fromThese[i]] for i in range(self.Drxn2node[rxn].fromNo)];
            self.Drxn2node[rxn].toThese=[self.Drxn2node[self.Drxn2node[rxn].toThese[i]] for i in range(self.Drxn2node[rxn].toNo)];
                    
        self.removed={};
        
    def pathForward(self,rxn):
        '''Traverses the reaction network starting at the given reaction in forward direction, 
        until a terminal is reached, or self.maxk steps are taken, in every branching path. Thus,
        a tree (with reactions as nodes) is created such that the input reaction is the top node 
        and terminal reactions are terminal nodes.
        Only one path is defined and maintained for each reachable reaction, which is the one 
        used when the reaction is first encountered. 
        Nodes in the tree are defined based on the simplenode class.
        
        OUTPUT
        
        self.topnode_forward: simplenode instance representing the given reaction, where all
        paths start. This is the top node in the tree created.
        
        self.Dfirst_forward: A dictionary that maps reaction ID to a tuple in the form of (node in the tree,
        distance to the reaction, path to the reaction). The keys of this dictionary
        define all reactions that are reachable from the original node.
        
        self.terminals_forward: List of terminal nodes as reaction IDs.
  
        USAGE
        
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.pathForward(rxn ID);
        
        '''
        
        reverseid=self.toggleID(rxn);
        if reverseid in self.Drxn2node:
            self.removeRxn(reverseid,False);
            replaceReverseRxn=True;
        else:
            replaceReverseRxn=False;
            
        terminals=[];
        k=0;
        topnode=simplenode(rxn,None,self.Drxn2node[rxn].toThese_ids,[]);
        Dfirst={};
        self.updateDfirst(topnode,k,Dfirst);
        Npre=-1;
        while len(Dfirst)!=Npre and k<self.maxk:
            k=k+1;
            #print k
            Npre=len(Dfirst);
            newparents=[];
            for x in Dfirst:
                if Dfirst[x][1]==k-1:
                    newparents.append(Dfirst[x][0]);
            for pa in newparents:
                #print 'parent '+pa.id
                for i in range(pa.childno):
                    childid=pa.children[i];
                    #print childid;
                    grandchildren=self.Drxn2node[childid].toThese_ids;
                    pa.children[i]=simplenode(childid,pa,grandchildren,pa.path);
                    if not grandchildren:
                        terminals.append(pa.children[i].id); 
                    self.updateDfirst(pa.children[i],k,Dfirst);
            #print len(Dfirst);
        
        self.topnode_forward=topnode;
        self.Dfirst_forward=Dfirst;
        self.terminals_forward=terminals;
        
        if replaceReverseRxn:
            self.replaceRemoved([reverseid]);

    def findForwardLoops(self,quiet=True):
        '''Determines loops in the current forward tree (self.Dfirst_forward)
        created by self.pathForward(). A loop is a case when a path from top node
        to a reaction includes a reversible reaction in both ways (i.e., as "f" and 
        "r" forms). Reactions whose paths have loops are stored in a dictionary 
        (self.forwardLoops) that maps reaction ID to a node in self.Dfirst_forward.
        
        USAGE: 
        
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.pathForward(rxn ID);
        
        Net.findForwardLoops(quiet=True);
        OR
        DforwardLoops=Net.findForwardLoops(quiet=False);
        '''
        if not hasattr(self,'Dfirst_forward'):
            raise Exception('You must create a forward tree using self.pathForward() function before using this function.');
        self.ForwardLoops={};
        reverseid=self.toggleID(self.topnode_forward.id);
        for rxn in self.Dfirst_forward:
            if not rxn==reverseid:
                if self.checkLoop_forward(rxn):
                    self.ForwardLoops[rxn]=self.Dfirst_forward[rxn][2];
        if not quiet:
            return self.ForwardLoops

    def pathBackward(self,rxn):
        '''Same as pathForward() except that the paths found are the shortest paths from
        the given reaction (<rxn>) to others.
        '''        
        reverseid=self.toggleID(rxn);
        if reverseid in self.Drxn2node:
            self.removeRxn(reverseid,False);
            replaceReverseRxn=True;
        else:
            replaceReverseRxn=False;
 
        terminals=[];
        k=0;
        topnode=simplenode(rxn,None,self.Drxn2node[rxn].fromThese_ids,[]);
        Dfirst={};
        self.updateDfirst(topnode,k,Dfirst);
        Npre=-1;
        while len(Dfirst)!=Npre and k<self.maxk:
            k=k+1;
            #print k
            Npre=len(Dfirst);
            newparents=[];
            for x in Dfirst:
                if Dfirst[x][1]==k-1:
                    newparents.append(Dfirst[x][0]);
            for pa in newparents:
                #print 'parent '+pa.id
                for i in range(pa.childno):
                    childid=pa.children[i];
                    #print childid;
                    grandchildren=self.Drxn2node[childid].fromThese_ids;
                    pa.children[i]=simplenode(childid,pa,grandchildren,pa.path);
                    if not grandchildren:
                        terminals.append(pa.children[i].id); 
                    self.updateDfirst(pa.children[i],k,Dfirst);
            #print len(Dfirst);
        
        self.topnode_backward=topnode;
        self.Dfirst_backward=Dfirst;
        self.terminals_backward=terminals;
        
        if replaceReverseRxn:
            self.replaceRemoved([reverseid]);

    def findBackwardLoops(self,quiet=True):
        '''Determines loops in the current backward tree (self.Dfirst_backward)
        created by self.pathBackward(). A loop is a case when a path to top node
        from a reaction includes a reversible reaction in both ways (i.e., as "f" and 
        "r" forms). Reactions whose paths have loops are stored in a dictionary 
        (self.backwardLoops) that maps reaction ID to a node in self.Dfirst_backward.
        
        USAGE: 
        
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.pathBackward(rxn ID);
        
        Net.findBackwardLoops(quiet=True);
        OR
        DforwardLoops=Net.findBackwardLoops(quiet=False);
        '''
        if not hasattr(self,'Dfirst_backward'):
            raise Exception('You must create a backward tree using self.pathBackward() function before using this function.');
        self.BackwardLoops={};
        reverseid=self.toggleID(self.topnode_backward.id);
        for rxn in self.Dfirst_backward:
            if not rxn==reverseid:
                if self.checkLoop_backward(rxn):
                    self.BackwardLoops[rxn]=self.Dfirst_backward[rxn][2];
        if not quiet:
            return self.BackwardLoops


    def updateDfirst(self,node,k,Dfirst):
        '''Subroutine of tree function (e.g., self.pathForward()).'''
        if not node.id in Dfirst:
            Dfirst[node.id]=(node,k,node.path);

    def mergeLists(self,LL,rxnid):
        '''Subroutine of main class.'''
        S=set();
        for L in LL:
            S=S.union(L);
        return list(S.difference([self.toggleID(rxnid)]))

    def set_maxk(self,maxk):
        '''Sets maxk, the maximum number of steps to advance from top node of reaction tree
        to a terminal. Reactions that cannot be reached within maxk steps are excluded from tree.
        If this parameter is not set, default value of self.maxk is used.'''
        self.maxk=maxk;
        
    def getPaths(self,direction='forward'):
        '''Returns a dictionary of paths in the current reaction tree. It maps reaction ID 
        to the path from top node to the reaction.'''
        Dpath={};
        for k in getattr(self,'Dfirst_'+direction):
            Dpath[k]=getattr(self,'Dfirst_'+direction)[k][2];
        return Dpath

    def removeRxn(self,rxnid,replaceOld=True):
        '''Removes a reaction from the network. Optionally, previously removed reactions
        can be replaced first (if applicable). The removed reaction is added to 
        self.removed dictionary, which maps reaction ID to a reaction node that is
        extracted from the network. 
        
        USAGE:
        
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.removeRxn(reaction ID, replaceOld=True);
        '''
        if replaceOld:
            self.replaceRemoved();
        self.removed[rxnid]=self.Drxn2node[rxnid];
        for x in self.removed[rxnid].fromThese:
            x.removeChild(rxnid);
        for x in self.Drxn2node[rxnid].toThese:
            x.removeParent(rxnid);        
        
    def replaceRemoved(self,Lrxn=None):
        '''Brings all or a subset of previously removed reactions (in self.removed) back 
        to the reaction network.  
        
        USAGE:
        
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.removeRxn(reaction ID, replaceOld=True);
        
        Net.removeRxn(reaction IDn, Lrxn=None); #Replaces all reactions in self.removed
        OR
        Net.removeRxn(reaction ID, Lrxn=[reaction ID1, reaction ID2,...]); #Order of reactions in the list does not matter.
        '''
        if Lrxn is None:
            Lrxn=[x for x in self.removed];
        for rxn in Lrxn:
            self.Drxn2node[rxn]=self.removed.pop(rxn); 
            for x in self.Drxn2node[rxn].fromThese:
                x.addChild(self.Drxn2node[rxn]);
            for x in self.Drxn2node[rxn].toThese:
                x.addParent(self.Drxn2node[rxn]);            
        
    def toggleID(self,name):
        '''Toggles the forward-reverse indicator ("f" and "r", respectively) in reaction ID.'''
        return name[:-1]+self.toggle[name[-1]]
        
    def checkLoop_forward(self,rxn):
        '''Checks if the path to a given reaction has loop in it. See also findForwardLoops().
        
        USAGE
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.pathForward(rxn ID);
        boolean=Net.checkLoop_forward(another rxn ID);        
        '''
        S=set([y[:-1] for y in self.Dfirst_forward[rxn][2]]);
        if not len(S)==len(self.Dfirst_forward[rxn][2]):
            return True
        else:
            return False

    def checkLoop_backward(self,rxn):
        '''Same as checkLoop_forward() for the backward direction.        
        '''
        S=set([y[:-1] for y in self.Dfirst_backward[rxn][2]]);
        if not len(S)==len(self.Dfirst_backward[rxn][2]):
            return True
        else:
            return False


    def findFirstLoop_forward(self,rxn):
        '''Finds the first reversible reaction in a path whose reversed version exists downstream
        in the same path. The path is defined from the top node in the current tree (self.Dfirst_forward)
        towards a given reaction.
        
        USAGE
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.pathForward(rxn ID);
        x=Net.findFirstLoop(another rxn ID); #If there is no loop in the path to the input reaction, None is returned.       
        
        '''
        if self.checkLoop_forward(rxn):
            Lrxnroot=[x[:-2] for x in self.Dfirst_forward[rxn][2]];
            i=0;
            while i<len(Lrxnroot)-1:
                rxnroot=Lrxnroot[i];
                if set(Lrxnroot[i+1:]).issuperset([rxnroot]):
                    rxnid=self.Dfirst_forward[rxn][2][i];
                    return rxnid;
                i=i+1;
        else:
            return None

    def findFirstLoop_backward(self,rxn):
        '''
        Finds the first reversible reaction in a backward path whose reversed version exists upstream
        in the same path. The path is defined from the top node in the current tree (self.Dfirst_backward)
        towards a given reaction.
        
        USAGE
        Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
        Net.pathBackward(rxn ID);
        x=Net.findFirstLoop_backward(another rxn ID); #If there is no loop in the path to the input reaction, None is returned.       
        
        '''
        if self.checkLoop_backward(rxn):
            Lrxnroot=[x[:-2] for x in self.Dfirst_backward[rxn][2]];
            i=0;
            while i<len(Lrxnroot)-1:
                rxnroot=Lrxnroot[i];
                if set(Lrxnroot[i+1:]).issuperset([rxnroot]):
                    rxnid=self.Dfirst_backward[rxn][2][i];
                    return rxnid;
                i=i+1;
        else:
            return None
            
class dist_data:
    '''This class is designed to efficiently use distance data in hdf format.'''    
   
    def __init__(self, filePath):
        '''
        <filePath> is the path to the distance h5ad file
        '''
        import h5py
        
        self.hdf_file = h5py.File(filePath, 'r')
        rows=[x.decode() for x in self.hdf_file['rows'][:]]
        cols=[x.decode() for x in self.hdf_file['cols'][:]]
        self.rows=['_'.join([x[:-1],x[-1]]) for x in rows]
        self.cols=['_'.join([x[:-1],x[-1]]) for x in cols]
        self.row_index={self.rows[i]:i for i in range(len(self.rows))}
        self.col_index={self.cols[i]:i for i in range(len(self.cols))}
        self.rowno,self.colno=len(self.rows),len(self.cols)

    def getRxns(self,Lrxn):    
        '''Creates and returns a reaction distance matrix for the reactions in the list <Lrxn>. The output
        matrix is an instant of matrix class in myCobra'''
        from myCobra import matrix
        
        if not isinstance(Lrxn,(tuple,list)):
            Lrxn=[Lrxn]
        rxns_windex=[]
        for rxn in Lrxn:
            for Dir in ('f','r'):
                rxn_dir='_'.join([rxn,Dir])
                if rxn_dir in self.row_index:
                    rxns_windex.append((self.row_index[rxn_dir],rxn_dir))
        rxns_windex_sorted = sorted(rxns_windex, key=lambda x: x[0])
        I = [pair[0] for pair in rxns_windex_sorted] 
        rxns = [pair[1] for pair in rxns_windex_sorted] 
        
        vals_from = self.hdf_file['data'][I, :]
        vals_to=self.hdf_file['data'][:, I]
        vals_to_reorder = vals_from.T
        vals_min = np.minimum(vals_to, vals_to_reorder)
        
        rows,cols=[x for x in self.rows],rxns
        m=matrix(vals_min,rows,cols)
        return m

    def close(self):
        self.hdf_file.close()    
        
###############################################################################
##################      FUNCTIONS      ########################################
###############################################################################        

def connectRxns2Metabolites(S,rxns,metabolites,lb,ub,byproducts=[]):
    '''
    Matches reactions and metabolites of a metabolic network model using four 
    dictionaries (see OUTPUT). Reactions are first converted to forward and (if applicaple) 
    reverse forms indicated by "f" and "r" in the end of reaciton ID, respectively. 
    Then reactants and products in the unidirectional reactions are parsed to build the maps. 
    
    INPUT
    
    S: S-matrix of the metabolic network model
    
    rxns: List of reactions (IDs) of the model in the order of appearance in the S matrix.
    
    metabolites: List of compounds (IDs) of the model in the order of appearance in the S matrix.
    
    lb and ub: List of lower (lb) and upper (ub) boundaries of reactions in rxns.
    
    byproducts: List of metabolites to be ignored.
    
    
    OUTPUT
    
    Drxn2reactants:     A dictionary that maps reaction ID to a list of reactants.
    
    Drxn2products:      A dictionary that maps reaction ID to a list of products.
    
    Dcpd2asreactant:    A dictionary that maps each metabolite to a list of reactions where it is used as a reactant.
    
    Dcpd2asproduct:     A dictionary that maps each metabolite to a list of reactions where it is used as a product.
    
    
    USAGE
    
    Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct=connectRxns2Metabolites(S,rxns,metabolites,lb,ub,byproducts=[]);
    
    '''
    #rxnid='RM01608f';
    #Given, S,rxns,metabolites,lb,ub,and byproducts
    n=len(rxns);
    m=len(metabolites);
    
    if byproducts:
        Sby=set(byproducts);
        for i in range(m):
            if Sby.issuperset([metabolites[i]]):
                for j in range(n):
                    S[i,j]=0.;
        
    Dcpd2asreactant={};
    Dcpd2asproduct={}; 
    Drxn2reactants={};
    Drxn2products={};
    rdirs=[];
    for i in range(n):
        if lb[i]>=0 and ub[i]>0:
            rdirs.append(1);
        elif lb[i]<0 and ub[i]<=0:
            rdirs.append(-1);
        elif lb[i]<0 and ub[i]>0:
            rdirs.append(0);
        elif lb==0 and ub==0:
            rdirs.append(None);
        else:
            rdirs.append('ERROR')
            print('WARNING: cannot resolve direction of '+rxns[i]+' with the boundaries information');

    for i in range(n):
        rxn=rxns[i];
        rdir=rdirs[i];
        if not rdir is None:
            for j in range(m):
                cpd=metabolites[j];
                coef=S[j,i];
                if coef:
                    if coef<0:
                        if rdir==1 or rdir==0:
                            addVal2DictList(Drxn2reactants,rxn+'_f',cpd);
                            addVal2DictList(Dcpd2asreactant,cpd,rxn+'_f');
                        if rdir==-1 or rdir==0:
                            addVal2DictList(Drxn2products,rxn+'_r',cpd);
                            addVal2DictList(Dcpd2asproduct,cpd,rxn+'_r'); 
                    else: #coef>0
                        if rdir==-1 or rdir==0:
                            addVal2DictList(Drxn2reactants,rxn+'_r',cpd);
                            addVal2DictList(Dcpd2asreactant,cpd,rxn+'_r');
                        if rdir==1 or rdir==0:
                            addVal2DictList(Drxn2products,rxn+'_f',cpd);
                            addVal2DictList(Dcpd2asproduct,cpd,rxn+'_f'); 
    fillInMissingItems(Drxn2reactants,Drxn2products);
    fillInMissingItems(Dcpd2asreactant,Dcpd2asproduct);
    
    return Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct


def establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts=[]):
    '''
    Builds a reaction network with nodes as unidirectional reactions (i.e., reversible reacitons
    are represented by two separate irreversible reactions) based on the metabolic network
    model and reaction boundaries provided. Due to conversion to unidirectional model,
    reaction names in the network have "f" or "r" in the end, which indicate that the reaction
    represents the forward or reverse direction of the original model reaction, respectively.  
    INPUT
    
    S: S-matrix of the metabolic network model
    
    rxns: List of reactions (IDs) of the model in the order of appearance in the S matrix.
    
    metabolites: List of compounds (IDs) of the model in the order of appearance in the S matrix.
    
    lb and ub: List of lower (lb) and upper (ub) boundaries of reactions in rxns.
    
    byproducts: List of metabolites to be ignored.
    
    USAGE
    
    Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts=[]);
    '''
    Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct=connectRxns2Metabolites(S,rxns,metabolites,lb,ub,byproducts);
    Net=rxnnetwork(Drxn2reactants,Drxn2products,Dcpd2asreactant,Dcpd2asproduct);
    return Net

def findPaths(rxnid,Net,quiet=True):
    '''
    Finds shortest paths from a given reaction to others in a given reaction network 
    from the rxnnetwork class. The rxnnetwork class can already define shortest paths, but some of these
    paths have loops (the use of at least one reversible reaction in both directions). This
    function coorects those loops by removing key reactions from analysis and reports
    the list of reactions for which the found path may not be the shortest as exhaustive
    corrections are not made.
    USAGE:
    
    Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts=[]); #See establishRxnNetwork() 
    Dpath,badrxns,errors=findPaths(reaction ID,Net);
    
    where, Dpath is the dictionary from reactions to path (from reaction ID to reaction), badrxns
    is the list of reactions for which the found path is not guaranteed to be the shortest, and 
    errors is the list of reactions for which no calculation was made because of errors. 
    
    Other options:
        -quiet: set this as False so that the loop fixes are printed on the screen.
    '''
    #STATIC
    kmax=1000;
    warningThreshold=25;
    
    Net.pathForward(rxnid);
    Net.findForwardLoops();
    Dpath=Net.getPaths('forward');
    Llooprxns=list(Net.ForwardLoops);
    if not quiet:
        N=len(Llooprxns);
        if N>warningThreshold:
            print('More than '+str(warningThreshold)+' loops were detected. Loop correction may take a few minutes.')
    Lbadrxns=[];
    Lerror=[];
    for looprxn in Llooprxns:
        isLoop=Net.checkLoop_forward(looprxn);
        path=Net.Dfirst_forward[looprxn][2];
        k=0;
        while path and isLoop and k<kmax:
            k=k+1;
            if not quiet:
                print('Reaction '+looprxn+' has a loop.')
            loopid=Net.findFirstLoop_forward(looprxn);
            Net.removeRxn(loopid,False);
            if not quiet:
                print('The following reactions are eliminated to fix the loop:'+str(list(Net.removed)));
            Net.pathForward(rxnid);
            if looprxn in Net.Dfirst_forward:
                path=Net.Dfirst_forward[looprxn][2];
            else:
                path=None;
            if path:
                isLoop=Net.checkLoop_forward(looprxn)
        if k<kmax:
            Dpath[looprxn]=path;
            if not len(Net.removed)==1 or not list(Net.removed)[0]==Net.toggleID(looprxn):
                Lbadrxns.append(looprxn);
        else:
            Lerror.append(looprxn);
        Net.replaceRemoved();
        Net.pathForward(rxnid);
    return Dpath,Lbadrxns,Lerror
    
def findPaths_backward(rxnid,Net,quiet=True):
    '''
    Same as findPaths() function except that it works on backward paths and loops.
    '''
    #STATIC
    kmax=1000;
    warningThreshold=25;
    
    Net.pathBackward(rxnid);
    Net.findBackwardLoops();
    Dpath=Net.getPaths('backward');
    Llooprxns=list(Net.BackwardLoops);
    if not quiet:
        N=len(Llooprxns);
        if N>warningThreshold:
            print('More than '+str(warningThreshold)+' loops were detected. Loop correction may take a few minutes.')
    Lbadrxns=[];
    Lerror=[];
    for looprxn in Llooprxns:
        isLoop=Net.checkLoop_backward(looprxn);
        path=Net.Dfirst_backward[looprxn][2];
        k=0;
        while path and isLoop and k<kmax:
            k=k+1;
            if not quiet:
                print('Reaction '+looprxn+' has a loop.')
            loopid=Net.findFirstLoop_backward(looprxn);
            Net.removeRxn(loopid,False);
            if not quiet:
                print('The following reactions are eliminated to fix the loop:'+str(list(Net.removed)));
            Net.pathBackward(rxnid);
            if looprxn in Net.Dfirst_backward:
                path=Net.Dfirst_backward[looprxn][2];
            else:
                path=None;
            if path:
                isLoop=Net.checkLoop_backward(looprxn)
        if k<kmax:
            Dpath[looprxn]=path;
            if not len(Net.removed)==1 or not list(Net.removed)[0]==Net.toggleID(looprxn):
                Lbadrxns.append(looprxn);
        else:
            Lerror.append(looprxn);
        Net.replaceRemoved();
        Net.pathBackward(rxnid);
    return Dpath,Lbadrxns,Lerror
    
def convertPaths2Distances(Dpath):
    '''
    Converts a path dictionary from findPaths() to a distance dictionary. 
    
    USAGE:
    
    Ddistance=convertPaths2Distances(Dpath);
    '''
    Ddistance={};
    for rxnid in Dpath:
        if not Dpath[rxnid] is None:
            Ddistance[rxnid]=len(Dpath[rxnid])-1;
    return Ddistance
    
def findDistances(rxnid,Net):
    '''
    Uses findPaths() function to calculate the shortest paths from a given reaction to others 
    in a given reaction network from the rxnnetwork class. Then uses convertPaths2Distances()
    function to convert these shortest paths to minimum distance.
    USAGE:
    
    Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts=[]); #See establishRxnNetwork() 
    Ddistance=findPaths(reaction ID,Net); #See establishRxnNetwork() 
    
    where, Ddistance is the dictionary from reactions to distance (from reaction ID to reaction).
    The last letter of Reaciton ID must indicate the direction ("f" or "r").
    '''  
    Dpath,_,_=findPaths(rxnid,Net);
    return convertPaths2Distances(Dpath)

    
def loadVarsFromFilenames(S,rxns,metabolites,lb,ub,byproducts,delimiterInSmatrix=None):
    '''
    Can be used to load all input variables of establishRxnNetwork() from filenames
    (the inputs of this function). The file for the S matrix must be a npy file generated by numpy 
    if <delimiterInSmatrix> is None. For text-formatted matrix files, indicate the delimiter using
    <delimiterInSmatrix). All other variables are one-dimensional lists. For each of these, 
    the file must be either a text file with all elements written line by line, 
    or a .pkl file of a pickled list.
    
    USAGE:
    
    Smatrix,reactions,metabolites,lb,ub,byproducts=loadVarsFromFilenames(fname_Smatrix, \
    fname_reactions,fname_metabolites,fname_lb,fname_ub,fname_byproducts)
    '''
    S=readMatrix(S,delimiter=delimiterInSmatrix);
    rxns=readList(rxns);
    metabolites=readList(metabolites);
    lb=readList(lb,isNumerical=True);
    ub=readList(ub,isNumerical=True);
    byproducts=readList(byproducts);
    return S,rxns,metabolites,lb,ub,byproducts   





    