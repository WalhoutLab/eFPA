import cobra
import pandas as pd
import numpy as np
import re

import myTable

from MetabolicDistance import *

import importlib
importlib.reload(myTable)

#from optlang.glpk_interface import Model, Variable, Constraint, Objective
from optlang import Model, Variable, Constraint, Objective
from optlang.symbolics import Basic, Zero

from copy import deepcopy

import time

try:
    import pickle5 as pickle
except:
    import pickle

### GLOBAL
pd.set_option('display.max_columns', None);
#Set Tolerances for the solver, using tols dictionary
#feasibility,optimality,integrality=1E-6,1E-6,1E-5; 
tols={'feasibility':1E-9,'optimality':1E-9,'integrality':1E-9}; #for Gurobi, including mixed integer programming
#tols={'feasibility':1E-6,'optimality':1E-6,'integrality':1E-5}; #default

### CLASSES
class matrix:
    '''Hashable 2D array object'''
    def __init__(self,array,rows,cols):
        self.vals=np.array(array);
        if not self.vals.size:
            if rows or cols:
                print('Empty array cannot have rows or columns, these will be converted to empty lists.');
                rows,cols=[],[];
        elif not len(rows)==self.vals.shape[0] or not len(cols)==self.vals.shape[1]:
            raise Exception('The shape of array not consistent with number of rows or columns provided');
        self.rows,self.cols=rows,cols;
        self.rowindex={rows[i]:i for i in range(len(self.rows))};
        self.colindex={cols[i]:i for i in range(len(self.cols))};
        
    def copy(self):
        return matrix(self.vals.copy(),self.rows.copy(),self.cols.copy());
        
    def getCols(self,cols):
        '''Returns matrix with only columns indicated in <cols>.'''
        I=[self.colindex[col] for col in cols]; 
        return matrix(self.vals[:,I],self.rows,cols)

    def getRows(self,rows):
        '''Returns matrix with only rows indicated in <rows>.'''
        I=[self.rowindex[row] for row in rows]; 
        return matrix(self.vals[I,:],rows,self.cols)        
        
    def appendRows(self,arr,newrows):
        '''Adds the np array <arr> rowwise, and expands row names by <newrows>. Returns the 
        modified matrix as a new object.'''
        Sint=set(newrows).intersection(self.rows);        
        if Sint:
            raise Exception('The following row names already exist: '+str(Sint));
        newvals=np.concatenate((self.vals,arr),axis=0);
        return matrix(newvals,self.rows+newrows,[x for x in self.cols])
        
    def appendCols(self,arr,newcols):
        '''Adds the np array <arr> columnwise, and expands column names by <newcols>. Returns the 
        modified matrix as a new object.'''
        Sint=set(newcols).intersection(self.cols);        
        if Sint:
            raise Exception('The following column names already exist: '+str(Sint));
        newvals=np.concatenate((self.vals,arr),axis=1);
        return matrix(newvals,[x for x in self.rows],self.cols+newcols)
        
    def getRow(self,row):
        '''Returns <row> as np.array'''
        try:
            return self.vals[self.rowindex[row],:];
        except:
            print('item not found')
            return None            

    def getCol(self,column):
        '''Returns <column> as np.array'''
        try:
            return self.vals[:,self.colindex[column]]
        except:
            print('item not found')
            return None              
        
    def isempty(self):
        if self.vals.size:
            return False
        else:
            return True
    
    def to_df(self):
        '''Returns a pandas dataframe representation of the matrix.'''
        return pd.DataFrame(self.vals,index=self.rows,columns=self.cols)

    @property
    def rowno(self):
        return len(self.rows)
        
    @property
    def colno(self):
        return len(self.cols)

    def __getitem__(self,ij):
        try:
            return self.vals[self.rowindex[ij[0]],self.colindex[ij[1]]]; 
        except:
            print('item not found')
            return None
            
    def __setitem__(self,ij,val):
        try:
            self.vals[self.rowindex[ij[0]],self.colindex[ij[1]]]=val;
        except:
            print('row ('+str(ij[0])+') or column ('+str(ij[1])+') not found')
            #print('row or column not found')  

class tensor3:
    '''Hashable 3D array object'''
    def __init__(self,array,I,J,K):
        self.vals=np.array(array);
        if not self.vals.size:
            if I or J or K:
                print('Empty array cannot have indices, all indices will be converted to empty lists.');
                I,J,K=[],[],[];
        elif not len(I)==self.vals.shape[0] or not len(J)==self.vals.shape[1] or not len(K)==self.vals.shape[2]:
            raise Exception('The shape of array not consistent with number of indices provided');
        self.I,self.J,self.K=I,J,K;
        self.Iindex={I[i]:i for i in range(len(self.I))};
        self.Jindex={J[i]:i for i in range(len(self.J))};
        self.Kindex={K[i]:i for i in range(len(self.K))};
    
    def getI(self,i):
        '''Returns element <i> in I as a J by K matrix.'''
        vals=self.vals[self.Iindex[i],:,:]; 
        return matrix(vals,self.J,self.K)

    def getJ(self,j):
        '''Returns element <j> in J as a I by K matrix.'''
        vals=self.vals[:,self.Jindex[j],:]; 
        return matrix(vals,self.I,self.K)

    def getK(self,k):
        '''Returns element <k> in K as a I by J matrix.'''
        vals=self.vals[:,:,self.Kindex[k]]; 
        return matrix(vals,self.I,self.J)

    @property
    def Ino(self):
        return len(self.I)
        
    @property
    def Jno(self):
        return len(self.J)

    @property
    def Kno(self):
        return len(self.K)
        
    def __getitem__(self,ijk):
        try:
            return self.vals[self.Iindex[ijk[0]],self.Jindex[ijk[1]],self.Kindex[ijk[2]]]; 
        except:
            print('item not found')
            return None
            
    def __setitem__(self,ijk,val):
        try:
            self.vals[self.Iindex[ijk[0]],self.Jindex[ijk[1]],self.Kindex[ijk[2]]]=val;
        except:
            print('I ('+str(ijk[0])+') or J ('+str(ijk[1])+') or K ('+str(ijk[2])+') not found')         
            

class gprnode:
    '''Node class for describing gene associations in gene-protein-reaction (GPR) rules.'''
    def __init__(self,index,parent,children,connector=None,gene=None):
        self.index=index;
        self.parent=parent;
        self.children=children;
        self.gene=gene;
        self.connector=connector;
        
    def appendChild(self,child):
        self.children.append(child);
        
    def setGene(self,gene):
        self.gene=gene;
    
    def setConnector(self,connector):
        self.connector=connector;

class opt():
    '''Interface for optimization using optlang.'''
    def __init__(self,rxns,modelname='myModel',irreversible=False):
        '''Establishes an optlang model with name <modelname> based on reactions 
        in <rxns>, a list of COBRA reactions. Reaction IDs and boundaries define
        the variables and their boundaries, respectively. Reaction stoichiometry
        defines the metabolite constraints, named after metabolite IDs.

        To make a model that has only irreversible reactions, set <irreversible>
        as True. 
        
        To set solver tolerances (i.e., to change the default values for feasibility,
        optimality, and integrality defined for this class), use setTolerances()
        function.
        '''
        self.Dcons={};
        self.Dvars={};
        self.rxn2varcoefs={};
        self.var2absvar={};
        self.var2binvars={};
        self.gene2intvars={};
        self.model=Model(modelname);
        if tols:
            self.setTolerances(tols);
        self.isIrreversible=irreversible;
        temp_cpd={};
        if not irreversible:
            for r in rxns:
                self.Dvars[r.id]=Variable(r.id, lb=r.lower_bound,ub=r.upper_bound);
                self.rxn2varcoefs[r.id]={r.id:1.};
                self.model.add(self.Dvars[r.id]);
                for met in r.metabolites:
                    name=met.id;
                    if not name in temp_cpd:
                        temp_cpd[name]={};
                    temp_cpd[name][self.Dvars[r.id]]=r.metabolites[met];
        else:
            for r in rxns:
                id_f=r.id+'_f';
                self.Dvars[id_f]=Variable(id_f, lb=max(0,r.lower_bound),ub=max(0,r.upper_bound));
                self.rxn2varcoefs[r.id]={id_f:1.};
                self.model.add(self.Dvars[id_f]);
                for met in r.metabolites:
                    name=met.id;
                    if not name in temp_cpd:
                        temp_cpd[name]={};
                    temp_cpd[name][self.Dvars[id_f]]=r.metabolites[met];
                if r.reversibility or r.lower_bound<0: #the latter is for <-- reactions that may be designated irreversible
                    id_r=r.id+'_r';
                    self.Dvars[id_r]=Variable(id_r, lb=max(0,-r.upper_bound),ub=max(0,-r.lower_bound));
                    self.rxn2varcoefs[r.id][id_r]=-1.;
                    self.model.add(self.Dvars[id_r]);
                    for met in r.metabolites:
                        name=met.id;
                        if not name in temp_cpd:
                            temp_cpd[name]={};
                        temp_cpd[name][self.Dvars[id_r]]=-r.metabolites[met];
        for met in temp_cpd:
            self.model.add(Constraint(Zero,lb=0.,ub=0.,name=met));
        for con in self.model.constraints:
            con.set_linear_coefficients(temp_cpd[con.name]);
            self.Dcons[con.name]=con;
            
    def setTolerances(self,tols):
        '''Sets model configuration tolerances according to <tols> dictionary 
        (tolerance -> value). Typical keys are feasibility, optimality, and integrality.'''
        for tol in tols:
            try:
                setattr(self.model.configuration.tolerances,tol,tols[tol]);
            except:
                print('Tolerance '+tol+' could not be set for this solver');
            
    def setObjective(self,objf,direction='max'):
        '''Sets the model objective using var -> coefficient pairs in <objf> dictionary.
        Any variable in the model can be used. Use <direction> to indicate if the
        objective is to maximize (direction="max") or minimize ("min") the objective
        function.'''
        self.model.objective=Objective(0,direction=direction);
        self.model.objective.set_linear_coefficients({self.Dvars[k]:objf[k] for k in objf});
        self.model.objective.direction=direction;  
        
    def addConstraint(self,D,name,lb,ub):
        '''Adds a constraint named <name> based on variable_name -> coefficient pairs in <D>, 
        with lower and upper boundaries <lb> and <ub>, respectively. Variables and coefficients
        in <D> define a linear expression.'''
        self.model.add(Constraint(Zero,lb=lb,ub=ub,name=name));
        self.Dcons[name]=getattr(self.model.constraints,name);
        self.Dcons[name].set_linear_coefficients({self.Dvars[k]:D[k] for k in D}); 
        
    def addAbsVar(self,varname):
        '''Adds an absolute value variable for the existing variable <varname>. The
        name of the absolute variable will be "<varname>_abs". 
        
        The constraints associated with this variable will be named 
        "<varname>_abs_1" and "<varname>_abs_2"'''
        if not varname in self.Dvars:
            raise Exception('No variable named '+varname+' to make an auxiliary to define the absolute function for!');
        absvarname=varname+'_abs';
        conname1,conname2=absvarname+'_1',absvarname+'_2';
        self.Dvars[absvarname]=Variable(absvarname, lb=0.);
        self.model.add(self.Dvars[absvarname]);
        self.model.add(Constraint(Zero,lb=0.,name=conname1));
        self.model.add(Constraint(Zero,lb=0.,name=conname2));
        self.Dcons[conname1]=getattr(self.model.constraints,conname1);
        self.Dcons[conname2]=getattr(self.model.constraints,conname2);
        self.Dcons[conname1].set_linear_coefficients({self.Dvars[varname]:-1.,self.Dvars[absvarname]:1.});
        self.Dcons[conname2].set_linear_coefficients({self.Dvars[varname]:1.,self.Dvars[absvarname]:1.});
        self.var2absvar[varname]=absvarname; 
        
    def addMomaVar(self,rxn,refFlux):
        '''Adds a MOMA absolute variable for <rxn> and its reference flux <refFlux>. If the MOMA variable
        for the given reaction is already set, then only boundaries will be changed based on <refFlux>.'''
        if not rxn in self.Dvars:
            raise Exception('Reaction not among current optimizer model variables!');
        varname='_'.join([rxn,'moma']);
        conname1,conname2=''.join([varname,'_1']),''.join([varname,'_2']);
        if not varname in self.Dvars:
            self.Dvars[varname]=Variable(varname, lb=0.);
            self.model.add(self.Dvars[varname]);
            self.model.add(Constraint(Zero,lb=-refFlux,name=conname1));
            self.model.add(Constraint(Zero,lb=refFlux,name=conname2));           
            self.Dcons[conname1]=getattr(self.model.constraints,conname1);
            self.Dcons[conname2]=getattr(self.model.constraints,conname2);
            self.Dcons[conname1].set_linear_coefficients({self.Dvars[rxn]:-1.,self.Dvars[varname]:1.});
            self.Dcons[conname2].set_linear_coefficients({self.Dvars[rxn]:1.,self.Dvars[varname]:1.});
        else:
            self.Dcons[conname1].lb=-refFlux;
            self.Dcons[conname2].lb=refFlux;
        return varname
        
        
    def addOnRxn(self,rxn,eps):
        '''Adds IMAT-type binary variables and relevant constraints to force 
        reaction <rxn> to carry flux of at least the possible values indicated by <eps>, 
        which is a tuple of flux thresholds in both directions of the reaction
        (eps_reverse,eps_forward).
        
        The added IMAT variable(s) will be named "<rxn>_yf" for forward and "<rxn>_yr"
        for reverse directions, when either is applicable.'''
        if self.Dvars[rxn].ub>0:        
            varname=rxn+'_yf';
            if not varname in self.Dvars:
                conname=varname;
                self.Dvars[varname]=Variable(varname,type='binary');
                self.model.add(self.Dvars[varname]);
                self.model.add(Constraint(Zero,lb=self.Dvars[rxn].lb,name=conname));
                self.Dcons[conname]=getattr(self.model.constraints,conname);
                self.Dcons[conname].set_linear_coefficients({self.Dvars[rxn]:1.,self.Dvars[varname]:self.Dvars[rxn].lb-eps[1]});
                if rxn in self.var2binvars:
                    self.var2binvars[rxn].append(varname);
                else:
                    self.var2binvars[rxn]=[varname];
        if self.Dvars[rxn].lb<0:        
            varname=rxn+'_yr';
            if not varname in self.Dvars:
                conname=varname;
                self.Dvars[varname]=Variable(varname,type='binary');
                self.model.add(self.Dvars[varname]);
                self.model.add(Constraint(Zero,ub=self.Dvars[rxn].ub,name=conname));
                self.Dcons[conname]=getattr(self.model.constraints,conname);
                self.Dcons[conname].set_linear_coefficients({self.Dvars[rxn]:1.,self.Dvars[varname]:self.Dvars[rxn].ub+eps[0]});     
                if rxn in self.var2binvars:
                    self.var2binvars[rxn].append(varname);
                else:
                    self.var2binvars[rxn]=[varname];

    def addOffRxn(self,rxn):
        '''Adds an IMAT-type binary variable and relevant constraints to prevent reaction <rxn> from carrying 
        any flux.
        
        The binary variable will be named "<rxn>_y" and constraints will be named
        "<rxn>_y_1" and "<rxn>_y_2"'''       
        varname=rxn+'_y';
        if not varname in self.Dvars:
            conname1,conname2=varname+'_1',varname+'_2';
            self.Dvars[varname]=Variable(varname,type='binary');
            self.model.add(self.Dvars[varname]);
            self.model.add(Constraint(Zero,lb=self.Dvars[rxn].lb,name=conname1));
            self.model.add(Constraint(Zero,ub=self.Dvars[rxn].ub,name=conname2));
            self.Dcons[conname1]=getattr(self.model.constraints,conname1);
            self.Dcons[conname2]=getattr(self.model.constraints,conname2);
            self.Dcons[conname1].set_linear_coefficients({self.Dvars[rxn]:1.,self.Dvars[varname]:self.Dvars[rxn].lb});
            self.Dcons[conname2].set_linear_coefficients({self.Dvars[rxn]:1.,self.Dvars[varname]:self.Dvars[rxn].ub});
            if rxn in self.var2binvars:
                self.var2binvars[rxn].append(varname);
            else:
                self.var2binvars[rxn]=[varname];        

    def addOnGene(self,gene,rxns):
        '''Adds an IMAT++ binary variable and relevant constraints to force 
        gene <gene> to carry flux at least one associated reaction in <rxns>. 
        Assumes that binary rxn variables have already been enetered.
        
        Technically, <gene> only defines the variable and constraint names and
        <rxns> define the actual contraints.

        The variable will be named "<gene>_y", and the constraints "<gene>_y_1"
        and "<gene>_y_2".         
        '''
        varname=gene+'_y';
        if not varname in self.Dvars:
            conname1,conname2=varname+'_1',varname+'_2';     
            self.Dvars[varname]=Variable(varname,type='integer');
            self.model.add(self.Dvars[varname]);
            self.model.add(Constraint(Zero,ub=0,name=conname1));
            self.model.add(Constraint(Zero,ub=1,name=conname2));
            self.Dcons[conname1]=getattr(self.model.constraints,conname1);
            self.Dcons[conname2]=getattr(self.model.constraints,conname2);
            dcoef={self.Dvars[varname]:1.};
            for rxn in rxns:
                for var in self.var2binvars[rxn]:
                    dcoef[self.Dvars[var]]=-1.;
            self.Dcons[conname1].set_linear_coefficients(dcoef);
            self.Dcons[conname2].set_linear_coefficients({self.Dvars[varname]:1.});
            self.gene2intvars[gene]=varname;        
        
    def constrainObjective(self,tol,name):
        '''Constrains the current objective within <tol> from current value. 
        
        The constraint will be named as <name>, and any present constraint with 
        this name will be removed first.'''
        lincoef=self.model.objective.get_linear_coefficients(self.model.objective.variables);
        if name in self.Dcons:
            self.model.remove(self.Dcons[name]);
        if self.model.objective.direction=='max':
            self.model.add(Constraint(Zero,lb=self.model.objective.value-tol,name=name)); ##
        else:
            self.model.add(Constraint(Zero,ub=self.model.objective.value+tol,name=name));
        self.Dcons[name]=getattr(self.model.constraints,name);
        self.Dcons[name].set_linear_coefficients(lincoef);
        
    def optimize(self):
        '''Solve the current problem for the current objective'''
        return self.model.optimize();   

    
        

class mnm():
    sizeRandomChoiceForUniqueMets=10; #number of metabolites to test for regular annotation of compartmentalization
    tol=1E-8;    
    lowTol=1E-12;
    toggledir={'f':'r','r':'f'};
    Dw_elements={'C': 12.011, 'H': 1.0079, 'O': 15.999, 'N': 14.007, 'P': 30.974, 'S': 32.065, 'Na': 22.99, 'Cl': 35.453, \
    'Fe': 55.845, 'Se': 78.971, 'K': 39.098, 'I': 126.9, 'Mo': 95.95, 'Co': 58.933, 'Br': 79.904, 'Zn': 65.38, 'Ba': 137.327, \
    'Ca': 40.078, 'Cu': 63.546, 'F': 18.998403, 'Li': 6.941, 'Mg': 24.305, 'X': 1000.0, 'R': 100.0, 'Y': 1000.0};
    def __init__(self,model=None,filetype=None,verbose=True):
        '''
        Metabolic network modeling class that harbors a Cobra model as self.model.
        Use the <model> (default=None) and filetype (default=None) entries to define the
        model in one of the following ways: 
        
        1. <model> is a cobra model, <filetype> not important.
        2. <model> is a filename (path) of a file accepted by the loadModel() function and
           <filetype> indicates the type of file (see loadModel() function).
        3. <model> is the filename (path) of a previously saved workspace file (.pkl) and
           <filetype> is "workspace" (or "ws"). Extension fo filename will be ignored.
        4. <model> is null (<filename> ignored). In this case, use loadModel() or loadWorkspace()
           functions to incorporate a model into the mnm object.
           
        Set <verbose> as False to avoid printing warning messages.   
        '''
        #Read model
        if not filetype:
            if model:
                if isinstance(model,cobra.core.model.Model):
                    self.model=model;
                    self.update();
                else:
                    raise Exception('<model> entry not understood. Please follow the instructions for inputs.');
            else:
                self.model=None;
                if verbose:
                    print('No model is given, so the mnm object is empty. Use loadModel or loadWorkspace functions to initiate a model.');
        elif filetype=='workspace' or filetype=='ws':
            self.loadWorkspace(model);
        else:
            self.loadModel(model,filetype);
        
    def loadModel(self,filename,filetype):
        '''Loads cobra model from <filename> assuming the file is of type indicated by 
        <filetype>. The following file types are valid: 
                "sbml"
                "json"
        '''
        if filetype=='sbml':
            self.model=cobra.io.read_sbml_model(filename);
        elif filetype=='json':
            self.model=cobra.io.load_json_model(filename);
        else:
            raise Exception('<filetype> not understood. It must be either None (default) or one of the valid filetypes.');        
        self.update(); 

    def loadWorkspace(self,filename):
        '''Loads a previously saved workspace (cobra model, optimizer [if applicable], and solution [if applicable]) 
        from <filename>.'''
        ws=self.op(filename); #workspace saved as a dictionary
        self.model=cobra.io.from_json(ws['model']);
        self.update();
        if 'optimizer' in ws:
            self.__setOptimizerFromDict(ws['optimizer']);
        for att in ('sln','fluxes'):
            if att in ws:
                setattr(self,att,ws[att]);        
        
    def update(self,verbose=True):
        '''Updates mnm object with respect to changes in .model. Sets .sln to nothing.'''
        self.mets={x.id:x for x in self.model.metabolites};
        self.rxns={x.id:x for x in self.model.reactions};
        self.genes,self.genes_id={},{};
        for x in self.model.genes:
            self.genes_id[x.id]=x;
            if x.name:
                self.genes[x.name]=x;
            else:
                self.genes[x.id]=x;
        self.uniqmets=self.__getUniqueMetabolites();
        self.irrevIDs=self.__getIrrevRxnIDs();
        self.sln=None;
        if hasattr(self,'optimizer'):
            delattr(self,'optimizer');
            if verbose:
                print('The current optimizer is deleted since the model is updated. Set the optimizer again if needed.');

    def saveModel(self,filename,filetype='json'):
        '''Saves the cobra model as given filetype to given filename. 
        filteype can be one of "json","sbml", or "mat" (case insensitive). 
        '''
        filetype=filetype.lower();
        if filetype=='json':
            cobra.io.save_json_model(self.model,filename);
        elif filetype=='sbml':
            cobra.io.write_sbml_model(self.model,filename);
        elif filetype=='mat':
            cobra.io.save_matlab_model(self.model,filename);
        else:
            raise Exception('filetype not understood. Choose one of json","sbml", or "mat".');

    def saveWorkspace(self,filename):
        '''Saves the current workspace (cobra model, optimizer [if applicable], and solution [if applicable]) 
        into <filename>.'''
        ws={'model':cobra.io.to_json(self.model)}; #workspace to be saved as a dictionary
        if hasattr(self,'optimizer'):
            ws['optimizer']=self.__optimizer2dict();
        ws['sln']=self.sln;
        if hasattr(self,'fluxes'):
            ws['fluxes']=self.fluxes;
        self.sv(ws,filename);

    def setTolerance(self,tol):
        '''Sets the tolerance threshold to define zero.'''
        self.tol=tol;

    def setConstraints(self,Dconstraints,zeroOtherBoundary=False,verbose=True):
        '''Imposes lower and upper boundary constraints on the model and updates
        self accordingly. Note that the optimizer will be deleted if set. 
        
        Input:
        -----
        Dconstraints: A dictionary that maps reactions to constraints in the form
                      rxnid:[LB,UB]
         
        zeroOtherBoundary: Boolean that indicates if the lower boundary of 
        boundary reactions (exchanges and sinks) not included in Dconstraints 
        AND NOT -1000 (or whatever is used to represent negative inf) should be
        made zero.
        '''
        for rxn in Dconstraints:
            self.rxns[rxn].lower_bound,self.rxns[rxn].upper_bound=Dconstraints[rxn];
        if zeroOtherBoundary:
            if not hasattr(self,'infinity'):
                self.__getInfinity();
            for robj in self.model.boundary:
                if not robj.id in Dconstraints:
                    if robj.lower_bound and not robj.lower_bound==-self.infinity:
                        robj.lower_bound=0.;
        self.update();
        #if hasattr(self,'optimizer') and verbose:
        #    print('These constraints did not update the current optimizer. Set constraints before setting the optimizer to see the effect or use editing methods like addConstraint');

    def setFluxSln(self,fluxDict,na=0.):
        '''Sets the solution as a Cobra solution based on fluxes in flux dictionary
        <fluxDict>. Objective function is arbitrarily set at 0, and status is "fictional".
        Rxns that are not represented in the provided flux dictionary take the value 
        in <na>.'''
        rxns,fluxes=[],[];        
        for rxn in self.rxns:
            rxns.append(rxn);
            if rxn in fluxDict:
                fluxes.append(fluxDict[rxn]);
            else:
                fluxes.append(na);                
        self.sln=cobra.Solution(0,'fictional',pd.Series(fluxes,rxns),None,None);

    def getSmatrix(self):
        '''Generates a stoichiometry matrix with rows as metabolites, columns as reactions, and
        elements as metabolite coefficients. This matrix is returned as a "matrix" object
        (see matrix class).'''
        S=matrix(np.zeros((len(self.mets),len(self.rxns))),list(self.mets),list(self.rxns));
        for rxn in S.cols:
            for met in self.rxns[rxn].metabolites:
                S[met.id,rxn]=self.rxns[rxn].metabolites[met];
        return S     

    def getMW(self,met,update_weights={}):
        '''Returns molecular weight of indicated metabolite <met> based on atomic weights in <<Dw_elements>> and 
        <update_weights>. The latter updates the atomic weight dictionary <<Dw_elements>> and must provide assumed mass
        for pseudo elements such as R and X unless the default values are to be used.'''
        dw_elements=self.Dw_elements.copy()
        dw_elements.update(update_weights)
        mw=0
        for el in self.mets[met].elements:
            mw=mw+self.mets[met].elements[el]*dw_elements[el]
        return mw

    def getMetaboliteSummary(self,met,display=True):
        '''Returns fluxes for metabolite <met> (cobra metabolite summary table) 
        for the current solution. If <display>, the table is printed on screen.'''
        if not self.sln:
            raise Exception('No solutions exist. Do some optimization first.');
        Tbl=self.mets[met].summary(self.sln);
        if display:
            print(Tbl);
        return Tbl

    def getInputOutputSummary(self,display=True,):
        '''Returns fluxes of active input/output reactions (cobra model summary table) 
        for the current solution. If <display>, the table is printed on screen.'''
        if not self.sln:
            raise Exception('No solutions exist. Do some optimization first.');
        Tbl=self.model.summary(self.sln);
        if display:
            try:
                print(Tbl);
            except:
                print('Cannot display cobra table')
        return Tbl

    def getTotalFlux(self):
        '''Returns the sum of absolute fluxes of reactions.'''
        if hasattr(self,"sln"):
            s=np.sum(np.abs(self.sln.fluxes));
        else:
            s=0;
            for r in self.model.reactions:
                s=s+abs(r.flux);
        return s     

    def getPerturbedRxns(self,genes,coef=0,exp=None):
        '''Generates a dictionary of reactions to perturbation coefficients given
        <genes> (in IDs) to be perturbed. In subsequent applications, these 
        perturbation coefficients mutiplied with reaction fluxes would give the pertubed
        fluxes. 
        
        In the default case when gene expressions are not provided, perturbation coefficients
        will be 0 or 1, depending on the GPR associations of reactions with perturbed genes.
        If gene expression levels are provided, a quantitative assessment of reaction expression
        will be done (using summation in the case of OR connections and minimization in the 
        case of AND) for unperturbed and perturbed states of genes and ratio of perturbed to 
        unperturbed will give the reaction coefficeints returned.

        INPUT
            <genes>  : List of gene IDs or a single gene ID.
            <coef>   : List of coefficients of perturbation for the corresponding genes in 
                       <genes> or a single value that will apply to all in <genes>. If
                       gene expression is provided, coefficient multiplied by expression
                       will give the perturbed expression level to be used in the calculation
                       of reaction coefficients. If gene expression is not provided,
                       only a coef of 0 is meaningful, as a quantitative assessment of
                       perturbation will not be done.
            <exp>    : Gene expression dictionary of the form {gene ID: expression}.
            
        OUTPUT
            rxn2coef : Dictionary of reactions to perturbation coeffcients, where the coefficients
                       indicate the ratio of reaction expressin in perturbed to unperturbed state.        
        '''
        if type(genes)==str:
            genes=[genes];
        N=len(genes);
        if not type(coef)==list:
            coef=N*[coef];
        if not exp and not coef==N*[0.]:
            print('WARNING: When expression levels are not provided, there is no use for non-zero coefficients for gene perturbation.')
        Srxns,Sgenes=set(),set();
        for gene in genes:
            Lrxn=self.genes_id[gene].reactions;
            for rxn in Lrxn:
                Srxns.add(rxn.id);
                Sgenes=Sgenes.union([g.id for g in rxn.genes]);
        gene2coef={genes[i]:coef[i] for i in range(N)};  
        if not exp:
            isQuantitative=False;
            exp={gene:1. for gene in Sgenes};
        else:
            if not Sgenes.issubset(exp):
                print(Sgenes.difference(exp));
                raise Exception('All genes associated with perturbed reactions must be included in the expression dictionary. The ones above are missing.');
            isQuantitative=True;
        exp_mod={}; 
        for gene in Sgenes:
            if gene in gene2coef:
                exp_mod[gene]=(exp[gene],exp[gene]*gene2coef[gene]);
            else:
                exp_mod[gene]=(exp[gene],exp[gene]);
        rxn2coef={};
        for rxn in Srxns:
            L=self.__getGPRtree(rxn);
            t=self.__traverseGPRtree_perturbation(L[0],exp_mod);
            if t[0]:
                rxn2coef[rxn]=t[1]/t[0];
            else:
                rxn2coef[rxn]=t[0]; #make this 1.0
        if not isQuantitative:
            for rxn in rxn2coef:
                if rxn2coef[rxn]:
                    rxn2coef[rxn]=1.;
        return rxn2coef
        
    def getTable(self,objects,attributes,annotations=None):
        '''Produces a dataframe for objects (i.e., metabolites, reactions, or genes) indicated by 
        <objects>, with columns as indicated by the <attributes> and (if applicable) <annotations> 
        lists which refer to object attributes and annotations, respectively. 
        
        Object IDs will form the row names. <annotations> will be seached in 'annotation' 
        attribute of the objects.
        
        E.g., let M be the mnm object in hand. Then:
            
            T=M.getTable(M.model.metabolites,['name', 'formula','charge'],['bigg.metabolite']);
            
            OR
            
            T=M.getTable(M.mets.values(),['id','name', 'formula','charge'],['bigg.metabolite']);
            
            both produce a table of all metabolites with columns indicating metabolite 
            name, formula,charge, and BIGG metabolite name.        
        '''
        return self.__getTable(objects,attributes,annotations)
        
    def getRxnTable(self,rxns,attributes=['name','reaction','flux','gene_name_reaction_rule','subsystem'],annotations=None):
        '''Returns a dataframe for reactions given in <rxns>, which can be a list of IDS, or a dataframe 
        with indices as IDs. getTable() function is used to generate the table. 
        See getTable() function for the usage of <attributes> and <annotations>.'''
        if isinstance(rxns,pd.DataFrame):
            rxns=[r for r in rxns.index];
        df=self.__getTable([self.rxns[rxn] for rxn in rxns],attributes,annotations);
        return df
        
    def addRxn(self,ID,attr,met_coef,met_attr={},verbose=True):
        '''Adds a reaction to cobra model and updates the mnm. 
            <ID>        :   Reaction ID
            <attr>      :   Reaction attributes other than ID
            <met_coef>  :   Dictionary of metabolite coefficients
            <met_attr>  :   Nested dictionary of metabolite attributes ({met: {attr:...}});
        '''
        rxnobj=cobra.Reaction(ID);
        for k in attr:
            setattr(rxnobj,k,attr[k]);
        met_coef_pr={}; 
        for met in met_coef:
            if met in self.mets:
                met_coef_pr[self.mets[met]]=met_coef[met];
            elif met in met_attr:
                metobj=cobra.Metabolite(met);
                for k in met_attr[met]:
                    setattr(metobj,k,met_attr[met][k]);
                met_coef_pr[metobj]=met_coef[met];
            else:
                raise Exception('Metabolite '+met+' was found neither in model nor in provided attribute dictionary (<met_attr>).')                               
        rxnobj.add_metabolites(met_coef_pr);
        self.model.add_reactions([rxnobj]);        
        self.update(verbose=verbose);
          
    def search(self,regex,where,exact=False,caseInsensitive=True):
        '''Searches regular expression <regex> in <where>, which can be one of
        "rxns", "genes", or "mets"; which stand for reactions, genes, and metabolites,
        respectively. To search for an exact word, enter a whole word as <regex>
        and set <exact> as True. By default, the search is case insensitive. To
        make it case sensitive, set <caseInsensitive>=False.
        
        The following fields are searched for each type of variable:
            reactions: name, id, annotation, reaction,gene_name_reaction_rule
            metabolites: name, id, annotation, formula
            genes: name, id, annotation
            
        Returns a dataframe of findings
        '''
        if not where in ['rxns','mets','genes']:
            try:
                where={'reactions':'rxns','reaction':'rxns','metabolites':'mets','metabolite':'mets','gene':'genes'}[where];
            except:
                raise Exception('<where> must be one of rxns, genes, or mets.');
        if exact:
            regex='\\b'+regex+'\\b';
        atts=['name','id'];
        if where=='mets':
            atts.append('formula');
        elif where=='rxns':
            atts.append('reaction');
            atts.append('gene_name_reaction_rule');
        #setting the search function
        if caseInsensitive:
            fsearch=lambda x: re.search(regex,x,re.IGNORECASE);
        else:
            fsearch=lambda x: re.search(regex,x);  
        #Searching
        L=[];    
        for x in getattr(self,where).values():
            for y in atts:
                s=getattr(x,y);
                if s:
                    if fsearch(s):
                        L.append((x.id,x.name,y,s));
            for k in x.annotation:
                s=str(x.annotation[k]);
                if s:
                    if fsearch(s):
                        L.append((x.id,x.name,k,s));
        #Wrapping up
        if L:
            df=pd.DataFrame(L);
            df.columns=['id','name','where','found'];
            df=df.set_index('id');
            return df
        else:
            print('Not found!');
            
    def show(self,what,attributes=['id','name','formula','charge','compartment','reaction','subsystem','gene_name_reaction_rule',\
            'gene_reaction_rule','lower_bound','upper_bound','annotation'],addAttributes=None):
        '''Shows the properties of <what> that are covered by <attributes> and <addAttributes> 
        lists together. <what> can be a metabolite, gene, or reaction. If <what> cannot be found as an 
        ID in these categories, then it will be searched as a whole word using the <search> function 
        and the first object that is found will be shown, along with what else was found as a note.
        
        In the case of genes and metabolites, associated reactions will also be tabulated. 
        
        Returns the attributes table and reactions table (if applicable). For a more detailed table of 
        reactions, use the getRxnTable() function with the returned dataframe of reactions.
        '''
        if addAttributes:
            attributes=attributes+addAttributes;
        rxns=None;
        comps=None;
        ids=None;
        obj=None;
        if what in self.uniqmets:
            obj=self.mets[self.uniqmets[what][0][0]];
            comps=','.join([t[1] for t in self.uniqmets[what]]);
            ids=','.join([t[0] for t in self.uniqmets[what]]);
            rxns=[];
            for t in self.uniqmets[what]:
                rxns=rxns+[r.id for r in self.mets[t[0]].reactions];
            rxns=list(set(rxns));
            cat='Metabolite';
        elif what in self.mets:
            obj=self.mets[what];
            cat='Metabolite';
        elif what in self.genes:
            obj=self.genes[what];
            cat='Gene';
        elif what in self.genes_id:
            obj=self.genes_id[what];
            cat='Gene';
        elif what in self.rxns:
            obj=self.rxns[what];
            cat='Reaction';
        else:
            print('Item not found in metabolite, gene, or reaction IDs. search() function will be used and the first item will be shown.');
            found=None;
            cats=('mets','genes','rxns');
            k=0;
            while found is None and k<3:
                where=cats[k];
                found=self.search(what,where,exact=True);
                cat={'rxns':'Reaction','genes':'Gene','mets':'Metabolite'}[where];
                k=k+1;
            if not found is None:
                itemid=found.index[0];
                if where=='genes':
                    where='genes_id';
                obj=getattr(self,where)[itemid];
                if len(found)>1:
                    print('Multiple items found in '+where+'. Showing the first of these. Others are '+','.join([x for x in found.index[1:]])+'\nThere may be other items in other categories');
                else:
                    print('Only one item found');
            else:
                'Entry not found in this model!';
        if obj:
            if hasattr(obj,'reactions') and not rxns:
                rxns=[r.id for r in getattr(obj,'reactions')];
            values=[cat];
            atts=['Category'];
            for att in attributes:
                if att=='compartment' and comps:
                    atts.append(att);
                    values.append(comps);
                elif att=='id' and ids:
                    atts.append(att);
                    values.append(ids);                    
                elif att=='annotation' and hasattr(obj,att):
                    d=getattr(obj,att);
                    for k in d:
                        atts.append(k);
                        values.append(d[k]);                        
                elif hasattr(obj,att):
                    atts.append(att);
                    values.append(getattr(obj,att));
            Tatt=pd.DataFrame({' ':values},index=atts);
            if rxns:
                Trxn=self.getRxnTable(rxns,['reaction','flux','subsystem','gene_name_reaction_rule']);
            else:
                Trxn=None;
            return Tatt,Trxn     
            
    def rshow(self,rxn,PrintSeparately=True):
        '''Same as show() for a reaction, except that a reaction string with metabolite names 
        will also be shown along with ids. Set <PrintSeparately> True to display rxn with IDs and names
        one under the other. so that you do not need to worry .'''
        if not rxn in self.rxns:
            raise Exception('This function only takes reaction IDs. The ID that you provided was not found in reactions.');
        df,_=self.show(rxn);
        s=self.rxns[rxn].reaction;
        for met in self.rxns[rxn].metabolites:
            s=re.sub(r'(^|\s){}($|\s)'.format(re.escape(met.id)), r'\1{}\2'.format(met.name), s);
        i=df.index.get_loc('reaction');
        df_new=pd.concat([df.iloc[:i+1],pd.DataFrame({' ':[s]},index=['reaction_wnames']),df.iloc[i+1:]]);
        if PrintSeparately:
            print('\n'+df_new.loc['reaction',' ']);
            print(df_new.loc['reaction_wnames',' ']);
        return df_new
            
    def fba_cobramodel(self,obj,Dir='max',parsimonious=False,verbose=True):
        '''
        Performs FBA on the cobra model (self.model) with a given objective.

        Syntax:
        ------
        self.fba('RM04432'); #To maximize the flux of reaction RM04432.
        self.fba('RM04432',Dir=min); #To minimize the flux of reaction RM04432 (maximizes flux in reverse direction, if the reaction is reversible).
        self.fba({'RM04432':1, 'RC00703':1}); #Reaction ID to coefficient dictionary #To maximize the sum of fluxes of the two reactions. For a weighted sum, change coefficients.

        Input
        -----
        obj (string or dict): Reaction ID or a dictionary of reaction IDs to coefficients for objective function. 
        
        Dir ('max' or 'min'): Direction of optimization indicating whether the objective function will be maximized or minimized.
        
        parsimonious (Boolean): If True, performs pFBA (total flux is minimized).
        
        verbose (Boolean): If True, prints the value of the objective function after optimization. This value is total minimum flux in the case of parsimonious FBA.
        
        Output:
        ------
        
        Reaction fluxes and sln attribute are updated with the optimization. Returns the objective value.
        '''
        self.__simpleObj(obj,Dir);
        if not parsimonious:
            self.sln=self.model.optimize();
        else:
            self.sln=cobra.flux_analysis.pfba(self.model);
        if verbose:
            print('Value of the objective function is '+str(self.sln.objective_value));
        return self.sln.objective_value  
            
    def setOptimizer(self,irreversible=False):
        '''Sets optlang interface using opt class. The Boolean <irreversible> indicates
        if the optimizer is to divide reactions into forward and reverse directions.
        '''
        self.optimizer=opt([self.rxns[r] for r in self.rxns],irreversible=irreversible);  
        self.objective_function,self.objective_status,self.objective_value=None,None,None;
        
    def addConstraint(self,name,rxn2coef,lb,ub):
        '''Adds a new constaint to the optimization problem with name <name>. Use a
        dictionary <rxn2coef> to define rxn name-coefficient pairs. Indicate
        the lower and upper boundaries for the sum of products from this dictionary
        as <lb> and <ub> respectively.'''
        unallowed=['pFBA','imat','imatplusplus_1','imatplusplus_2','imatplusplus_3','imatplusplus_4','fpa_sum'];
        self.__checkOptimizer();
        if name in unallowed:
            raise Exception('The given constraint name is not valid. The following names are taken and not allowed for custom constraints:\n'+str(unallowed));
        if name[:4]=='imat':
            raise Exception('A contraint name cannot start with "imat". Derivatives of this phrase are used throughout. Try another name.');
        if name in self.optimizer.Dcons:
            raise Exception('A constraint with the given name already exists. Use a different name or try editing methods.');
        if self.optimizer.isIrreversible:
            D={};
            for rxn in rxn2coef:
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    D[rxn_dir]=rxn2coef[rxn]*self.optimizer.rxn2varcoefs[rxn][rxn_dir];
            self.optimizer.addConstraint(D,name,lb,ub);
        else:
            self.optimizer.addConstraint(rxn2coef,name,lb,ub);
            
    def removeConstraint(self,con):
        '''Removes a constraint (<con>) from current optimizer.'''
        if con in self.optimizer.Dcons:
            self.optimizer.model.remove(self.optimizer.Dcons[con]);
            _=self.optimizer.Dcons.pop(con);   
            
    def changeRxnBoundaries(self,rxn2bounds):
        '''Uses a dictionary (<rxn2bounds>) in the form of {rxnid:[lb,ub]} to edit
        the variable boundaries in the current optimizer. These changes are temporary 
        and will be gone when optimizer is reset (e.g., to change the model
        reversibility type). A dictionary of the original boudaries is returned which
        can be later used to restore the original boundaries.'''
        self.__checkOptimizer();
        Dori=self.__getCurrentBoundaries(rxn2bounds.keys());
        if self.optimizer.isIrreversible:
            for rxn in rxn2bounds:
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    if rxn_dir[-1]=='f':
                        self.optimizer.Dvars[rxn_dir].lb=max(0,rxn2bounds[rxn][0]);
                        self.optimizer.Dvars[rxn_dir].ub=max(0,rxn2bounds[rxn][1]);
                    else:
                        self.optimizer.Dvars[rxn_dir].lb=max(0,-rxn2bounds[rxn][1]);
                        self.optimizer.Dvars[rxn_dir].ub=max(0,-rxn2bounds[rxn][0]);                        
        else:
            for rxn in rxn2bounds:
                self.optimizer.Dvars[rxn].lb=rxn2bounds[rxn][0];
                self.optimizer.Dvars[rxn].ub=rxn2bounds[rxn][1]; 
        return Dori
         
    def constrainObjective(self,name,tol=None,relTol=None):
        '''Constrains current objective function within a tolerance of current 
        objective value. The tolerance can be provided either as an absolute or 
        relative value.
        
        INPUT
            <name>   : Name of constraint.
            <tol>    : Absolute value of the tolerable difference from the objective.
            <relTol> : Relative value of the tolerable difference with respect
                       to the current value (e.g., 0.1 for up to 10% difference to be tolerated).
        '''
        try:
            absobj=abs(self.optimizer.model.objective.value);
        except:
            raise Exception('The optimizer must be set and have an objective value for this function to work.');
        if tol and relTol:
            print('WARNING: Since <tol> is given, <relTol> will be ignored.');
        if not tol is None:
            tolf=tol;
        elif not relTol is None:
            tolf=relTol*absobj
        else:
            raise Exception('You must provide either <tol> or <relTol>.');  
        self.optimizer.constrainObjective(tolf,name);              

    def FBA(self,objfo,direction='max',update=True):
        '''Performs flux balance analysis with the current optimizer using objective
        function <objfo>. The objective is maximized or minimimized as indicated
        by <direction>. The optimized objective value is returned. 
        
        INPUT
            <objfo>       : Objective function as dictionary (rxn -> coefficient) or a
                            single reaction ID.
            <direction>   : 'max' or 'min'
            <update>      : Boolean that indicates whether a solution should be generated
                            and incorporated.
        '''
        self.__checkOptimizer();
        conv={'max':1.,'min':-1.}[direction]; #converter of sign for minimization
        if not type(objfo)==dict:
            objfo={objfo:1.};
        if self.optimizer.isIrreversible:
            objf={};
            for rxn in objfo:
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    objf[rxn_dir]=conv*self.optimizer.rxn2varcoefs[rxn][rxn_dir];
        else:
            objf={k:objfo[k]*conv for k in objfo};
        self.__optimize(objf,direction='max',conv=conv);
        if update:
            self.sln=self.__getSolution();
        return self.objective_value
        
    def FVA(self,rxn):
        '''Performs flux variability analysis of reaction <rxn> and returns a
        tuple of minimum flux followed by maximum flux.'''
        vmax=self.FBA(rxn,'max',update=False);
        vmin=self.FBA(rxn,'min',update=False);
        return (vmin,vmax)
        
    def pFBA(self,objf,direction='max',update=True,leaveConstraint=False):
        '''Performs parsimonious flux balance analysis with the current optimizer
        model. Similar to FBA() function except that total flux is minimized in the
        end.

        INPUT
            <objf>        : Objective function as dictionary (rxn -> coefficient) or a
                            single reaction ID.
            <direction>   : 'max' or 'min', to maximize or minimize the objective, respectively.
            <update>      : Boolean that indicates whether a solution should be generated
                            and incorporated.        
        <leaveConstraint> : Boolean that indicates whether the constraint used to constrain
                            the objective before flux minimization should be left on.  
        '''
        self.__checkOptimizer();
        self.removeConstraint('pFBA');
        xo=self.FBA(objf); #FBA solution
        self.optimizer.constrainObjective(self.tol,'pFBA');
        vmin_total=self.minimizeFluxSum(self.rxns.keys());  
        if not leaveConstraint:
            self.removeConstraint('pFBA'); 
        if update:
            self.sln=self.__getSolution(False);                                    
        return (xo,vmin_total)
    
    def IMAT(self,high,low,eps=0.01,rxn2eps={},relaxedHigh=False,getSummary=True,\
            minimizeFluxSum=True,constrain=True,tol=0,flux_thre=1E-4,eps_tol=0.01,bin_tol=0.001):
        '''Performs IMAT integration according to Shlomi et al., 2008 (PMID: 18711341). 
        Briefly, first, highly and lowly expressed genes are converted to highly and lowly
        expressed reactions, designated as "ON" and "OFF" reactions, respectively.
        Then binary variables are used to maximize the total number of ON reactions
        that carry flux and OFF reactions that do not. The goal of optimization is
        to maximize the sum of these binary variables (ysum). The maximum indicates a
        best fit between gene expression and flux state.

        INPUT
            <high>            : Highly expressed genes as gene IDs.
            <low>             : Lowly expressed genes as gene IDs.
            <eps>             : Generic flux threshold to impose on ON reactions.
            <rxn2eps>         : Dictionary of specific flux thresholds in the form
                                {rxn:(eps_reverse,eps_forward)}. Reactions not included 
                                in this dictionary will be given <eps> in both directions.
            <relaxedHigh>     : Boolean that indicates whether relaxed rules should be used 
                                when determining ON reactions. This means, any reaction 
                                associated with a highly expressed gene will be an ON
                                reaction, unless it is an OFF reaction because of lowly 
                                expressed genes.
            <getSummary>      : If True, then a summary of results will be produced. The 
                                summary contains ON and OFF reactions and their flux
                                states. The evaluations are based on both binary variables
                                that are associated with each reaction, and on the flux
                                values in comparison to flux thresholds and zero flux tolerance.
            <minimizeFluxSum> : If True, total flux will be minimized while the best fit
                                is maintaned with a tolerence as indicated by <tol>. 
            <constrain>       : Boolean that indicates whether the optimizer model should be 
                                left constrained by ysum based on a tolerance indicated 
                                by <tol> (see below). This is useful for subsequent operations like FVA.
                                If flux sum is minimized (<minimizeFluxSum>=True), then this
                                parameter is set True and the entry will be ignored.
            <tol>             : Tolerance for difference from best fit. Can be an integer,
                                thus indicating ysum-<tol> reactions will have to fit
                                expression data, or a ratio in which case ysum*(1-<tol>) will
                                be the lower limit for the number of fitting reactions.
            <flux_thre>       : Flux threshold that indicates a non-zero flux in a reaction 
                                after fitting is done. Only used in making the summary and 
                                has no effect on the outcome.
            <eps_tol>         : Tolerance for flux thresholds used when evaluating if flux-carrying
                                reactions carry the minimum imposed flux. Only used in making the summary and 
                                has no effect on the outcome.
            <bin_tol>         : Tolerance for binary 1 used when evaluating if a reaction is
                                in ON or OFF state after fitting. Only used in making the summary and 
                                has no effect on the outcome.
                                
        OUTPUT                        
            The following are returned in the respective order:
                ysum, total flux, summary
                None is returned instead of summary if it is not available.
        '''
        self.__checkOptimizer();
        if self.optimizer.isIrreversible:
            raise Exception('Optimizer is set as an irreversible model. IMAT runs only with regular models.');
        self.__removeIMATConstraints(); #in case there are prior IMAT or IMAT++ related constraints   
        print('Preparing solver');
        #Determine off and on reactions
        Roff=self.__getONorOFFrxns(low,'OFF');
        if relaxedHigh:
            Ron=set();
            for gene in high:
                for rxn in self.genes_id[gene].reactions:
                    Ron.add(rxn.id);
            Ron=Ron.difference(Roff);
        else:
            Ron=self.__getONorOFFrxns(high,'ON');  
        #Set flux thresholds for on reactions
        Ron2eps={};
        for rxn in Ron:
            if rxn in rxn2eps:
                Ron2eps[rxn]=rxn2eps[rxn];
            else:
                Ron2eps[rxn]=(eps,eps);
        #Add reaction variables        
        for rxn in Ron:
            self.optimizer.addOnRxn(rxn,Ron2eps[rxn]);
        for rxn in Roff:
            self.optimizer.addOffRxn(rxn);
        #Set the objective
        objf={};
        for rxn in self.optimizer.var2binvars:
            for var in self.optimizer.var2binvars[rxn]:
                objf[var]=1.;
        #optimize (and constrain if necessary)
        print('Optimizing Ron and Roff');
        self.__optimize(objf,direction='max');
        ysum=self.objective_value;
        if minimizeFluxSum or constrain:
            if type(tol)==int:
                tolf=tol+self.tol;
            elif tol<1.:
                tolf=ysum*(tol)+self.tol;
            else:
                print('<tol> must either be an integer or a ratio smaller than 1. Given value is ignored, zero tolerance (tol-0) assumed when constraining integer variables.')
                tolf=self.tol;
            if minimizeFluxSum:
                self.optimizer.constrainObjective(tolf,'imat');  
                print('Minimizing flux sum');
                vtotal=self.minimizeFluxSum(self.rxns.keys(),getDifferentials=False);
            else:
                self.sln=self.__getSolution(getDifferentials=False);
                vtotal=self.getTotalFlux();
        else:
            self.sln=self.__getSolution(getDifferentials=False);
            vtotal=self.getTotalFlux();
        #summarize results    
        if getSummary:
            summary={'on':{'success':{},'failure':{},'actual_failure':{},'actual_success':{}},\
            'off':{'success':{},'failure':{},'actual_failure':{},'actual_success':{}}};
            summary['on']['rxns']=[rxn for rxn in Ron.intersection(self.optimizer.var2binvars)];
            summary['off']['rxns']=[rxn for rxn in Roff.intersection(self.optimizer.var2binvars)];
            for rxn in summary['on']['rxns']: 
                binvals=[self.optimizer.Dvars[x].primal for x in self.optimizer.var2binvars[rxn]];
                var_val=list(zip(self.optimizer.var2binvars[rxn],binvals));
                v=self.optimizer.Dvars[rxn].primal;
                if max(binvals)>(1.-bin_tol):
                    summary['on']['success'][rxn]=[v,var_val];
                    if v>-(Ron2eps[rxn][0]*(1.-eps_tol)) and v<Ron2eps[rxn][1]*(1.-eps_tol):
                        summary['on']['actual_failure'][rxn]=v;
                else:
                    summary['on']['failure'][rxn]=[self.optimizer.Dvars[rxn].primal,var_val];
                    if abs(v)>=flux_thre:
                        summary['on']['actual_success'][rxn]=v;
            for rxn in summary['off']['rxns']:
                binval=self.optimizer.Dvars[self.optimizer.var2binvars[rxn][0]].primal;
                v=self.optimizer.Dvars[rxn].primal;
                if binval>(1.-bin_tol):
                    summary['off']['success'][rxn]=self.optimizer.Dvars[rxn].primal;
                    if abs(v)>=flux_thre:
                        summary['off']['actual_failure'][rxn]=v;
                else:
                    summary['off']['failure'][rxn]=self.optimizer.Dvars[rxn].primal;
                    if abs(v)<self.tol:
                        summary['off']['actual_success'][rxn]=v;
        else:
            summary=None;   
        if constrain and not minimizeFluxSum:
            print('Constraining best fit for subsequent applications');
            self.optimizer.constrainObjective(tolf,'imat'); 
        return (ysum,vtotal,summary)

    def IMATplusplus(self,high,low,rare,eps=0.01,rxn2eps={},relaxedHigh=False,getSummary=True,\
            minimizeFluxSum=True,rescueLatent=False,tol=0,lowflux_tol=1E-5,allflux_tol_ratio=0.05,\
            flux_thre=1E-4,eps_tol=0.99,bin_tol=0.999,addLowRxns=set(),addOffRxns=set(),neverLowOrOff=set(),weights=None):
        '''Performs IMAT++ integration according to Yilmaz,Li et al., 2020 (PMID: 33022146). 
        Briefly, first, highly and rarely expressed genes are convreted to highly and rarely
        expressed reactions, designated as "ON" and "OFF" reactions, respectively. These reactions
        are associated with binary variables. Then, each highly expressed gene, designated an "ON" gene,
        is associated with an integer variable equal to less than the sum of all "ON" reaction binaries 
        associated with that gene. Gene variables are also constrained to be less than 1, which means
        it is sufficient to have one of the multiple reactions associated with a gene (if applicable)
        carry flux to activate the highly expressed gene; this is different from IMAT since IMAT demands all
        associated reactions of a highly expressed gene to carry flux. To achieve optimal agreement between
        flux distribution and highly and rarely expressed genes ("best fit"), the sum of ON gene variables and OFF 
        reaction variables (ysum) are maximized. Then keeping this optimal state, the flux sum of lowly expressed
        reactions ("LOW" reactions, defined based on lowly and rarely expressed genes) is minimized. Then, 
        keeping this minimum flux state for lowly expressed reactions as well as the ON/OFF agreement, 
        total flux is minimized to obtain a preliminary flux distribution (PDF). Finally, latent reactions 
        (i.e., ON reactions without flux whose reactants [metabolites] are available [active in the current 
        flux distribution] and which are associated with only highly expressed genes], are rescued 
        (made carry flux) iteratively, by maximizing the sum of binary variables associated with such reactions. 
        Total flux is minimized in the end which yields the optimal flux distribution (OFD).
        
        To obtain a PFD, set <rescueLatent>=False, <minimizeFluxSum>=True (default option)
        To obtain an OFD, set <rescueLatent>=True

        INPUT
            <high>            : Highly expressed genes as gene IDs.
            <low>             : Lowly expressed genes as gene IDs.
            <rare>            : Rarely expressed genes as gene IDs.
            <eps>             : Generic flux threshold to impose on ON reactions.
            <rxn2eps>         : Dictionary of specific flux thresholds in the form
                                {rxn:(eps_reverse,eps_forward)}. Reactions not included 
                                in this dictionary will be given <eps> in both directions.
            <relaxedHigh>     : Boolean that indicates whether relaxed rules should be used 
                                when determining ON reactions. This means, any reaction 
                                associated with a highly expressed gene will be an ON
                                reaction, unless it is an LOW reaction because of lowly 
                                and rarely expressed genes.
            <getSummary>      : If True, then a summary of results will be produced. The 
                                summary contains ON genes, OFF reactions and their flux
                                states. These evaluations are based on both binary variables
                                that are associated with each reaction, and on the flux
                                values in comparison to flux thresholds and zero flux tolerance.
                                If applicable (see below), the status of detected latent reactions 
                                are shown as rescued and not rescued. Zero flux latent reactions 
                                that have been rescued are those which appear during iteration but later
                                lose their latent status as the flux distribution evolves.
            <minimizeFluxSum> : If True, total flux will be minimized while the best fit
                                (sum of ON genes and OFF reactions) and minimized sum of fluxes of 
                                LOW reactions are maintaned with tolerences as indicated by 
                                <tol> and <lowflux_tol>, respectively. If <rescueLatent> is
                                True (see below), this will be automatically True and the entry
                                will be ignored.
            <rescueLatent>    : If True, latent reactions will be iteratively detected and rescued.
                                At each step, total flux will be minimized to obtain a valid
                                flux distribution and constrained using <allflux_tol_ratio> (see below).
            <tol>             : Tolerance for difference from best fit. Can be an integer,
                                thus indicating ysum-<tol> genes/reactions will have to fit
                                expression data, or a ratio in which case ysum*(1-<tol>) will
                                be the lower limit for the number of fitting reactions.
            <lowflux_tol>     : Tolerance for sum of fluxes of LOW reactions. This sum will be
                                upper bound by x+<lowflux_tol>, where x is the current value.
           <allflux_tol_ratio>: Tolerance for sum of fluxes of all reactions. This sum will be
                                upper bound by x*(1+<allflux_tol_ratio>), where x is the current value.
                                This is a critical parameter for rescuing latent reactions.
                                If <minimizeFluxSum> is True and <rescueLatent> is False, set
                                this variable as None to avoid constraining total flux.
            <flux_thre>       : Flux threshold that indicates a non-zero flux in a reaction.
                                Used when detecting latent reactions and when preparing summary. 
            <eps_tol>         : Tolerance for flux thresholds used when evaluating if flux-carrying
                                reactions carry the minimum imposed flux. Only used in making the summary and 
                                has no effect on the outcome.
            <bin_tol>         : Tolerance for binary 1 used when evaluating if a reaction is
                                in ON or OFF state after fitting. Only used in making the summary and 
                                has no effect on the outcome.
            <addLowRxns>      : Set of additional rxns to be treated in "LOW" state even if they are not associated
                                with lowly or rarely expressed genes or any genes.
            <addOffRxns>      : Set of additional rxns to be treated in "OFF" state even if they are not associated
                                with rarely expressed genes or any genes.  
            <neverLowOrOff>   : Set of rxns that are NOT to be treated in "LOW" or "OFF" state even if they are associated
                                with lowly or rarely expressed genes.
            <weights>         : If provided as a dictionary of reactions to weights, the algorithm will minimize a
                                weighted sum of fluxes instead of treating all fluxes equally. All weights not included in
                                this disctionary will be taken as 1.0, as in the default case of flux minimization.
            .                                
        OUTPUT                        
            The following are returned in the respective order:
                ysum, total flux of LOW reactions,total flux of all reactions, summary
                None is returned instead of summary if it is not available.
                 
        '''
        self.__checkOptimizer();
        if self.optimizer.isIrreversible:
            raise Exception('Optimizer is set as an irreversible model. IMAT++ runs only with regular models.');
        Dweight={rxn:1. for rxn in self.rxns};
        if weights:
            Lthrow=[];
            for rxn in weights:
                if rxn in Dweight:
                    Dweight[rxn]=weights[rxn];
                else:
                    Lthrow.append(rxn);
            if Lthrow:
                print(Lthrow);
                print('Weights were provided for the reactions above but they are ingnored since the model does not have these reactions');
        self.__removeIMATConstraints(); #in case there are prior IMAT or IMAT++ related constraints 
        print('Preparing solver');
        #Determine off,low, on, and (if applicable) superhigh reactions
        print('Determining OFF reactions');        
        Roff=self.__getONorOFFrxns(rare,'OFF');
        Roff=Roff.union(addOffRxns);
        Roff=Roff.difference(neverLowOrOff);
        print('Determining Low reactions');        
        Rlow=self.__getONorOFFrxns(set(low).union(rare),'OFF');
        Rlow=Rlow.union(addLowRxns);
        Rlow=Rlow.difference(neverLowOrOff);
        print('Determining ON reactions');        
        if relaxedHigh:
            Ron=set();
            for gene in high:
                for rxn in self.genes_id[gene].reactions:
                    Ron.add(rxn.id);
            Ron=Ron.difference(Rlow);
        else:
            Ron=self.__getONorOFFrxns(high,'ON');  
        if rescueLatent: 
            Rsuper=set();
            for rxn in Ron:
                if np.all([gene.id in high for gene in self.rxns[rxn].genes]):
                    Rsuper.add(rxn);    
        #Set flux thresholds for on reactions
        Ron2eps={};
        for rxn in Ron:
            if rxn in rxn2eps:
                Ron2eps[rxn]=rxn2eps[rxn];
            else:
                Ron2eps[rxn]=(eps,eps);
        #Add reaction variables and ON genes        
        for rxn in Ron:
            self.optimizer.addOnRxn(rxn,Ron2eps[rxn]);
        for rxn in Roff:
            self.optimizer.addOffRxn(rxn);
        Gon=self.__getONgenes(high,Ron);
        for gene in Gon:
            self.optimizer.addOnGene(gene,Gon[gene]);            
        #Set the objective
        objf={};
        for rxn in Roff:
            objf[self.optimizer.var2binvars[rxn][0]]=1.;
        for gene in Gon:
            objf[self.optimizer.gene2intvars[gene]]=1.;
        print('Optimizing Gon and Roff');
        #optimize integers
        self.__optimize(objf,direction='max');
        ysum=self.objective_value;
        #constrain integers
        if type(tol)==int:
            tolf=tol+self.tol;
        elif tol<1.:
            tolf=ysum*(tol)+self.tol;
        else:
            print('<tol> must either be an integer or a ratio smaller than 1. Given value is ignored, zero tolerance (tol-0) assumed when constraining integer variables.')
            tolf=self.tol;
        self.optimizer.constrainObjective(tolf,'imatplusplus_1');
        print('Minimizing low flux');
        #minimize low
        if minimizeFluxSum or rescueLatent:
            vmin_low=self.minimizeFluxSum(Rlow,update=False,getDifferentials=False);
        else:
            vmin_low=self.minimizeFluxSum(Rlow,update=True,getDifferentials=False);
        #constrain
        if minimizeFluxSum or rescueLatent:
            #constrain fluxes
            self.optimizer.constrainObjective(lowflux_tol,'imatplusplus_2');
            print('Minimizing total flux');
            #minimize total flux
            vtotal=self.minimizeFluxSum(Dweight.copy(),getDifferentials=False);
        else:
            vtotal=self.getTotalFlux();
        #rescue latent reactions
        if rescueLatent:
            print('Rescuing latent reactions');
            latentrxns=self.__getLatentRxns(Rsuper,flux_thre);
            print(latentrxns);            
            alllatentrxns=latentrxns.copy();
            N,Npre,k=len(latentrxns),-1,0;
            while latentrxns and N>Npre:
                k=k+1;
                print('\tIteration #'+str(k)+' to rescue latent.');
                print('\tConstraining total flux and forcing latent reactions to carry flux.');
                self.removeConstraint('imatplusplus_3');
                self.constrainObjective('imatplusplus_3',relTol=allflux_tol_ratio);                
                Npre=N;
                objf={var:1 for var in alllatentrxns};              
                self.__optimize(objf,direction='max');
                yfit_latent=self.objective_value;
                print('\t'+str(int(yfit_latent))+' out of '+str(N)+' latent reactions are rescued according to solver output.');
                self.removeConstraint('imatplusplus_4');                
                self.constrainObjective('imatplusplus_4',0.);
                print('\tMinimizing total flux and recalculating latent reactions');
                vtotal=self.minimizeFluxSum(Dweight.copy(),getDifferentials=False);
                print('\tTotal flux is '+str(vtotal)); 
                latentrxns=self.__getLatentRxns(Rsuper,flux_thre);               
                alllatentrxns=alllatentrxns.union(latentrxns);
                N=len(alllatentrxns);  
                print(latentrxns);
                print('\t'+str(N-Npre)+' new latent reactions found.')
        #summarize results    
        if getSummary:
            print('Summarizing');
            summary={'on':{'success':{},'failure':{},'actual_failure':{},'actual_success':{}},\
            'off':{'success':{},'failure':{},'actual_failure':{},'actual_success':{}}};
            summary['on']['genes']=[gene for gene in Gon];
            summary['on']['rxns']=[rxn for rxn in Ron];
            summary['off']['rxns']=[rxn for rxn in Roff.intersection(self.optimizer.var2binvars)];
            summary['low']={rxn:self.optimizer.Dvars[rxn].primal for rxn in Rlow};
            if rescueLatent:
                summary['latent']={'rescued':{},'notrescued':set()};
                summary['latent']['notrescued']=set([x[:-3] for x in latentrxns]);
                for rxn in alllatentrxns.difference(latentrxns):
                    summary['latent']['rescued'][rxn[:-3]]=self.optimizer.Dvars[rxn[:-3]].primal;
            for gene in summary['on']['genes']: 
                intval=self.optimizer.Dvars[self.optimizer.gene2intvars[gene]].primal;
                dfailure={};
                dsuccess={};                
                if intval>=bin_tol:
                    for rxn in Gon[gene]:
                        v=self.optimizer.Dvars[rxn].primal;
                        if v>-(Ron2eps[rxn][0]*eps_tol) and v<Ron2eps[rxn][1]*eps_tol:
                            binvals=[self.optimizer.Dvars[x].primal for x in self.optimizer.var2binvars[rxn]];
                            dfailure[rxn]=[v,list(zip(self.optimizer.var2binvars[rxn],binvals))];
                        else:
                            dsuccess[rxn]=v;
                    if not dsuccess:
                        summary['on']['actual_failure'][gene]=dfailure.copy();
                    summary['on']['success'][gene]=dsuccess.copy();
                else:
                    for rxn in Gon[gene]:
                        v=self.optimizer.Dvars[rxn].primal;                    
                        if abs(v)>=flux_thre:
                            dsuccess[rxn]=v;
                        else:
                            dfailure[rxn]=v;
                    if dsuccess:
                        summary['on']['actual_success'][gene]=dsuccess.copy();
                    summary['on']['failure'][gene]=dfailure.copy();
            for rxn in summary['off']['rxns']:
                binval=self.optimizer.Dvars[self.optimizer.var2binvars[rxn][0]].primal;
                v=self.optimizer.Dvars[rxn].primal;
                if binval>bin_tol:
                    summary['off']['success'][rxn]=self.optimizer.Dvars[rxn].primal;
                    if abs(v)>=flux_thre:
                        summary['off']['actual_failure'][rxn]=v;
                else:
                    summary['off']['failure'][rxn]=self.optimizer.Dvars[rxn].primal;
                    if abs(v)<self.tol:
                        summary['off']['actual_success'][rxn]=v;
        else:
            summary=None;  
        if (allflux_tol_ratio and minimizeFluxSum) or rescueLatent:
            print('Constraining total flux for subsequent applications')
            self.constrainObjective('imatplusplus_3',relTol=allflux_tol_ratio);                    
        return (ysum,vmin_low,vtotal,summary)
        
    def setFPA(self,expTable,byproducts=None,networkTable=None,specialPenalties={},\
               inherit_array=None,inherit_rows=None,inherit_cols=None, noGenePenalty=1.,nanPenalty=1.,\
               plus=1.0, nandist=100.,mergeNetworks=False,setDistanceNet=True,verbose=True):
        '''Sets the variables and parameters for FPA (flux potential analysis).
        
        First, an expression penalty table of reactions is set. 
        Expression penalties are calculated as a function of GPR rules and gene 
        expression levels following previously established rules, unless otherwise 
        stated in the Input.
        
        Second,  a network of reaction nodes (instances of rxnnode class from
        MetabolicDistance) and (if provided) a pre-assembled (inherited) table 
        of reaction-reaction distances are established. The network or the table 
        is used in FPA to determine distances of all model reactions to the target 
        reaction. The availability of an inherited table of distances will reduce 
        the computation time of FPA.
        
        Third, a network of reactions that can carry flux in each condition is set
        as a table.
        
        Finally, the optimizer is set as irreversible if not already that way. Meaning,
        the model is converted to an irreversible model with originally reversible
        reactions replaced by two reactions each representing one of the directions.

        INPUT        
        --Related to Expression--        
            <expTable>        : Pandas dataframe of gene expression values with gene IDs as rows and 
                                conditions as columns. NaN is permitted as a value. Missing genes (model
                                genes not represented in the table) will be detected and NaN will be
                                inserted for all conditions for the expression level of these genes.
            <specialPenalties>: Dictionary of reactions that are assigned penalties not conforming to
                                the regular rules. The format of the dictionary is either as 
                                {rxnd:penalty (floating number)}, or as
                                {rxnd:{condition:penalty}}, or as a mixture of these, where, rxnd
                                indicates reaction ID with directionality (i.e., "_f" or "_r" in the end). 
                                In the first option the number representing the penalty will be distributed
                                to all conditions, while the second allows condition-specific penalties.
            <noGenePenalty>   : Penalty that will be assigned to reactions without any gene 
                                association (unless indicated in <specialPenalties>).
            <nanPenalty>      : Replaces NaN in penalty calculations. This may be due to missing
                                expression values.
            <plus>            : Positive number to be added to every value in the expression table. Make this zero
                                only if this addition was done apriori as zero expression values will
                                cause infinite penalties.
                                
        --Related to Distance--                                
            <byproducts>    : Metabolites that are to be ignored when setting inter-reaction
                              distances. Can be a list of unique metabolite names or 
                              compartmentalized metabolite IDs.
            <inherit_array> : A table (numpy array, matrix object, or data frame) of
                              previously determined distances. Rows must cover all 
                              irreversible reaction IDs according to the .irrevIDs attribute,
                              and cols include the targets for which distances have been calculated.
                              All reaction names must be provided with directionality (i.e., using "_f" and "_r"
                              suffixes to denote forward and reverse reactions, respectively).
                              If a numpy array is provided, then rows and columns need
                              to be provided separately as <inherit_rows> and <inherit_cols>.
            <inherit_rows>  : The rows of the <inherit_array>, if it is provided as a numpy array.
                              Must cover all reactions in the reaction network, meaning,
                              all reactions of the model with directionality according to
                              the .irrevIDs attribute. Use "_f" and "_r" suffixes to 
                              denote forward and reverse directions of a reaction, respectively. 
            <inherit_cols>  : Same as <inherit_rows> except that it identifies the columns of
                              the numpy array. These are the target reactions for which distances to
                              other reactions are available.
            <nandist>       : Distance value to fill in for missing distances. A large
                              number (>>10) is recommended, to represent infinite distance.
            <setDistanceNet>: Boolean that indicates if a distance network should be established
                              at all. Set this to False if all distances are provided in 
                              <inherit_array>.
        --Related to Network--                                
            <networkTable>  : Binary table (of ones and zeros) in the form of a
                              dataframe with rows indicating reactions with directionality
                              (i.e., with suffixes "_f" and "_r" for forward and reverse directions,
                              respectively) and columns the conditions. If a reaction
                              can carry flux in a condition, the corresponding value
                              must be 1, otherwise 0. Default value is none, in which
                              case, all reactions will be assumed to be able to carry 
                              flux in all conditions. The IDs of directional reactions 
                              must cover all of those in .irrevIDs attribute. Conditions 
                              must be the same as those provided by expression penalties table (see above).
            <mergeNetworks> : If true, columns of the <networkTable> will be merged into
                              one column using maximum function over each row and this 
                              combined network will be used for all conditions during
                              FPA. In other words, if a reaction can carry flux in
                              at least one condition, it will be part of the network
                              used in FPA.
                              
        --Miscelleaneous--                              
            <verbose>       : Boolean that indicates whether important messages will be displayed.

        OUTPUT
            Nothing is returned but the following objects are formed:
            --Related to Expression-- 
                .fpa_expTable    : A matrix object that represents the gene expression table modified with NaNs 
                                   for missing genes and addition of a number to every value (if applicable) (see above).
                .fpa_expPenalties: A matrix object of penalties with reactions as rows and conditions as columns. 
                .fpa_conditions  : List of conditions.  

        --Related to Distance-- 
                .fpa_distNet       : The network of reaction nodes that will be used to
                                     calculate distances from any target.
                .fpa_calculatedDist: A matrix object of available distance calculations.
                .fpa_nandist       : Scalar indicating what to replace missing distances with or None. 
                .fpa_byproducts    : List of by product metabolite IDs.  

        --Related to Network-- 
                .fpa_rxnNet   : matrix equivalent of the network table given or modified 
                                according to the input.                    

        --Other-- 
                .fpa_results  : Dictionary of FPA results that is created as an empty dictionary. See FPA() 
                                for details.
                .fpa_slns     : Dicionary of FPA solutions that is created as an empty dictionary. See FPA() 
                                for details.
                .fpa_activeTarget : Represents the ID of target reaction which current results belong to. Set
                                    initially as None.
        '''
        #set the optimizer right
        if hasattr(self,'optimizer'):
            if not self.optimizer.isIrreversible:
                self.setOptimizer(irreversible=True);
                if verbose:
                    print('WARNING: The optimizer is reset as irreversible. All new constraints added by addConstraints() function are lost.');
        else:
            self.setOptimizer(irreversible=True);
        self.__setExpressionPenalties(expTable,specialPenalties,noGenePenalty,nanPenalty,plus,verbose);
        self.__setDistanceCalculator(byproducts,inherit_array,inherit_rows,inherit_cols,nandist,setDistanceNet);
        self.__setConditionNetworks(networkTable,mergeNetworks);
        self.fpa_results,self.fpa_slns,self.fpa_activeTarget={'rFP':{},'FP':{}},{},None;

    def FPA(self,target,dist_order,allowance=1.0,conditions='all',specialDistances={},\
            correctDistanceLoops=True,blockEntryOf=None,blockExitOf=None,specialBoundaries=None,doSuperCond=True,\
            removeFPAconstraint=True,saveSolutions=True,reportingInstructions=None,verbose=True,**kwargs):
        '''Performs flux potential analysis for <target> reaction. FPA maximizes the flux of the <target> under specific
        constraints that are a function of relative gene expression levels. The maximized flux is called flux potential (FP).
        FP for each condition analyzed is eventually divided by that of the "super condition", a hypothetical condition 
        where relative gene expression is always the highest and all network reactions can carry flux, to obtain relative 
        FP (rFP) values. See PMID: 33022146 for the formulation of FPA.
        
        INPUT
            <target>      : The ID of the target reaction with directionality (i.e., with a suffix of "_f" or "_r").
                            The objective function of FPA is maximization of <target> flux and metabolic distances
                            are all calculated with respect to the <target>.
            <dist_order>  : Scalar used in decay function which reduces the expression penalty with metabolic distance.
                            The formula is 1/[(distance+1)^dist_order]. When multiplied by expression penalty, this formula 
                            determines final coefficient of each network reaction in the weighted sum constraint of FPA.
            <allowance>   : Flux allowance. Determines the upper limit for the weighted sum constraint in FPA.
            <conditions>  : Conditions for which FPA will be calculated. Use a list of a subset of conditions or
                            simply indicate 'all' for doing all conditions (default). "super_cond" should not be a part of 
                            this list.
        <specialDistances>: Dictionary of distances (from the target) that will overwrite calculated distances only 
                            during the current run.The format is {rxn_dir: distance}, where rxn_dir is reaction ID with 
                            directionality (i.e., with a suffix of "_f" or "_r").  
        <correctDistanceLoops>: Boolean indicating if deficient pathways encountered during metabolic distance calculations
                            should be checked and corrected. These are cases where the shortest path found to or from
                            the <target> from or to a reaction includes reversible reactions that were traveresed twice, 
                            once in each direction. Correction of these cases require an iterative process which may take
                            a long time. The benefit (or, the compromise of not correcting) is typically very small, so this
                            parameter can be made False to save on time. If a distance matrix is provided during the setting 
                            of FPA (see setFPA()), then this parameter is null, unless the matrix does not cover the <target>.
        <blockEntryOf>    : A list of metabolites or a single metabolite whose entry into the network through exchange, demand 
                            or sink reactions will be temprarily blocked by appropriate constraints. Both unique and specific
                            metabolite IDs are allowed. Useful for metabolite-level FPA.
        <blockExitOf>     : A list of metabolites or a single metabolite whose exit from the network through exchange, demand 
                            or sink reactions will be temprarily blocked by appropriate constraints. Both unique and specific
                            metabolite IDs are allowed. Useful for metabolite-level FPA.
        <specialBoundaries>: A dictionary of reaction boundaries to be applied temporarily during FPA. The format is:
                             {rxn ID: [LB,UB]}, where rxn ID must be from original model and not from irreversible model of the
                             optimizer.
        <saveSolutions>   : Boolean indicating whether flux solutions for each condition will be saved in .fpa_slns dictionary 
                            (as cobra.solution object).
        <reportingInstructions>: If provided, a report of FPA analyis will be generated and written to indicated files. 
                            The instructions must be provided as a tuple in the form of (file_path,pfx,sfx,dlm), where,
                            file_path is the path to the folder where files will be saved, pfx and sfx are prefix and suffix
                            that will be used in filename (the full name is /pfx/condition/_/<target>/sfx/), and dlm
                            is the delimiter used in the output files. The report will include IDs, stoichiometries,
                            fluxes, distances, expression penalties, overall coefficients (used in the weighted sum), and 
                            contributions to allowance of every reaction that carries flux in the pertaining condition. 
        <doSuperCond>     : Boolean indicating if FPA should be done on the super condition. This will be forcefully done regardless
                            of the input when the target is different than the active target (.fpa_activeTarget), which is the
                            reaction previously analysed without any interruption by other reactions. But if the active target
                            remains the same, then this parameter can be made False to save time. 
        <removeFPAconstraint>: Boolean indicating if the weighted sum constraint that limits total flux by <allowance> should be
                            removed at the end. Making this parameter False (default is True) is useful only when analyzing one
                            condition and wanting to use the allowance-restrained network for some subsequent optimization.
            <verbose>     : Boolean that indicates whether important messages will be displayed.  
            kwargs        : Other arguments related to eFPA (see eFPA() function), which include <dist_boundary> and <base>

        OUTPUT
              FP and rFP results in .fpa_results dictionary will be returned as a copy. In addition to .fpa_results, the 
            following will be modified with the results: 
                .fpa_slns : Dictionary of solutions (if applicable) (see above in <saveSolutions>).
                .fpa_activeTarget : Replaced by <target>. 
              Further, if <reportingInstructions> is provided, a report of FPA in each condition will be written as 
            instructed (see above in <reportingInstructions>).

        '''
        #check input
        if not hasattr(self,'fpa_expPenalties'):
            raise Exception('First, you need to set the variables needed for FPA using setFPA() function.')
        if not target in self.irrevIDs:
            raise Exception('The ID of the target reaction must be directional (i.e., with a suffix of "_f" or "_r") and be part of .irrevIDs attribute.');
        if set(specialDistances).difference(self.irrevIDs):
            print(set(specialDifferences).difference(self.irrevIDs));
            raise Exception('The above reaction IDs in <specialDistances> input were not recognized. Reaction IDs must be included in .irrevIDs attribute.');
        blockedMets={};
        if blockEntryOf:
            blockedMets.update({'entry':self.__getMetaboliteList(blockEntryOf,'blockEntryOf input error')});
        if blockExitOf:
            blockedMets.update({'exit':self.__getMetaboliteList(blockExitOf,'blockExitOf input error')});
        if specialBoundaries:
            S=set(specialBoundaries).difference(self.rxns);
            if S:
                print(S);
                raise Exception('The above reactions from <specialBoundaries> entry could not be recognized. Make sure to use original reaction IDs of the model without directionality.')
        #FPA vs eFPA
        is_efpa=False;
        if dist_order is None:
            dist_boundary=kwargs.get('dist_boundary', None);
            base=kwargs.get('base', None);
            hardBoundary=kwargs.get('hardBoundary', None);
            if dist_boundary is None or base is None or hardBoundary is None:
                raise Exception('''At least one eFPA parameter is missing although <dist_order> is not provided,
                and hence the function is automatically switched to eFPA. Here are the eFPA parameters you provided:
                dist_boundary: {}, hardBoundary: {}, base: {}'''.format(dist_boundary, hardBoundary, base))
            is_efpa=True;
        elif 'dist_boundary' in kwargs:
            print('WARNING: both FPA and eFPA parameters detected. eFPA parameters will be ignored since <dist_order>, an FPA parameter, is provided and FPA will be performed in the original manner.');            
        #set conditions to run
        if conditions=='all':
            conds=[x for x in self.fpa_conditions];
        elif type(conditions)==str:
            conds=[conditions];
        else:
            conds=[x for x in conditions];
        Sdiff=set(conds).difference(self.fpa_conditions);
        if Sdiff:
            print(Sdiff)
            raise Exception('The condition names above were not recognized. Conditions must be consistent with .fpa_conditions attribute.');
        if doSuperCond or not self.fpa_activeTarget==target:
            conds.append('super_cond');
        if not doSuperCond and not self.fpa_activeTarget==target and verbose:
            print('Super condition added to list of conditions to be analyzed despite <doSuperCond> being False, because the active target changed.');
        #set distance from target
        if not target in self.fpa_calculatedDist.cols:
            if self.fpa_distNet is None:
                raise Exception('The target is not in provided distance array as a column and setting of network was also waived, so there is no way of getting distances for this target. Reset FPA with either a distance array including this target or setting <setDistanceNet> as True');
            if target in self.fpa_distNet.Drxn2node:
                if correctDistanceLoops:
                    if verbose: print('Finding forward distances from target.');
                    Dpath,_,_=findPaths(target,self.fpa_distNet,quiet=not verbose);
                    Df=convertPaths2Distances(Dpath);
                    if verbose: print('Finding backward distances from target.');
                    Dpath,_,_=findPaths_backward(target,self.fpa_distNet,quiet=not verbose);
                    Dr=convertPaths2Distances(Dpath);
                else:
                    if verbose: print('Finding forward distances from target.');
                    self.fpa_distNet.pathForward(target);
                    Df={rxn:self.fpa_distNet.Dfirst_forward[rxn][1] for rxn in self.fpa_distNet.Dfirst_forward};
                    if verbose: print('Finding backward distances from target.');
                    self.fpa_distNet.pathBackward(target);
                    Dr={rxn:self.fpa_distNet.Dfirst_backward[rxn][1] for rxn in self.fpa_distNet.Dfirst_backward};  
                self.__appendCol2DistanceMatrix(target);
                for rxnd in self.irrevIDs:
                    if rxnd in Df:
                        if rxnd in Dr:
                            self.fpa_calculatedDist[rxnd,target]=min(Df[rxnd],Dr[rxnd]);
                        else:
                            self.fpa_calculatedDist[rxnd,target]=Df[rxnd];
                    elif rxnd in Dr:
                        self.fpa_calculatedDist[rxnd,target]=Dr[rxnd];
            else:
                self.__appendCol2DistanceMatrix(target);
                if verbose: print('WARNING: the target reaction '+target+' is unreachable from anywhere in the network and therefore all distances are preset as nandist before incorporating exceptions (that is, if exceptions are provided)');
        Ddist_this={rxnd:self.fpa_calculatedDist[rxnd,target] for rxnd in self.irrevIDs}; #extracted to incorporate special distances
        for rxnd in specialDistances:
            Ddist_this[rxnd]=specialDistances[rxnd];
        #FPA
        #block indicated metabolites (if applicable) and enforce special boundaries (if applicable)
        unblocker0={};
        if 'entry' in blockedMets:
            for met in blockedMets['entry']:
                unblocker0.update(self.__blockMetabolite(met,True,False,exclude=target[:-2]));
        if 'exit' in blockedMets:
            for met in blockedMets['exit']:
                unblocker0.update(self.__blockMetabolite(met,False,True,exclude=target[:-2]));
        if specialBoundaries:
            unblocker0.update(self.changeRxnBoundaries(specialBoundaries));
        #block other direction (if applicable), keep original constraint
        otherDir=target[:-1]+self.toggledir[target[-1]];
        unblocker1=None;        
        if otherDir in self.optimizer.Dvars:
            unblocker1=(self.optimizer.Dvars[otherDir].lb,self.optimizer.Dvars[otherDir].ub);
            self.optimizer.Dvars[otherDir].lb,self.optimizer.Dvars[otherDir].ub=0.,0.;
        for cond in conds:
            if self.fpa_rxnNet[target,cond]:
                if verbose:print('Now doing '+cond);   
                self.removeConstraint('fpa_sum'); #remove older constraint, if any
                #calculate reaction coefficients and add the weighted sum constraint               
                Dcoef={}; #dictionary of overall penalty for every reaction
                if not is_efpa:
                    for rxnd in self.irrevIDs:
                        Dcoef[rxnd]=self.fpa_expPenalties[rxnd,cond]/((Ddist_this[rxnd]+1)**dist_order);  
                elif not hardBoundary:
                    for rxnd in self.irrevIDs:
                        Dcoef[rxnd]=self.fpa_expPenalties[rxnd,cond]/(1+(base**(Ddist_this[rxnd]-dist_boundary)));                     
                else:
                    for rxnd in self.irrevIDs:
                        d=Ddist_this[rxnd];
                        if d>dist_boundary:
                            Dcoef[rxnd]=0.;
                        else:
                            Dcoef[rxnd]=self.fpa_expPenalties[rxnd,cond];                    
                self.optimizer.addConstraint(Dcoef,'fpa_sum',0.,allowance);
                #block inactive reactions, keep original boundaries
                unblocker2={};
                for rxnd in self.fpa_rxnNet.rows:
                    if not self.fpa_rxnNet[rxnd,cond]:
                        unblocker2[rxnd]=(self.optimizer.Dvars[rxnd].lb,self.optimizer.Dvars[rxnd].ub);
                        self.optimizer.Dvars[rxnd].lb,self.optimizer.Dvars[rxnd].ub=0.,0.;
                #optimize
                self.__optimize({target:1.},'max',1.);
                #save results and (if applicable) solutions
                self.fpa_results['FP'][cond]=self.objective_value;
                if saveSolutions:
                    self.fpa_slns[cond]=self.__getSolution();
                #write report (if applicable)
                if reportingInstructions:
                    file_path,pfx,sfx,dlm=reportingInstructions;
                    if not file_path[-1]=='/':file_path=file_path+'/';
                    fo=open(file_path+pfx+cond+'_'+target+sfx,'w');
                    fo.write(dlm.join(['ID','Stoichiometry','Flux','Distance','ExpPenalty','Coefficient','Contribution'])+'\n');
                    for rxnd in self.irrevIDs:
                        v=self.optimizer.Dvars[rxnd].primal;
                        if v>self.lowTol:
                            L=[rxnd,self.rxns[rxnd[:-2]].reaction,v,Ddist_this[rxnd],\
                               self.fpa_expPenalties[rxnd,cond],Dcoef[rxnd],Dcoef[rxnd]*v];
                            fo.write(dlm.join([str(x) for x in L])+'\n');
                    fo.close();
                #restore original boundaries
                for rxnd in unblocker2:
                    self.optimizer.Dvars[rxnd].lb,self.optimizer.Dvars[rxnd].ub=unblocker2[rxnd];                
            else:
                if verbose:print('The target is inactive in '+cond+' and only a zero flux potential will be saved. No solution exists.');
                self.fpa_results['FP'][cond],self.fpa_slns[cond]=0.,None;
        #calculate and save rFP
        for cond in conds:
            try:
                self.fpa_results['rFP'][cond]=self.fpa_results['FP'][cond]/self.fpa_results['FP']['super_cond'];
            except:
                self.fpa_results['rFP'][cond]=np.nan;
        #UNBLOCKS: IMPORTANT unblockers should follow the reverse order of blocking, to avoid 
        #changing original boundaries when they share the same reaction(s).
        #unblock blocked direction 
        if unblocker1:
            self.optimizer.Dvars[otherDir].lb,self.optimizer.Dvars[otherDir].ub=unblocker1; 
        #unblock metabolites and restore changed special boundaries
        if unblocker0:
            print(unblocker0);
            _=self.changeRxnBoundaries(unblocker0);  
            print(_);      
        #remove fpa constraint unless stated otherwise
        if removeFPAconstraint:
            self.removeConstraint('fpa_sum');
        #register target as active target
        self.fpa_activeTarget=target;
        return deepcopy(self.fpa_results)


    def eFPA(self, target, dist_boundary, hardBoundary=False, base=2, allowance=1.0, conditions='all', 
            specialDistances={},correctDistanceLoops=True, blockEntryOf=None, blockExitOf=None, 
            specialBoundaries=None, doSuperCond=True, removeFPAconstraint=True, saveSolutions=True, 
            reportingInstructions=None, verbose=True):
        '''Performs enhanced flux potential analysis for <target> reaction. Same as FPA() except for the decay function
        which reduces expression penalty of reactions. eFPA decay function uses different parameters (<dist_boundary> and 
        <hardBoundary or <base>) instead of only one parameter (dist_order) of the original FPA() function. 
        The formula of the decay function can be one of the following:
            
                1/[1+(base)^(distance-dist_boundary)]         (1)
                
                distance>dist_boundary => 0
                distance<=dist_boundary => 1                  (2)

            Equation 1 is used when <hardBoundary> is False (default) and Equation 2 when True.
            
        INPUT
            <target>      : The ID of the target reaction with directionality (i.e., with a suffix of "_f" or "_r").
                            The objective function of eFPA is maximization of <target> flux and metabolic distances
                            are all calculated with respect to the <target>.
            <dist_boundary> : Scalar used in decay function as shown above. Indicates the distance where expression penalty
                            goes through a transition to remarkably lower values.
            <hardBoundary>: Boolean that indicates whether Equation 1 (False) or 2 (True) is going to be used in 
                            calculating the decay (see above for equations).                            
            <base>        : Parameter of Equation 1 (see above), which determines how steep the transition from high to
                            low penalties is during the distance decay.
                            
        All other arguments and OUTPUT are as in FPA() function. See help in FPA() for further information.
        '''
        results=self.FPA(target, dist_order=None, dist_boundary=dist_boundary, hardBoundary=hardBoundary,base=base, allowance=allowance, 
            conditions=conditions, specialDistances=specialDistances, correctDistanceLoops=correctDistanceLoops,  
            blockEntryOf=blockEntryOf, blockExitOf=blockExitOf, specialBoundaries=specialBoundaries, doSuperCond=doSuperCond, 
            removeFPAconstraint=removeFPAconstraint, saveSolutions=saveSolutions, reportingInstructions=reportingInstructions,
            verbose=verbose);
        return results

        
    def setCompass_minus(self,expTable, noGenePenalty=1.,nanPenalty=1.,plus=1.0, fpaPenalties=False,verbose=True):
        '''Sets the variables and parameters for Compass (PMID: 34216539). These include setting an expression penalty table of reactions 
        and making the optimizer irreversible if it is not already set that way. The latter means the model is converted to an 
        irreversible model with originally reversible reactions replaced by two reactions each representing one of the directions. 
        
        The penalty rule for expression is based on PMID: 34216539, unless it is otherwise stated with the <fpaPenalties> parameter.
        According to this, the penalty is inverse of expression, after the addition of <plus> to expression values. This corresponds
        to the penalty = 1/(1+x) formula in the original publication, with x interpreted as total expression value for all genes
        associated with the reaction for which the penalty is being calculated. 
        
        NOTE: Smoothing (merging data from neighboring cells to mitigate technical noise) is not included in this version of Compass algorithm,
        hence the name "Compass_minus".

        INPUT        
            <expTable>        : Pandas dataframe of gene expression values with gene IDs as rows and 
                                conditions as columns. NaN is permitted as a value. Missing genes (model
                                genes not represented in the table) will be detected and NaN will be
                                inserted for all conditions for the expression level of these genes.
            <noGenePenalty>   : Penalty that will be assigned to reactions without any gene 
                                association.
            <nanPenalty>      : Replaces NaN in penalty calculations. This may be due to missing
                                expression values.
            <plus>            : Positive number to be added to every value in the expression table. Make this zero
                                only if this addition was done apriori as zero expression values will
                                cause infinite penalties.  
            <fpaPenalties>    : Boolean that indicates whether penalty calculations should be based on FPA algorithm (see above).
                                If False, then default Compass penalties will apply.                    
            <verbose>         : Boolean that indicates whether important messages will be displayed.

        OUTPUT
            Nothing is returned but the following objects are formed:

         .compass_expTable    : A matrix object that represents the gene expression table modified with NaNs 
                                for missing genes and addition of a number to every value (if applicable) (see above).
         .compass_expPenalties: A matrix object of penalties with reactions as rows and conditions (e.g., cells) as columns. 
         .compass_conditions  : List of conditions (e.g., cells).  
                   
            .compass_results  : Dictionary of Compass results that is created as an empty dictionary. See compass_minus() 
                                for details.
            .compass_slns     : Dicionary of Compass solutions that is created as an empty dictionary. See compass_minus() 
                                for details.
       .fcompass_activeTarget : Represents the ID of target reaction which current results belong to. Set
                                initially as None.
        '''
        #set the optimizer right
        if hasattr(self,'optimizer'):
            if not self.optimizer.isIrreversible:
                self.setOptimizer(irreversible=True);
                if verbose:
                    print('WARNING: The optimizer is reset as irreversible. All new constraints added by addConstraints() function are lost.');
        else:
            self.setOptimizer(irreversible=True);
        self.__setExpressionPenalties_compass(expTable,noGenePenalty,nanPenalty,plus,fpaPenalties,verbose);
        self.compass_results,self.compass_slns,self.compass_activeTarget={},{},None;

        
    def compass_minus(self,target,w=0.95,conditions='all',saveSolutions=True,reportingInstructions=None,verbose=True):
        '''Performs Compass (PMID:34216539) to calculate the resistance for the <target> reaction. This resistance is a minimized sum of
        products of reaction penalties and fluxes. The penalties should be set a priori by setCompass_minus() function. 
        
        The difference of the compass_minus() algorithm from the original (PMID:34216539) is the absence of the merging of single cells (here: conditions)
        with their nearest neighbors using a smoothing parameter. Instead, compass_minus is limited to evaluating individual conditions without
        any merging.
        
        INPUT
            <target>      : The ID of the target reaction with directionality (i.e., with a suffix of "_f" or "_r").
            <w>           : "Optimality slack parameter". The maximum value the flux (vmax) of the target reaction will be calculated first and then
                            the reaction flux will be constrained to this maximal value multiplied by <w>. 
            <conditions>  : Conditions to which compass_minus will be applied to. Use a list of a subset of conditions or
                            simply indicate 'all' for doing all conditions (default). 
        <saveSolutions>   : Boolean indicating whether flux solutions for each condition will be saved in .compass_slns dictionary 
                            (as cobra.solution object).
   <reportingInstructions>: If provided, a report of compass analyis will be generated and written to indicated files. 
                            The instructions must be provided as a tuple in the form of (file_path,pfx,sfx,dlm), where,
                            file_path is the path to the folder where files will be saved, pfx and sfx are prefix and suffix
                            that will be used in filename (the full name is /pfx/condition/_/<target>/sfx/), and dlm
                            is the delimiter used in the output files. The report will include IDs, stoichiometries,
                            fluxes, expression penalties, and contributions to sum (<w>*flux) of every reaction that carries flux 
                            in the pertaining condition. 
            <verbose>     : Boolean that indicates whether important messages will be displayed.  

        OUTPUT
            Results in .compass_results dictionary, which has the form {condition:resistance} will be returned as a copy. In addition, the 
            following will be modified with the results: 
                .compass_slns : Dictionary of solutions (if applicable) (see above in <saveSolutions>).
                .compass_activeTarget : Replaced by <target>. 
            Further, if <reportingInstructions> is provided, a report of Compass in each condition will be written as 
            instructed (see above in <reportingInstructions>).
            
        NOTE: The constraint on the target reaction (<w>*vmax) and the constraint on its reverse (if applicable; the algorithm constraints the
        exact reverse of the target reaction, if it exists, to zero flux to avoid a flux loop that involves the forward and reverse forms
        of the target and bypasses the rest of the network) will be removed when execution of compass_minus is finished.

        '''
        #check input
        if not hasattr(self,'compass_expPenalties'):
            raise Exception('First, you need to set the variables needed for compass_minus using setCompass_minus() function.')
        if not target in self.irrevIDs:
            raise Exception('The ID of the target reaction must be directional (i.e., with a suffix of "_f" or "_r") and be part of .irrevIDs attribute.');
        #set conditions to run
        if conditions=='all':
            conds=[x for x in self.compass_conditions];
        elif type(conditions)==str:
            conds=[conditions];
        else:
            conds=[x for x in conditions];
        Sdiff=set(conds).difference(self.compass_conditions);
        if Sdiff:
            print(Sdiff)
            raise Exception('The condition names above were not recognized. Conditions must be consistent with .compass_conditions attribute.');
        #Compass
        #block other direction (if applicable), keep original constraint
        otherDir=target[:-1]+self.toggledir[target[-1]];
        unblocker1=None;        
        if otherDir in self.optimizer.Dvars:
            unblocker1=(self.optimizer.Dvars[otherDir].lb,self.optimizer.Dvars[otherDir].ub);
            self.optimizer.Dvars[otherDir].lb,self.optimizer.Dvars[otherDir].ub=0.,0.;
        vmax=abs(self.FBA(target[:-2],{'_f':'max','_r':'min'}[target[-2:]],False))
        if vmax:
            unblocker2=(self.optimizer.Dvars[target].lb,self.optimizer.Dvars[target].ub)
            self.optimizer.Dvars[target].lb=w*vmax;
            for cond in conds:
                if verbose:print('Now doing '+cond);   
                #calculate reaction coefficients and add the weighted sum constraint               
                Dcoef={}; #dictionary of overall penalty for every reaction
                for rxnd in self.irrevIDs:
                    Dcoef[rxnd]=self.compass_expPenalties[rxnd,cond]  
                #optimize
                self.__optimize(Dcoef,'min',1.);
                #save results and (if applicable) solutions
                self.compass_results[cond]=self.objective_value;
                if saveSolutions:
                    self.compass_slns[cond]=self.__getSolution();
                #write report (if applicable)
                if reportingInstructions:
                    file_path,pfx,sfx,dlm=reportingInstructions;
                    if not file_path[-1]=='/':file_path=file_path+'/';
                    fo=open(file_path+pfx+cond+'_'+target+sfx,'w');
                    fo.write(dlm.join(['ID','Stoichiometry','Flux','ExpPenalty','Contribution'])+'\n');
                    for rxnd in self.irrevIDs:
                        v=self.optimizer.Dvars[rxnd].primal;
                        if v>self.lowTol:
                            L=[rxnd,self.rxns[rxnd[:-2]].reaction,v,Dcoef[rxnd],Dcoef[rxnd]*v];
                            fo.write(dlm.join([str(x) for x in L])+'\n');
                    fo.close();
            #restore original boundaries
            self.optimizer.Dvars[target].lb,self.optimizer.Dvars[target].ub=unblocker2;                
        else:
            if verbose:print('The target is inactive in the presumed direction and NaN resistance values will be assigned to every condition. No solutions exist.');
            for cond in conds:
                self.compass_results[cond]=np.nan
                self.compass_slns[cond]=None
        #unblock blocked direction 
        if unblocker1:
            self.optimizer.Dvars[otherDir].lb,self.optimizer.Dvars[otherDir].ub=unblocker1; 
        #register target as active target
        self.compass_activeTarget=target;
        return deepcopy(self.compass_results)


    def __setExpressionPenalties_compass(self,expTable,noGenePenalty,nanPenalty,plus,fpaPenalties,verbose):
        '''
        Subroutine of setCompass(). See setCompass() instructions related to expression.               
        '''
        if fpaPenalties:
            self.__setExpressionPenalties(expTable,{},noGenePenalty,nanPenalty,plus,verbose);
            for att in ('expTable','expPenalties','conditions'):
                compass_name,fpa_name='compass_'+att,'fpa_'+att;                
                setattr(self,compass_name,getattr(self,fpa_name))
                delattr(self,fpa_name)
        else:
            #quality check of the table and modifications  
            if isinstance(expTable,pd.DataFrame):
                self.compass_expTable=matrix(expTable.to_numpy(),list(expTable.index),list(expTable.columns));
            elif isinstance(expTable,matrix):
                self.compass_expTable=expTable.copy();
            else:            
                raise Exception('Expressions must be provided either as a dataframe or as a matrix (not numpy) object');
            if not np.issubdtype(self.compass_expTable.vals.dtype,np.number):
                raise Exception('Data is not numeric. The values of the expression table should constitute a numerical numpy array.');
            #conditions
            self.compass_conditions=[col for col in self.compass_expTable.cols];
            Ncond=len(self.compass_conditions);
            #genes
            missingGenes=set(self.genes_id).difference(self.compass_expTable.rows);
            Nmiss=len(missingGenes);
            if Nmiss==len(self.genes):
                raise Exception('None of the genes in model were found in the provided expression table. Please make sure the row indices are gene IDs.');
            #value addition
            if plus<=0:
                if verbose:
                    print('WARNING: the variable <plus> must be a positive number. The given number is ignored and the default value of 1 is used for <plus>.')
                plus=1.;
            #convert missing gene values to nans
            if Nmiss:
                nans=np.empty((Nmiss,self.compass_expTable.colno));
                nans[:]=np.nan;            
                self.compass_expTable=self.compass_expTable.appendRows(nans,list(missingGenes));
            #calculate regular penalties
            rxn2pen={};
            gpr2pen={};
            for rxn in self.rxns:
                if self.rxns[rxn].genes:
                    gprrule=self.rxns[rxn].gene_reaction_rule;
                    if gprrule in gpr2pen:
                        rxn2pen[rxn]=gpr2pen[gprrule];
                    else:
                        rxn2pen[rxn]=[]
                        for cond in self.compass_conditions:
                            c=plus;
                            for gene in self.rxns[rxn].genes:
                                c=c+self.compass_expTable[gene.id,cond]
                            rxn2pen[rxn].append(1./np.nan_to_num(c,nan=nanPenalty));
                    gpr2pen[gprrule]=[x for x in rxn2pen[rxn]];
                else:
                    rxn2pen[rxn]=Ncond*[noGenePenalty];
            dummyvals=np.empty((len(self.irrevIDs),len(self.compass_conditions)));
            self.compass_expPenalties=matrix(dummyvals,self.irrevIDs,self.compass_conditions);
            for i in range(self.compass_expPenalties.rowno):
                rxnd=self.compass_expPenalties.rows[i];
                self.compass_expPenalties.vals[i,:]=rxn2pen[rxnd[:-2]];
        
        
    def minimizeFluxSum(self,rxns,update=True,getDifferentials=True):
        '''Minimizes the sum of absolute value of fluxes of the reactions in <rxns>.

        INPUT
            <rxns>        : List of reactions or a dictionary of reactions as {rxn_id:weight} with rxn IDs including directionality
                            in the case of irreversible models.
            <update>      : Boolean that indicates whether a solution should be generated
                            and incorporated.         
        <getDifferentials>: Boolean that indicates whether shadow prices and reduced costs 
                            should be part of the solution. Set this as False when 
                            mixed-integer programming is involved, as an error will
                            likely be encountered otherwise.
        '''
        if not isinstance(rxns, dict):
            rxns={rxn:1. for rxn in rxns};
        objf_minflux={}; #variables -> coefficients in flux minimization objective
        if self.optimizer.isIrreversible:
            for rxn in rxns:
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    objf_minflux[rxn_dir]=rxns[rxn_dir];
        else:
            Snonneg,Sother=self.__getRxnDirs();
            for rxn in Snonneg.intersection(rxns):
                objf_minflux[rxn]=rxns[rxn];
            for rxn in Sother.intersection(rxns):
                rxn_abs=rxn+'_abs';
                if not rxn_abs in self.optimizer.Dvars:
                    self.optimizer.addAbsVar(rxn);
                objf_minflux[rxn_abs]=rxns[rxn];
        self.__optimize(objf_minflux,direction='min'); 
        if update:
            self.sln=self.__getSolution(getDifferentials=getDifferentials); 
        return self.objective_value
        
    def moma(self,refFluxes,weights={},ignore=[],getDifferentials=True):
        '''Performs Minomization of Metabolic Adjustment (MOMA). A flux distribution that
        minimally deviates from a reference distribution is calculated, using linearized 
        absolute values.

        INPUT
            <refFluxes>   : Dictionary of reference fluxes ({rxn id: flux value}). Deviation from this
                            distribution will be minimized. The minimized variable is the sum of absolute
                            value of differences from the fluxes in the dictionary. Coefficient of every flux
                            in this sum is 1, unless stated otherwise in <weights> (see below).
            <weights>     : Dictionary of coefficients of reaction fluxes in the minimized sum (see description of <refFluxes>).
                            The format is {reaction id: coefficient}.The default value (1) will be overwritten by the numbers 
                            in this dictionary, when applicable. In regular MOMA, all weights are 1, hence this dictionary is empty.
            <ignore>      : Reactions that should be ignored (equivalent to setting their weights to zero in <weights>) during the
                            minimization of the sum of deviations from the original flux distribution. List and list-like containers will work.
        <getDifferentials>: Boolean that indicates whether shadow prices and reduced costs 
                            should be part of the solution. 
        '''
        self.__checkOptimizer();
        objf={}; #variables -> coefficients in moma objective
        for rxn in refFluxes:
            if not rxn in ignore:
                varname=self.optimizer.addMomaVar(rxn,refFluxes[rxn]);
                objf[varname]=1.;
        for rxn in weights:
            rxnvar='_'.join([rxn,'moma']);
            if rxnvar in objf:
                objf[rxnvar]=weights[rxn];
            else:
                print('Reaction '+rxn+' is in weights dictionary but not in flux dictionary or is also in ignored list, so the given weight for this reaction is ignored');
        self.__optimize(objf,direction='min'); 
        self.sln=self.__getSolution(getDifferentials=getDifferentials); 
        deltaFlux={rxn:self.fluxes[rxn]-refFluxes[rxn] for rxn in refFluxes};
        return self.objective_value,deltaFlux
        
    def getRxnWeights(self,exp,specialWeights=None,noGeneWeight=0.1,nanExp=9.,plus=1.):
        '''
        Given expression values (e.g., as TPM) in <expTable> and other parameters, determines weights for every reaction and
        returns as a {rxn id : weight} dictionary.
        
        INPUT
            <exp>             : Dictionary of gene expression values as {gene ID : expression value}. 
                                NaN is permitted as a value. Missing genes (model genes not represented 
                                in the dictionary) will be detected and NaN will be inserted.
            <specialWeights>  : Dictionary of reactions that are assigned weights not conforming to
                                the regular rules. The format of the dictionary is  {rxn : weight}
            <noGeneWeight>    : Weight that will be assigned to reactions without any gene 
                                association (unless indicated in <specialWeights>).
            <nanExp>          : Replaces NaN in weight calculations. This may be due to missing
                                expression values.
            <plus>            : Positive number to be added to every value in the expression table. Make this zero
                                only if this addition was done apriori as zero expression values will
                                cause infinite weights.
        '''
        #quality check and pre=processing 
        if not isinstance(exp,dict):
            raise Exception('Expressions must be provided as a dictionary from genes to values.');
        nothingFound=True;
        for gene in self.genes_id:
            if gene in exp:
                nothingFound=False;
                if np.isnan(exp[gene]):
                    exp[gene]=nanExp+plus;
                else:
                    exp[gene]=exp[gene]+plus;
            else:
                exp[gene]=nanExp+plus;
        if nothingFound:
            raise Exception('None of the genes in model were found in the provided expression table. Please make sure the row indices are gene IDs.');

        if specialWeights:
            Sdiff=set(specialWeights).difference(self.rxns);
            if Sdiff:
                if len(Sdiff)==len(specialWeights):
                    raise Exception('None of the reaction keys in specialWeights dictionary is recognized. Use appropriate reaction IDs from .rxns attribute.');
                else:
                    print(Sdiff);
                    print('WARNING: The above reaction IDs are not part of model reactions and will be ignored.');
        #calculate regular penalties
        rxn2weight={};
        gpr2weight={};
        for rxn in self.rxns:
            if self.rxns[rxn].genes:
                gprrule=self.rxns[rxn].gene_reaction_rule;
                if gprrule in gpr2weight:
                    rxn2weight[rxn]=gpr2weight[gprrule];
                else:
                    gprtree=self.__getGPRtree(rxn);
                    w=self.__traverseGPRtree_weight(gprtree[0],exp)
                    rxn2weight[rxn],gpr2weight[gprrule]=w,w;
            else:
                rxn2weight[rxn]=noGeneWeight;
        #add special weights
        if specialWeights:
            for rxn in specialWeights:
                rxn2weight[rxn]=specialWeights[rxn];
        return rxn2weight

    def __optimizer2dict(self):
        '''Converts the current optimizer object into a dictionary that can be saved for future use.'''
        Dsave={};
        for att in self.getAttributes(self.optimizer):
            if not att=='Dcons' and not att=='Dvars':
                if att=='model':
                    Dsave[att]=self.optimizer.model.to_json();
                else:
                    Dsave[att]=getattr(self.optimizer,att);
        return Dsave
        
    def __setOptimizerFromDict(self,Dsaved):
        '''Recreates the optimizer object from a dictionary that was previously saved.'''
        self.optimizer=opt([]);
        for att in self.getAttributes(self.optimizer):
            if not att=='Dcons' and not att=='Dvars' and not att=='model':
                setattr(self.optimizer,att,Dsaved[att]);
        self.optimizer.model=Model.from_json(Dsaved['model']);
        self.optimizer.Dcons={x.name:x for x in self.optimizer.model.constraints};
        self.optimizer.Dvars={x.name:x for x in self.optimizer.model.variables}; 

    def __getONorOFFrxns(self,genes,ONorOFF):
        '''Uses IMAT rules to determine "ON" or "OFF" reactions as indicated by
        <ONorOFF> and based on the gene list in <genes>.'''
        rxns=set();
        bln={'ON':('True','False',True),'OFF':('False','True',False)}[ONorOFF];
        for gene in genes:
            for rxn in self.genes_id[gene].reactions:
                s=' '+rxn.gene_reaction_rule+' ';
                for x in rxn.genes:
                    if x.id in genes:
                        s=re.sub('([ \(])'+x.id+'([ \)])','\\1'+bln[0]+'\\2',s); #re.sub(' '+x.id+' ',' '+bln[0]+' ',s); #s=re.sub('\\b'+x.id+'\\b',bln[0],s)
                    else:
                        s=re.sub('([ \(])'+x.id+'([ \)])','\\1'+bln[1]+'\\2',s);
                if eval(s.strip())==bln[2]:
                    rxns.add(rxn.id);
        return rxns

    def __getONgenes(self,high,Ron):
        '''Uses IMAT++ rules to determine ON genes given the set of highly expressed 
        genes (<high>) and the set of ON reactions (<Ron>). Returns a dictionary of
        ON genes (as [gene -> ON reactions] for each gene).'''
        Gon={};
        Ronactual=Ron.intersection(self.optimizer.var2binvars);
        for gene in high:
            rxns=[x.id for x in self.genes_id[gene].reactions];
            onrxns=Ronactual.intersection(rxns);
            if onrxns:
                Gon[gene]=onrxns;
        return Gon

    def __checkOptimizer(self):
        '''Subroutine of functions like FBA and pFBA to check if optimizer is set.
        If not, then it is set with default parameters.'''
        if not hasattr(self,'optimizer'):
            print('optimizer has not been set. Setting optimizer with default options.');
            self.setOptimizer();

    def __simpleObj(self,obj,Dir):
        '''
        Sets objective function according to the needs of functions such as fba();
        '''
        if type(obj)==dict:
            obj={self.rxns[x]:obj[x] for x in obj};
        self.model.objective=obj;
        self.model.objective.direction=Dir;
            
        
    def __getInfinity(self):
        '''Determines the large number used in representing infinite flux in reaction boundaries.'''
        V=np.abs([self.rxns[r].lower_bound for r in self.rxns]+[self.rxns[r].upper_bound for r in self.rxns]);
        V.sort();
        if V[-1]==V[-2] and V[-2]==V[-3]:
            self.infinity=V[-1];
        else:
            self.infinity=1000.;
            
    def __getRxnDirs(self):
        '''Subroutine of functions like pFBA(), which use absolute values for reactions
        that may proceed backwards. Returns sets of reactions that may ("Sother") 
        or may not ("Snonneg")  carry negative flux from current optimizer.'''
        Snonneg,Sother=set(),set();
        for rxn in self.rxns:
            if self.optimizer.Dvars[rxn].lb<0:
                Sother.add(rxn);
            else:
                Snonneg.add(rxn);
        return Snonneg,Sother        
        
    def __getSolution(self,getDifferentials=True):
        '''Converts results of current optimization into a solution object and a
        flux discionary.'''
        fluxes,reduced_costs,shadow_prices={},{},{};
        if getDifferentials:
            try:
                _=getattr(self.optimizer.model,'reduced_costs');
            except:
                getDifferentials=False;
        if getDifferentials:
            for rxn in self.optimizer.rxn2varcoefs:
                fluxes[rxn]=0.;
                reduced_costs[rxn]=0.;
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    fluxes[rxn]=fluxes[rxn]+\
                    self.optimizer.rxn2varcoefs[rxn][rxn_dir]*self.optimizer.Dvars[rxn_dir].primal;
                    reduced_costs[rxn]=reduced_costs[rxn]+\
                    self.optimizer.rxn2varcoefs[rxn][rxn_dir]*self.optimizer.Dvars[rxn_dir].dual;
            for met in self.mets:
                shadow_prices[met]=self.optimizer.Dcons[met].dual;
        else:
            for rxn in self.optimizer.rxn2varcoefs:
                fluxes[rxn]=0.;
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    fluxes[rxn]=fluxes[rxn]+\
                    self.optimizer.rxn2varcoefs[rxn][rxn_dir]*self.optimizer.Dvars[rxn_dir].primal;  
        if getDifferentials:
            sln=cobra.Solution(self.objective_value,self.objective_status,\
            myTable.pandas.Series(fluxes.values(),fluxes.keys()),\
            myTable.pandas.Series(reduced_costs.values(),reduced_costs.keys()),\
            myTable.pandas.Series(shadow_prices.values(),shadow_prices.keys()));
        else:
            sln=cobra.Solution(self.objective_value,self.objective_status,\
            myTable.pandas.Series(fluxes.values(),fluxes.keys()));
        self.fluxes=fluxes.copy();        
        return sln

    def __getUniqueMetabolites(self):
        '''Determines unique metabolites based on regular annotation (<metname>_<compartment letter>) or, if this is not possible,
        based on metabolite name. A dictionary (self.uniqmets) is generated which maps each metabolite to its id's and compartments for all
        appearances.'''
        Lrand=np.random.choice(self.model.metabolites,self.sizeRandomChoiceForUniqueMets);
        isRegular=True;
        i=0;
        while isRegular and i<len(Lrand):
            L=re.split('_',Lrand[i].id);
            if not L[-1]==Lrand[i].compartment:
                isRegular=False;
            i=i+1;
        if isRegular:
            print('Compartmentalization convention is detected in metabolite IDs (<metabolite id>_<compartment letter>). Unique metabolites will be named according to ID without the compartment suffix.')
            namefnc=lambda x:x.id[:-(len(x.compartment)+1)];
        else:
            print('Metabolite IDs do not follow standard annotation for compartments. Unique metabolites will be named after metabolite names rather than IDs.') 
            namefnc=lambda x:x.name;
        Duniqmets={};
        for met in self.model.metabolites:
            name=namefnc(met);
            if name in Duniqmets:
                Duniqmets[name].append((met.id,met.compartment));
            else:
                Duniqmets[name]=[(met.id,met.compartment)];
        return Duniqmets              
            
    def __optimize(self,objf,direction='max',conv=1.):
        '''Subroutine for functions that perform optimization. Optimizes the problem 
        in the optimizer and sets objective-related attributes. Objective function
        in dictionary form (<objf>) and a sign converter (<conv>) need to be provided.'''
        self.optimizer.setObjective(objf,direction=direction);
        self.objective_function=objf.copy();
        self.objective_status=self.optimizer.optimize();
        try:
            self.objective_value=conv*self.optimizer.model.objective.value;  
        except:
            self.objective_value=np.nan;
            
    def __getTable(self,objects,attributes,annotations):
        '''Makes an attributes table for given <objects> and <attributes> lists.'''
        data={};
        Lid=[];
        if not annotations:
            annotations=[];
        elif type(annotations)==str:
            annotations=[annotations];
        for att in attributes + annotations:
            data[att]=[];
        for obj in objects:
            for att in data:
                if att=='flux' and self.sln:
                    data[att].append('{:.6f}'.format(self.sln.fluxes[getattr(obj,'id')]));
                elif att=='flux':
                    try:
                        data[att].append('{:.6f}'.format(getattr(obj,att)));
                    except:
                        data[att].append(np.nan);
                elif att in annotations and hasattr(obj,'annotation'):
                    d=getattr(obj,'annotation');
                    if att in d:
                        data[att].append(d[att]);
                    else:
                        data[att].append(np.nan);
                elif hasattr(obj,att):
                    data[att].append(getattr(obj,att));
                else:
                    data[att].append(np.nan);
            Lid.append(obj.id);
        df=pd.DataFrame(data,index=Lid);
        df.dropna(how='all', axis=1, inplace=True);
        return df

    def __getLatentRxns(self,Rsuper,flux_thre):
        '''Determines latent rections, which are defined as reactions in the given
        <Rsuper> set of highly expressed reactions whose reactants (or products, if the reaction is reversible)
        consist of metabolites that all carry flux based on a threshold of total (metabolite-level) flux
        sum designated as <flux_thre>. Since directinality is important, and since this subroutine is used to
        activate ON reactions in Rsuper, latent reactions are returned as a list of IMAT variables in relevant
        directions (rxnid_yf for forward and rxnid_yr for reverse).'''
        activemets,inactivemets,latentrxns=set(),set(),set();
        for rxn in Rsuper:
            r=self.rxns[rxn];
            if abs(self.optimizer.Dvars[r.id].primal)<self.tol:
                if r.upper_bound>0. and r.reactants:
                    islatent,n,k=True,len(r.reactants),0;
                    while k<n and islatent:
                        m=r.reactants[k];
                        if m.id in inactivemets:
                            islatent=False;
                        elif not m.id in activemets:
                            if sum([self.optimizer.Dvars[x.id].primal for x in m.reactions])>=flux_thre:
                                activemets.add(m.id);
                            else:
                                inactivemets.add(m.id);
                                islatent=False
                        k=k+1
                    if islatent:
                        latentrxns.add(r.id+'_yf');
                if r.lower_bound<0. and r.products:
                    islatent,n,k=True,len(r.products),0;
                    while k<n and islatent:
                        m=r.products[k];
                        if m.id in inactivemets:
                            islatent=False;
                        elif not m.id in activemets:
                            if sum([self.optimizer.Dvars[x.id].primal for x in m.reactions])>=flux_thre:
                                activemets.add(m.id);
                            else:
                                inactivemets.add(m.id);
                                islatent=False
                        k=k+1;
                    if islatent:
                        latentrxns.add(r.id+'_yr');
        return latentrxns
        
    def __removeIMATConstraints(self):
        '''Removes a constraints related to IMAT applications ("imat" or "imatplusplus_i", where i=1-4) 
        from current optimizer.'''
        Lcons=[con for con in self.optimizer.Dcons];
        for con in Lcons:
            if re.search('\\bimat(plusplus_[1-4]){0,1}\\b',con):
                self.optimizer.model.remove(self.optimizer.Dcons[con]);
                _=self.optimizer.Dcons.pop(con);     

    def __setExpressionPenalties(self,expTable,specialPenalties,noGenePenalty,nanPenalty,plus,verbose):
        '''
        Subroutine of setFPA(). See setFPA() instructions related to expression.               
        '''
        #quality check of the table and modifications  
        if isinstance(expTable,pd.DataFrame):
            self.fpa_expTable=matrix(expTable.to_numpy(),list(expTable.index),list(expTable.columns));
        elif isinstance(expTable,matrix):
            self.fpa_expTable=expTable.copy();
        else:            
            raise Exception('Expressions must be provided either as a dataframe or as a matrix (not numpy) object');
        if not np.issubdtype(self.fpa_expTable.vals.dtype,np.number):
            raise Exception('Data is not numeric. The values of the expression table should constitute a numerical numpy array.');
        #conditions
        if 'super_cond' in self.fpa_expTable.cols:
            raise Exception('One of the conditions is named "super_cond". You cannot use this name as it is taken to identify the super condition of FPA.');
        self.fpa_conditions=[col for col in self.fpa_expTable.cols];
        Ncond=len(self.fpa_conditions);
        #genes
        missingGenes=set(self.genes_id).difference(self.fpa_expTable.rows);
        Nmiss=len(missingGenes);
        if Nmiss==len(self.genes):
            raise Exception('None of the genes in model were found in the provided expression table. Please make sure the row indices are gene IDs.');
        #special penalties
        Sdiff=set(specialPenalties).difference(self.irrevIDs);
        if Sdiff:
            if len(Sdiff)==len(specialPenalties):
                raise Exception('None of the reaction keys in specialPenalties dictionary is recognized. Use appropriate IDs from .irrevIDs attribute. Note that rxn keys must be directional, meaning, reaciton IDs followed by a "_f" or "_r" suffix');
            else:
                print(Sdiff);
                print('WARNING: The above reaction IDs are not part of irrevIDs list and will be ignored.');
        #value addition
        if plus>0:
            self.fpa_expTable.vals=self.fpa_expTable.vals+plus;
        elif plus<0 and verbose:
            print('WARNING: the variable <plus> must be a positive number or zero. The given number is ignored.')        
        #convert missing gene values to nans
        if Nmiss:
            nans=np.empty((Nmiss,self.fpa_expTable.colno));
            nans[:]=np.nan;            
            self.fpa_expTable=self.fpa_expTable.appendRows(nans,list(missingGenes));
        #calculate regular penalties
        rxn2pen={};
        gpr2pen={};
        for rxn in self.rxns:
            if self.rxns[rxn].genes:
                gprrule=self.rxns[rxn].gene_reaction_rule;
                if gprrule in gpr2pen:
                    rxn2pen[rxn]=gpr2pen[gprrule];
                else:
                    c=np.nan_to_num(self.__getRxnPenalty(rxn),nan=nanPenalty);
                    rxn2pen[rxn],gpr2pen[gprrule]=c,c;
            else:
                rxn2pen[rxn]=Ncond*[noGenePenalty];
        dummyvals=np.empty((len(self.irrevIDs),len(self.fpa_conditions)));
        self.fpa_expPenalties=matrix(dummyvals,self.irrevIDs,self.fpa_conditions);
        for i in range(self.fpa_expPenalties.rowno):
            rxnd=self.fpa_expPenalties.rows[i];
            self.fpa_expPenalties.vals[i,:]=rxn2pen[rxnd[:-2]];
        #Keep the original penalties for the record
        self.fpa_originalExpPenalties=self.fpa_expPenalties.copy();
        #add special penalties
        rxnd2pen_super={}; #holds special rules regarding super condition
        for rxnd in specialPenalties:
            if type(specialPenalties[rxnd])==dict:
                for cond in specialPenalties[rxnd]:
                    self.fpa_expPenalties[rxnd,cond]=float(specialPenalties[rxnd][cond]);
                rxnd2pen_super[rxnd]=min(1.,np.min(self.fpa_expPenalties.getRow(rxnd))); #super condition rule
            else:
                val=float(specialPenalties[rxnd]);
                for cond in self.fpa_conditions:
                    self.fpa_expPenalties[rxnd,cond]=val;
                rxnd2pen_super[rxnd]=val;
        #add super condition
        self.fpa_expPenalties=self.fpa_expPenalties.appendCols(np.ones([self.fpa_expPenalties.rowno,1]),['super_cond']);
        for rxnd in rxnd2pen_super:
            self.fpa_expPenalties[rxnd,'super_cond']=rxnd2pen_super[rxnd];
                                
    def __setDistanceCalculator(self,byproducts,inherit_array,inherit_rows,inherit_cols,nandist,setDistanceNet):
        '''
        Subroutine of setFPA(). See setFPA() instructions related to distance.                                           
        '''
        #check and set byproducts
        if byproducts:
            byproducts=self.__getByProducts(byproducts);
        else:
            byproducts=[];
        #establish the reaction network
        if setDistanceNet:
            S=self.getSmatrix();
            Net=establishRxnNetwork(S.vals,S.cols,S.rows,[self.rxns[x].lower_bound for x in S.cols],\
                                    [self.rxns[x].upper_bound for x in S.cols],byproducts);
        else:
            Net=None;
        #check and set the inherited distance calculations 
        if inherit_array:
            #check and allocate
            if isinstance(inherit_array,matrix):
                Srows,Scols=set(inherit_array.rows),set(inherit_array.cols);
            elif isinstance(inherit_array,pd.DataFrame):
                Srows,Scols,vals=set(inherit_array.index),set(inherit_array.columns),inherit_array.to_numpy();
            elif isinstance(inherit_array,np.ndarray):
                try:
                    Srows,Scols=set(inherit_rows),set(inherit_cols);
                except:
                    raise Exception('<inherit_rows> and <inherit_cols> should both be provided as lists to define the dimensions of <inherit_array> numpy matrix.');
                vals=inherit_array.copy();
            else:
                raise Exception('<inherit_array> must be in the form of either matrix (from myCobra), or pandas DataFrame, or numpy.ndarray.');
            for rxn_dir in Scols:
                if not rxn_dir in self.irrevIDs: #previously --if not rxn_dir in Net.Drxn2node--
                    rxn=rxn_dir[:-2];
                    if not rxn in self.rxns:
                        raise Exception('The reaction ID root (all letters except the last two) of '+str(rxn_dir)+' in columns of <inherit_array> is not recognized as a model reaction.');
                    else: 
                        _dir=rxn_dir[-2:];
                        if not _dir=='_f' and not _dir=='_r':
                            raise Exception(str(rxn_dir)+' in columns of <inherit_array> does not indicate a valid direction for the pertaining reaction. The last two letters of reaction identifiers must be either "_f" (forward direction) or "_r" (reverse direction).')
            self.__checkIrrevIDsInRows(Srows,'<inherit_array>');
            #set
            if isinstance(inherit_array,pd.DataFrame):
                self.fpa_calculatedDist=matrix(np.nan_to_num(vals,nan=nandist),list(inherit_array.index),list(inherit_array.columns));
            elif isinstance(inherit_array,np.ndarray):
                self.fpa_calculatedDist=matrix(np.nan_to_num(vals,nan=nandist),inherit_rows,inherit_cols);
            else:
                self.fpa_calculatedDist=matrix(np.nan_to_num(inherit_array.vals,nan=nandist),inherit_array.rows,inherit_array.cols);
        else:
            if Net is None:
                raise Exception('If <setDistanceNet> is False, then an <inherit_array> must be provided.')
            self.fpa_calculatedDist=matrix([],[],[]); #empty matrix that allows row search in related code
        #set the remaining variables related to metabolic distance calculations
        self.fpa_distNet,self.fpa_byproducts,self.fpa_nandist=Net,byproducts,nandist;
        
    def __setConditionNetworks(self,networkTable,mergeNetworks):
        '''
        Subroutine of setFPA(). See setFPA() instructions related to network.                                        
        '''
        if not networkTable is None:
            if isinstance(networkTable,pd.DataFrame):
                rows,cols,vals=list(networkTable.index),list(networkTable.columns),networkTable.to_numpy();
            elif isinstance(networkTable,matrix):
                rows,cols,vals=networkTable.rows,networkTable.cols,networkTable.vals;        
            #check rows and cols of given table
            self.__checkIrrevIDsInRows(rows,'<networkTable>');
            if not set(cols)==set(self.fpa_conditions):
                print('You provided the following conditions in network table:\n'+str(networkTable.cols));
                print('But the following conditions are given in expression table:\n'+str(self.fpa_conditions));
                raise Exception('conditions provided for network table and expression table should be the same.');
            #check values and form a matrix object
            if not ((vals==0) | (vals==1)).all():
                raise Exception('Some values of the network Table are not binary (0 or 1). Only zeros and ones allowed in this input.');
            if mergeNetworks:
                for i in range(vals.shape[0]):
                    vals[i,:]=np.max(vals[i,:]);
        else:
            rows,cols=[x for x in self.irrevIDs],[x for x in self.fpa_conditions];
            vals=np.ones((len(rows),len(cols)));
        self.fpa_rxnNet=matrix(vals,rows,cols);
        #add super condition
        self.fpa_rxnNet=self.fpa_rxnNet.appendCols(np.ones([self.fpa_rxnNet.rowno,1]),['super_cond']);

    def __appendCol2DistanceMatrix(self,target):
        '''Adds column of .fpa_nandist values with title <target> to the calculated distance array (.calculatedDist).'''
        distcol=self.fpa_nandist*np.ones([len(self.irrevIDs),1]);  
        if self.fpa_calculatedDist.isempty():
            self.fpa_calculatedDist=matrix(distcol,self.irrevIDs,[target]);
        else:
            self.fpa_calculatedDist=self.fpa_calculatedDist.appendCols(distcol,[target]);
        
    def __getMetaboliteList(self,given,errorname):
        '''Converts <given> into a valid list of model metabolites. <given> can be None,
        a string, or a list. Both unique and specific metabolite names are covered. 
        If a given name is not found, an exception will be thrown using <errorname> as identifier.'''
        L,Lbad=[],[];
        if given:
            if type(given)==str:
                given=[given];
            for k in given:
                if k in self.uniqmets:
                    L=L+[x[0] for x in self.uniqmets[k]];
                elif k in self.mets:
                    L.append(k);
                else:
                    raise Exception(errorname+': Metabolite '+str(k)+' is not found');
        return L        

    def __blockMetabolite(self,met,blockEntry,blockExit,exclude=[]):
        '''Blocks metabolite <met> from entry (if <blockEntry> is True) or exit (if blockExit is True)
        or both (if both true). This is done by properly constraining demand, sink and exchange
        reactions associated with the metabolite. Reactions indicated by <exclude> 
        (as string for one reaction or as a list [or the like] of reactions as IDs) are omitted.'''
        if type(exclude)==str:
            exclude=set([exclude]);
        else:
            exclude=set(exclude);
        Srxns=set([]);
        for rxn in self.mets[met].reactions:
            if rxn in self.model.exchanges or rxn in self.model.demands or rxn in self.model.sinks:
                Srxns.add(rxn.id);
        Srxns=Srxns.difference(exclude);
        if Srxns:
            Dori=self.__getCurrentBoundaries(Srxns);
            if blockExit:
                for rxn in Srxns:
                    if self.rxns[rxn].metabolites[self.mets[met]]<0:
                        _=self.changeRxnBoundaries({rxn:[min(0,Dori[rxn][0]),0.]});
                    else:
                        _=self.changeRxnBoundaries({rxn:[0,max(0,Dori[rxn][1])]});
            if blockEntry:
                for rxn in Srxns:
                    if self.rxns[rxn].metabolites[self.mets[met]]<0:
                        _=self.changeRxnBoundaries({rxn:[0,max(0,Dori[rxn][1])]});
                    else:
                        _=self.changeRxnBoundaries({rxn:[min(0,Dori[rxn][0]),0.]});
        else:
            Dori={};
        return Dori
        
    def __getRxnPenalty(self,rxn):
        '''Generates and returns expression penalty for reaction named <rxn>.'''
        genes=[gene.id for gene in getattr(self.model.reactions,rxn).genes];
        exp={gene:self.fpa_expTable.getRow(gene) for gene in genes};
        gprtree=self.__getGPRtree(rxn);   
        _,c=self.__traverseGPRtree(gprtree[0],exp);
        return c        
   
    def __getGPRtree(self,rxn):
        '''Converts GPR rule of <rxn> to a list of nodes in gprnodes class.This 
        list is returned.'''
        s=' '+self.rxns[rxn].gene_reaction_rule+' ';
        #preparation of s for parsing
        for gene in self.rxns[rxn].genes:
            s=re.sub('([ \(])'+gene.id+'([ \)])','\\1('+gene.id+')\\2',s);
        for t in  [(' or ',' | '),(' and ',' & ')]:
            s=re.sub(t[0],t[1],s);
        s=s.strip();
        #definition of nodes by parsing s
        root=gprnode(0,None,[]);
        Lnodes,current,i=[root],0,0;
        word='';
        while i<len(s):
            x=s[i];
            if x=='(':
                Lnodes.append(gprnode(len(Lnodes),Lnodes[current],[]));
                Lnodes[current].appendChild(Lnodes[-1]);
                current=Lnodes[-1].index;
            elif x==')':
                if word:
                    Lnodes[current].setGene(word);
                    word='';
                current=Lnodes[current].parent.index;
            elif x=='|' or x=='&':
                if Lnodes[current].connector is None:
                    Lnodes[current].setConnector(x);
                elif not x==Lnodes[current].connector:
                    print('WARNING: Bad GPR detected for rxn '+rxn+'. The switch in connector from '+Lnodes[current].connector+' to '+x+' is ignored.');                    
                    #raise Exception('Inconsistent connectors in the GPR rule '+s);
            elif not x==' ':
                word=word+x;
            i=i+1;
        return Lnodes
        
    def __checkIrrevIDsInRows(self,rows,sobject):
            Sdiff=set(self.irrevIDs).difference(rows);
            if Sdiff:
                print(Sdiff);
                raise Exception('The row IDs in '+sobject+' do not cover the above model reactions with directionality. The rows of this object must cover all irreversible reactions indicated by .irrevIDs attribute.');

    def __getIrrevRxnIDs(self):
        '''Returns reaction IDs with directionality determined based on lower boundaries.
        All reactions have a forward direction ("_f") suffix and some (set as reversibility=True
        or lower_boundary<0) also have a reverse direction ("_r" suffix).'''
        L=[];
        for rxn in self.rxns:
            L.append(rxn+'_f');
            if self.rxns[rxn].reversibility or self.rxns[rxn].lower_bound<0:
                L.append(rxn+'_r');
        return L
        
    def __getByProducts(self,Lin):
        '''Returns the list of metabolite IDs that match unique metabolite names or IDs in <Lin>.'''
        Lout,ignored=[],[];
        for met in Lin:
            if met in self.uniqmets:
                for t in self.uniqmets[met]:
                    Lout.append(t[0]);
            elif met in self.mets:
                Lout.append(met);
            else:
                ignored.append(met);
        if ignored:
            print('The following metabolite names in byproducts were not found and ignored:\n'+str(ignored));
        return Lout
        
    def __getCurrentBoundaries(self,rxns):
        '''Returns current boundaries (i.e., those in optimizer) for reactions in <rxns>.'''
        D={};        
        if not self.optimizer.isIrreversible:
            for rxn in rxns:
                D[rxn]=[self.optimizer.Dvars[rxn].lb,self.optimizer.Dvars[rxn].ub];
        else:
            for rxn in rxns:
                L=[0,0,0,0];
                for rxn_dir in self.optimizer.rxn2varcoefs[rxn]:
                    if rxn_dir[-1]=='f':
                        L[2],L[3]=[self.optimizer.Dvars[rxn_dir].lb,self.optimizer.Dvars[rxn_dir].ub];
                    else:
                        L[0],L[1]=[-self.optimizer.Dvars[rxn_dir].ub,-self.optimizer.Dvars[rxn_dir].lb];
                if L[2] and L[3]:
                    D[rxn]=[L[2],L[3]];
                elif L[0] and L[1]:
                    D[rxn]=[L[0],L[1]];
                elif L.count(0)>=2: #!can be changed to "else" if this is certain
                    D[rxn]=[L[0],L[3]];
                else:
                    raise Exception('Unexpected boundaries encountered for '+rxn+'. Code of __getCurrentBoundaries needs to be checked. Please report this incidence.');
        return D
        

    @staticmethod
    def __extendTable(table,objects,attributes,refcol,position):
        '''Subroutine of functions like extendRxnTable().'''
        if type(attributes)==str:
            attributes=[attributes];
        for att in attributes:
            vals=[getattr(getattr(objects,obj),att) for obj in table.index];
            table.insertColumn(att,vals,refcol,position=position);
            refcol=att;

    @staticmethod
    def __traverseGPRtree(node,exp):
        '''Traverses a gpr tree starting at gprnode <node> to derive expression penalty
        for FPA based on gene expression values in <exp> and genes in the bottom nodes.
        Expression sum and expression penalty are returned in the respective order.'''
        if node.connector:
            L=[mnm.__traverseGPRtree(x,exp) for x in node.children];
            S=np.stack([t[0] for t in L]);            
            C=np.stack([t[1] for t in L]);
            if node.connector=='&':
                s=np.nanmin(S,axis=0);
                c=np.nanmax(C,axis=0);                
            elif node.connector=='|':
                s=np.nansum(S,axis=0);
                c=np.nanmax(s)/s;
            else:
                raise Exception('ERROR: unknown connector.')
        elif node.gene:
            s=exp[node.gene];
            c=np.nanmax(s)/s;
        else: #the root in the single gene case
            if len(node.children)>1:
                raise Exception('ERROR: bad GPR node with no connectors and multiple children'); 
            s,c=mnm.__traverseGPRtree(node.children[0],exp);                
        return s,c
        
    @staticmethod        
    def __traverseGPRtree_perturbation(node,Dexp):
        '''Traverses a gpr tree starting at gprnode <node> to derive reaction expression
        based on dual expression values in <Dexp> which reflect expression of genes
        in unperturbed and perturbed conditions, respectively ({gene : (expression,perturbed expression)}).
        A tuple of reaction expression for unperturbed and perturbed cases is returned.'''
        if node.connector:
            L=[mnm.__traverseGPRtree_perturbation(x,Dexp) for x in node.children];
            if node.connector=='&':
                t=(min([x[0] for x in L]),min([x[1] for x in L]));                
            elif node.connector=='|':
                t=(sum([x[0] for x in L]),sum([x[1] for x in L]));
            else:
                raise Exception('ERROR: unknown connector.')
        elif node.gene:
            t=Dexp[node.gene];
        else: #the root in the single gene case
            if len(node.children)>1:
                raise Exception('ERROR: bad GPR node with no connectors and multiple children'); 
            t=mnm.__traverseGPRtree_perturbation(node.children[0],Dexp);                
        return t
        
    @staticmethod
    def __traverseGPRtree_weight(node,exp):
        '''Traverses a gpr tree starting at gprnode <node> to derive reaction weights
        based on gene expression values in <exp> and genes in the bottom nodes.
        Weight is returned.'''
        if node.connector:
            L=[mnm.__traverseGPRtree_weight(x,exp) for x in node.children];
            if node.connector=='&':
                w=np.nanmax(L);             
            elif node.connector=='|':
                w=np.nanmin(L);
            else:
                raise Exception('ERROR: unknown connector.')
        elif node.gene:
            w=1./exp[node.gene];
        else: #the root in the single gene case
            if len(node.children)>1:
                raise Exception('ERROR: bad GPR node with no connectors and multiple children'); 
            w=mnm.__traverseGPRtree_weight(node.children[0],exp);                
        return w
            
    @staticmethod
    def isnum(x):
        '''Returns whether the input varible is numerical or not.'''
        try:
            _=float(x);
            return True
        except:
            return False
            
    @staticmethod
    def sv(obj,filename):
        '''Pickles and saves <obj> in <filename>.'''
        fo=open(filename,'wb');
        pickle.dump(obj, fo, pickle.HIGHEST_PROTOCOL)
        fo.close();
    
    @staticmethod    
    def op(filename):
        '''Opens the object pickled and saved in <filename>.'''
        fo=open(filename, 'rb');
        obj=pickle.load(fo);
        fo.close();
        return obj

    @staticmethod    
    def getAttributes(obj):
        '''Returns all attributes of an object, excluding methods.'''
        return {name: attr for name, attr in obj.__dict__.items()
                if not name.startswith("__")
                and not callable(attr)
                and not type(attr) is staticmethod}
                



if False:
    M=mnm('../cobra/iCEL1314.json','json');
    d=M.getPerturbedRxns(Dn2w['sams-1'],0.,Dexp);
    print(d);

if False:
    M=mnm('../cobra/iCEL1314.json','json');
    fo=open('/data/yilmazl/singlecell/temp2.txt','r');
    Dref=eval(fo.read());
    fo.close();
    fo=open('/data/yilmazl/singlecell/temp3.txt','r');
    Dref2=eval(fo.read());
    fo.close();
    M.changeRxnBoundaries({'RM00479':[0.,0.]})
    x,dF=M.moma(Dref);
    print(x);
    x,dF=M.moma(Dref2);  
    print(x);
    M.pFBA('RM04432');
    x,dF=M.moma(Dref);
    print(x);

if False:
    #M=mnm('e_coli_core.json','json');
    M=mnm('../cobra/iCEL1314.json','json');
    M.pFBA('RCC0005');
    print(M.fluxes['BIO0010']);
    D,F={},{};
    for rxn in M.rxns:
        D[rxn]=[M.optimizer.Dvars[rxn].lb,M.optimizer.Dvars[rxn].ub];
        F[rxn]=M.fluxes[rxn];
    fo=open('/data/yilmazl/singlecell/temp.txt','w');
    fo2=open('/data/yilmazl/singlecell/temp3.txt','w');
    fo.write(str(D));
    fo2.write(str(F));
    fo.close();
    fo2.close();
    
if False:
    rows=[x for x in M.irrevIDs];
    cols=['a','b','c','d','e'];
    df=pd.DataFrame(np.random.random([len(rows),len(cols)]),index=rows,columns=cols);
    mat=matrix(df.to_numpy(),rows,cols);
    
    to=time.time();
    df2=pd.DataFrame(np.zeros([len(rows),len(cols)]),index=rows,columns=cols);
    for row in df.index:
        for col in df.columns:
            df2.loc[row,col]=df.loc[row,col]**0.5;
    print(time.time()-to);
    to=time.time();
    mat2=matrix(np.zeros([len(rows),len(cols)]),rows,cols);
    for row in mat.rows:
        for col in mat.cols:
            mat[row,col]=mat[row,col]**0.5;
    print(time.time()-to);

        
#s='(abc) | (((def) & (egf)) & ((hij) | (klm)) & (nop))' 


#fnc(Lnodes[3])        

#M=mnm('../cobra/iCEL1314.json','json')
#print(M.getMW('pyr_c')) 

if False:      
    if True:
        M=mnm('iCEL1314.json','json');
        #M.setOptimizer(irreversible=False);
        #R=mnm('Recon3DModel_301.json','json');
        #M.FBA('BIO0010');
    #T=R.search('beta-Alanine','mets')
    
    t1=time.time();
    ysum,vmin_low,vmin_total,summary=M.IMATplusplus([M.genes[x].id for x in ['acdh-1','bcat-1','Y44A6D.5','pcca-1','bas-1','amx-2']],[M.genes[x].id for x in ['sdhb-1','hphd-1']],[M.genes[x].id for x in ['icl-1']],\
    relaxedHigh=True,minimizeFluxSum=True,rescueLatent=True);
    #ysum,vmin_total,summary=M.IMAT([M.genes[x].id for x in ['acdh-1','hphd-1']],[M.genes[x].id for x in ['pcca-1']],relaxedHigh=True,minimizeFluxSum=True,constrain=True);
    
    
    #print(M.FVA('RM04432'));
    #M.optimizer.model.objective=Objective(M.optimizer.Dvars['RC01061'], direction='min');
    #M.optimizer.model.optimize();
    if False:
        tpl=M.pFBA('BIO0010','max',leaveConstraint=True);
        print(tpl);
    
    t2=time.time();
        
    print(t2-t1);    
    print(time.time()-to); 
    
    
    if False:
        M.setConstraints(loadobj('conEco_LB'),True);
        M.fba('BIOMASS_Ec_iJO1366_core_53p95M',parsimonious=True);
    
    if False:
        Mo=mnm('iCEL1314.json','json');
        M=mnm('iCEL1314.json','json');
        Mo.fba('RCC0005',parsimonious=False);
        M.fba('RCC0005',parsimonious=True);
        #T=M.getFluxes(0.5,1.5,absolute=True)
#print(time.time()-to);             