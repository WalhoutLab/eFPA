import h5py
from scipy.sparse import csr_matrix
from myCobra import matrix as Mat
import pandas as pd
import numpy as np

class hdf:
    '''This class is designed to process HDF5 file objects, facilitating the management of sparse matrices and data exploration.'''    
    def __init__(self, filePath, setPaths=True, setSparse=True, dataMap={'data':'X/data', 'indices':'X/indices', 'indptr':'X/indptr'}):
        '''Initializes an hdf5 file located at <filePath>. Constructs objects for "data", "indices", and "indptr" based on specified <dataMap>.
        
        - <setPaths>: When True, generates a list of object paths within the file.
        - <setSparse>: When True, constructs a sparse matrix from the data provided.
        - Follow guidance on argument adjustments if initial inputs are incompatible with the specified h5 file.
        '''
        self.hdf_file = h5py.File(filePath, 'r')
        self.setData(dataMap)
        if setPaths:
            self.setPaths()
        if setSparse:
            self.setSparse()
    
    def setPaths(self):
        '''Generates a list of all paths within the HDF5 file.'''
        self.paths = []
        self.hdf_file.visit(self.paths.append)
        
    def setData(self, dataMap):
        '''Sets data attributes based on the provided <dataMap> dictionary.'''
        try:
            for k in dataMap:
                setattr(self, k, self.hdf_file[dataMap[k]][:])
        except KeyError as e:
            print(f'KeyError: {e} The specified key was not found in the HDF5 file. Check your dataMap keys.')

    def setObs(self,path,catPath=None):
        '''Sets observations based on object of hdf file in <path>. If this is only codes for a category, then the
        path to the categories must also be specified using <catPath>'''
        self.obs=self.__getList(path,catPath)
        self.obs_index={self.obs[i]:i for i in range(self.Nobs)}

    def setVar(self,path,catPath=None):
        '''Sets variables based on object of hdf file in <path>. If this is only codes for a category, then the
        path to the categories must also be specified using <catPath>'''
        self.var=self.__getList(path,catPath)
        self.__checkSparse()
        self.var_index={self.var[i]:i for i in range(self.Nvar)}
            
    def setSparse(self):
        '''Generates a sparse matrix using the "data", "indices", and "indptr" attributes. Requires prior execution of setData.'''
        if not hasattr(self, 'data'):
            print('Warning: Data object not found. Ensure data is set using setData before initializing sparse matrix construction.')
        else:
            shape = (self.Nobs, self.Nvar) if (self.Nobs and self.Nvar) else None
            self.sparse = csr_matrix((self.data, self.indices, self.indptr), shape=shape)
            
    def setObsTable(self,pathMap,indexName=None,pfx='obs/'):
        '''Generates a dataframe of observations using the getDataFrame() function. See getDataFrame for guidance. 
        If obs attribute is available and <indexName> is not provided,it will be used as index. In reverse, 
        if obs has not yet been formed but <indexName> is provided, obs will be formed as this index.'''
        self.obsTable=self.getDataFrame(pathMap,indexName,pfx)
        if not indexName and hasattr(self,'obs'):
            self.obsTable.index=self.obs
        if indexName and not hasattr(self,'obs'):
            self.obs=self.obsTable.index.tolist()
            self.obs_index={self.obs[i]:i for i in range(self.Nobs)}

    def setVarTable(self,pathMap,indexName=None,pfx='var/'):
        '''Generates data frame of variables using the getDataFrame() function. See getDataFrame for guidance.
        If var attribute is available and <indexName> is not provided,it will be used as index. In reverse, 
        if var has not yet been formed but <indexName> is provided, var will be formed as this index.'''
        self.varTable=self.getDataFrame(pathMap,indexName,pfx)  
        if not indexName and hasattr(self,'var'):
            self.varTable.index=self.var
        if indexName and not hasattr(self,'var'):
            self.var=self.varTable.index.tolist()
            self.var_index={self.var[i]:i for i in range(self.Nvar)}
        self.__checkSparse()
 
    def setObsMask(self, obsList=None, cols=None, vals=None, addBy=None):
        """
        Generates a mask for the observations DataFrame (obsTable), if available. This method can be based
        on a list of observation values (<obsList>) or on column names (<cols>) and a list of value tuples (<vals>)
        for nested masking on the obsTable DataFrame.

        Parameters:
        - obsList: List of observation identifiers.
        - cols: List of column names from the obs DataFrame for masking.
        - vals: List of value tuples, each matching the length of <cols>, providing values for masking.
        - addBy: Operator for combining with a pre-existing mask, if any. Accepts "|" (OR) or "&" (AND).

        If <obsList> is provided, <cols> and <vals> are ignored. If both <cols> and <vals> are provided,
        they are used to generate a mask via nested masking logic (see getDataFrameMask for details).
        The <addBy> parameter determines how this new mask is combined with any existing mask.
        """
        if obsList:
            if cols or vals:
                print('Note: <obsList> is provided; <cols> and <vals> will be ignored.')
            if not hasattr(self, 'obs'):
                raise Exception('Error: obs attribute not found. Please set observations first.')
            obs_mask = np.zeros(self.Nobs, dtype=bool)
            for obs in obsList:
                obs_mask[self.obs_index[obs]] = True
        elif cols and vals:
            if not hasattr(self, 'obsTable'):
                raise Exception('Error: obsTable attribute not found. Please set the obsTable first.')
            obs_mask = self.getDataFrameMask(self.obsTable, cols, vals)
        else:
            raise Exception('Error: Insufficient input. Provide either <obsList> or both <cols> and <vals>.')

        # Combine the new mask with the existing mask if specified
        if addBy and hasattr(self, 'obs_mask'):
            if addBy == '|':
                self.obs_mask |= obs_mask.copy()
            elif addBy == '&':
                self.obs_mask &= obs_mask.copy()
            else:
                print('Warning: <addBy> value not recognized ("|" or "&" expected). The obs_mask will be overwritten.')
                self.obs_mask = obs_mask.copy()
        else:
            self.obs_mask = obs_mask.copy()
            
    def setVarMask(self, varList=None, cols=None, vals=None, addBy=None):
        """
        Generates a mask for the variables DataFrame (varTable), if available. This method can be based
        on a list of variable values (<varList>) or on column names (<cols>) and a list of value tuples (<vals>)
        for nested masking on the varTable DataFrame.

        Parameters:
        - varList: List of variable identifiers.
        - cols: List of column names from the var DataFrame for masking.
        - vals: List of value tuples, each matching the length of <cols>, providing values for masking.
        - addBy: Operator for combining with a pre-existing mask, if any. Accepts "|" (OR) or "&" (AND).

        If <varList> is provided, <cols> and <vals> are ignored. If both <cols> and <vals> are provided,
        they are used to generate a mask via nested masking logic (see getDataFrameMask for details).
        The <addBy> parameter determines how this new mask is combined with any existing mask.
        """
        if varList:
            if cols or vals:
                print('Note: <varList> is provided; <cols> and <vals> will be ignored.')
            if not hasattr(self, 'var'):
                raise Exception('Error: var attribute not found. Please set variables first.')
            var_mask = np.zeros(self.Nvar, dtype=bool)
            for var in varList:
                var_mask[self.var_index[var]] = True
        elif cols and vals:
            if not hasattr(self, 'varTable'):
                raise Exception('Error: varTable attribute not found. Please set the varTable first.')
            var_mask = self.getDataFrameMask(self.varTable, cols, vals)
        else:
            raise Exception('Error: Insufficient input. Provide either <varList> or both <cols> and <vals>.')

        # Combine the new mask with the existing mask if specified
        if addBy and hasattr(self, 'var_mask'):
            if addBy == '|':
                self.var_mask |= var_mask.copy()
            elif addBy == '&':
                self.var_mask &= var_mask.copy()
            else:
                print('Warning: <addBy> value not recognized ("|" or "&" expected). The var_mask will be overwritten.')
                self.var_mask = var_mask.copy()
        else:
            self.var_mask = var_mask.copy()        
    
    def getDataFrame(self,pathMap,indexName=None,pfx=''):
        '''Generates and returns a data frame of keys in <pathMap> with key names as columns. 
        
        - <pathMap>: A dictionary of the form {name: (path to values, path to categories)}, where path to values
        indicates the path of the variable to be listed, and if available, path to categories indicate the categorical
        data for an indexed variabled. If not available, path to categories should be None.
        - <indexPath>: If indicated, this column (name in <pathMap>) will be used as the index of the data frame
        - <pfx>: This will precede all paths in <pathMap>, thus avoiding repetition of terms like "obs/" and "var/".'''
        if pfx and not pfx.endswith('/'):
            pfx=pfx+'/';
        df=pd.DataFrame()
        for name, (codePath, catPath) in pathMap.items():
            try:
                data = self.hdf_file[pfx + codePath][:]
                if catPath:  # If categorical path is specified
                    categories = self.hdf_file[pfx + catPath][:]
                    if isinstance(categories[0],bytes):
                        df[name]=pd.Categorical.from_codes(codes=data,categories=[x.decode() for x in categories])
                    else:
                        df[name]=pd.Categorical.from_codes(codes=data,categories=categories)
                elif isinstance(data[0],bytes):
                    df[name] = [x.decode() for x in data]
                else:
                    df[name]=data
            except KeyError as e:
                print(f'Warning: Failed to retrieve data for "{name}" due to missing path: {e}')

        if indexName and indexName in df.columns:
            df.set_index(indexName, inplace=True)
        elif indexName:
            print(f'Warning: Index name "{indexName}" not found in columns. Index not set.')
        
        return df
    
    def getMaskedSparse(self):
        '''Returns the masked data as sparse matrix.'''
        if not hasattr(self,'obs_mask') and not hasattr(self,'var_mask'):
            raise Exception('No mask found for observations or variables. Use setObsMask and/or setVarMask first')
        elif not hasattr(self,'sparse'):
            raise Exception('No sparse data was found. Use setSparse first')
        elif not hasattr(self,'var_mask'):
            try:
                self.var_mask = np.ones(self.Nvar, dtype=bool)
                print('No mask found for variables. Defaulted to all variables.')
            except:
                raise Exception('Variables not found. Use setVar first.')
        elif not hasattr(self,'obs_mask'):
            try:
                self.obs_mask = np.ones(self.Nobs, dtype=bool)
                print('No mask found for observations. Defaulted to all observations.')
            except:
                raise Exception('Observations not found. Use setObs first.')        
        return self.sparse[self.obs_mask, :][:, self.var_mask]
    
    def writeMaskedSparse(self,filepath,renameObs=None,renameVar=None,indexName_obs=None,indexName_var=None):
        '''Writes masked sparse data and masked observations and variables tables (or just observations
        and variables, when tables are not available) to file <filepath>.
        
        Use column renaming maps <renameObs> or <renameVar> to rename columns of obs or var data frames
        that cause conflicts with anndata reserved names. E.g.: renameVar={'_index':'feature'}.
        '''
        import anndata
        data=self.getMaskedSparse()
        if hasattr(self,'obsTable'):
            df_obs=self.obsTable.loc[self.obs_mask].copy()
        else:
            df_obs=pd.DataFrame(index=np.array(self.obs)[self.obs_mask])
        if hasattr(self,'varTable'):
            df_var=self.varTable.loc[self.var_mask].copy()
        else:
            df_var=pd.DataFrame(index=np.array(self.var)[self.var_mask])
        
        if renameObs:
            df_obs.rename(renameObs,axis=1,inplace=True)

        if indexName_obs:
            df_obs.index.name=indexName_obs                        
            
        if renameVar:
            df_var.rename(renameVar,axis=1,inplace=True)

        if indexName_var:
            df_var.index.name=indexName_var  
        
        adata=anndata.AnnData(X=data,obs=df_obs,var=df_var)
        adata.write_h5ad(filepath)
                                           
    def close(self):
        '''Closes the HDF5 file.'''
        self.hdf_file.close()
        
    def __getList(self,path,catPath):
        data=self.hdf_file[path][:]
        if catPath:  # If categorical path is specified
            categories = self.hdf_file[catPath][:]
            data = [categories[i] for i in data]
        if isinstance(data[0],bytes):
            return [x.decode() for x in data]
        else:
            return [x for x in data] 
                
    def __checkSparse(self):
        if hasattr(self,'sparse'):
            if not self.Nvar==self.sparse.shape[1]:
                print('WARNING: The number of rows in the original sparse matrix do not match the the number of variables. Make sure var or varTable are formed from correct paths and then reset the sparse matrix using setSparse() function.')

    def __getitem__(self, ij):
        if not hasattr(self, 'obs_index') or not hasattr(self, 'var_index'):
            raise Exception('You need both obs and var set with the appropriate set functions to use the default getitem method.')
        if not hasattr(self, 'sparse'):
            print('No sparse matrix found to get items from. Using setSparse() function to build one.')
            self.setSparse()  # Ensure this calls the correct method to initialize the sparse matrix

        obs_indices, var_indices = ij

        try:
            # Handle slice objects or single indices for observations
            if isinstance(obs_indices, slice):
                obs_slice = range(len(self.obs))[obs_indices]
            elif isinstance(obs_indices, list):
                obs_slice = [self.obs_index[i] for i in obs_indices]
            else:  # Single index
                obs_slice = self.obs_index[obs_indices]

            # Handle slice objects or single indices for variables
            if isinstance(var_indices, slice):
                var_slice = range(len(self.var))[var_indices]
            elif isinstance(var_indices, list):
                var_slice = [self.var_index[i] for i in var_indices]
            else:  # Single index
                var_slice = self.var_index[var_indices]
        except KeyError:
            # If the name doesn't exist, return None or an empty array
            print('One or more specified observation or variable names do not exist.')
            return np.array([]) 

        # Fetch and return the submatrix
        vals=self.sparse[obs_slice, :][:, var_slice].toarray()
        if vals.shape==(1,1):
            return vals[0][0]
        else:
            return vals
                
    @property
    def Nvar(self):
        if hasattr(self,'var'):
            return len(self.var)
        elif hasattr(self,'varTable'):
            return self.varTable.shape[0]
        elif hasattr(self,'sparse'):
            return self.sparse.shape[1]
        else:
            return None
    
    @property
    def Nobs(self):
        if hasattr(self,'indptr'):
            return len(self.indptr) - 1
        elif hasattr(self,'obs'):
            return len(self.obs)
        elif hasattr(self,'obsTable'):
            return self.obsTable.shape[0]
        else:
            return None
        
    @staticmethod
    def getDataFrameMask(df,cols,vals):
        '''
        Creates and returns a nested mask on DataFrame <df> based on column/value pairings specified in <cols> and <vals>.

        Parameters:
        - <cols>: A list or a single column name specifying the DataFrame columns to be used for masking. The columns 
        are used in the order they are provided for nested masking.
        - <vals>: A list of tuples, where each tuple corresponds to the conditions used for masking the columns 
        specified in <cols>. Each tuple can contain single values or lists of values that are matched against 
        the columns. A single value is used to create a mask where the DataFrame column equals the value. 
        A list of values is used to create a mask where the DataFrame column matches any of the values in the list.

        Examples:
        - Single column, single condition:
          cols = 'organ_tissue'
          vals = ('Skin', ['epithelial cell', 'melanocyte'])
          This will mask rows where 'organ_tissue' is 'Skin' AND 'cell_ontology_class' includes either 
          'epithelial cell' or 'melanocyte'.

        - Multiple columns, multiple conditions (nested AND-OR logic):
          cols = ['organ_tissue', 'cell_ontology_class']
          vals = [('Skin', ['epithelial cell', 'melanocyte']), ('Liver', 'hepatocyte')]
          This creates a mask for rows where 'organ_tissue' is 'Skin' AND 'cell_ontology_class' includes either 
          'epithelial cell' or 'melanocyte', OR where 'organ_tissue' is 'Liver' AND 'cell_ontology_class' is 'hepatocyte'.

          cols = ['organ_tissue', 'cell_ontology_class']
          vals = [(['Skin', 'Thymus'], 'macrophages'), ('Liver', ['hepatocyte', 'endocytes'])]
          This masks rows where 'organ_tissue' includes either 'Skin' or 'Thymus' AND 'cell_ontology_class' 
          is 'macrophages', OR where 'organ_tissue' is 'Liver' AND 'cell_ontology_class' includes either 
          'hepatocyte' or 'endocytes'.

        Operation:
        - Within a tuple in <vals>, conditions are combined using AND logic to create a mask for each column
        specified in <cols>.
        - Between tuples in <vals>, the resulting masks are combined using OR logic. This allows for complex, 
        nested conditional masking of the DataFrame based on multiple criteria.

        The function returns a boolean mask that can be used to filter the DataFrame based on the specified conditions.
        '''
        cols = [cols] if isinstance(cols, str) else cols
        # Wrap vals in a list if it's a single condition not in a list
        if isinstance(vals, tuple) and not isinstance(vals[0], (list, tuple)):
            prevals=[vals]
        elif not isinstance(vals,(list,tuple)):
            prevals=(vals,)
        else:
            prevals=vals
        # Preprocess vals to ensure all elements are properly formatted as tuples
        vals=[]
        for val in prevals:
            if isinstance(val,tuple):
                vals.append(val)
            else:
                vals.append((val,))
        
        or_mask = np.zeros(len(df), dtype=bool)
        for val_tuple in vals:
            and_mask = np.ones(len(df), dtype=bool)
            for col, val in zip(cols, val_tuple):
                # Apply AND mask for each condition within a tuple
                if isinstance(val, (list, tuple)):
                    and_mask &= df[col].isin(val).values
                else:
                    and_mask &= (df[col] == val).values
            # Combine conditions between tuples with OR
            or_mask |= and_mask
        
        return or_mask