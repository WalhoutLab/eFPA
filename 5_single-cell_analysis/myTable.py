import pandas
import numpy as np

class pdtable(pandas.core.frame.DataFrame):

    def printRow(self,rowname,rowcol=None):
        '''Prints the attributes of a row indicated by <rowname> either in the index (default) or in <rowcol>.'''
        if rowcol:
            x=self.loc[(self[rowcol]==rowname)];
        else:
            x=self.loc[[rowname]];
        multirow=len(x.index)>1;
        for row in x.index:
            if multirow: print('Row:\t'+str(row));
            for col in x.columns:
                print(col+'\t'+str(x.loc[row,col]));
            if multirow: print('~~end of row~~\n');     
        
    def insertColumn(self,colname,values,refcol,position='before'):
        '''Inserts a column named <colname> with values <values> before or after
        <refcol> as indicated by <position>.'''
        if not position.lower() in ['after','before']:
            raise Exception('<position> must be "before" or "after". '+position+' found.')
        i=self.columns.get_loc(refcol);
        if position.lower()=='before':
            i=i;
        else:
            i=i+1;
        self.insert(i,colname,values);   

    ##############SEARCH##############
    def searchcols(self,regexp,cols,exact=False,printscreen=True,onlysearched=True):
        '''Searches <regexp> in columns <cols>. Returns a dataframe of rows where the searched value were found.
        Arguments
            regexp: Regular expression.
            cols: A string that indicates a single column name or a list of column names. In the case of a list
                  ":" operator is allowed.
            exact: If true, <regexp> is searched as a value.
            printscreen: If true, the dataframe of rows with searched values is printed.
            onlysearched: If true, only the searched columns are included in the resulting dataframe.
            '''
        if not type(cols)==str:
            I=self.__resolveVars(cols);
            cols=[self.columns[i] for i in I];
        else:
            cols=[cols];
        Sindex=set();
        if not exact:
            for col in cols: 
                df=self.__searchsinglecol(regexp,col);
                Sindex=Sindex.union(df.index);
        else:
            for col in cols:
                df=self.__searchsinglecol_exact(regexp,col);
                Sindex=Sindex.union(df.index);
        if len(cols)>1:
            Lindex=list(Sindex);
            Lindex.sort();
            df=self.loc[Lindex,];
        if onlysearched:
            df=df[cols];
        if printscreen:
            print(df);
        return pdtable(df)


    def __searchsinglecol(self,regexp,col):
        '''Subroutine of searchcols.'''
        return self[self[col].str.contains(regexp,na=False)];

    def __searchsinglecol_exact(self,regexp,col):
        '''Subroutine of searchcols.'''
        Io=np.where(self[col].values==regexp)[0];
        I=[self.index[i] for i in Io];
        return self.loc[I,];  
        
    #############PRIVATE################                   
    def __resolveVars(self,vars,returnNames=False):
        '''Converts variables as column names in the string list <vars> to a list of column indices. 
        ":" operator is allowed. If <returnNames>, indices are converted to a list of colnames before
        returning.'''
        I=[];
        if not type(vars)==list:
            vars=[vars];
        for x in vars:
            L=re.split(':',x);
            if len(L)==2:
                I=I+list(range(self.columns.get_loc(L[0]),self.columns.get_loc(L[1])+1));
            elif len(L)==1:
                I=I+[self.columns.get_loc(L[0])];
            else:
                raise Exception('At least one string in '+str(vars)+' seems to have either a null or an element with multiple ":" symbols.');    
        if returnNames:
            return [self.columns[i] for i in I]
        else:
            return I        

###################################################################################
############################## MATRIX #############################################
###################################################################################

class matrix:
    
    def __init__(self,vals,rows=[],cols=[]):
        '''Forms a numpy mattrix that can be edited and accessed based on row names and column  names.
        <vals>: 2D list or numpy array of values, or a tuple of form (rowno,columnno,val). If tuple is provided,
                an array of the single value val will be generated with the dimensions indicated (rowno and columnno).
        <rows>: Row names. If not provided, integers of row order will be used.
        <cols>: Column names. If not provided, integers of column order will be used.
        '''
        #Form the value array
        if type(vals)==list:
            self.array=np.array(vals);
        elif type(vals)==tuple:
            if vals[2]==0:
                self.array=np.zeros((vals[0],vals[1]));
            else:
                self.array=np.ones((vals[0],vals[1]));
                self.array=self.array*vals[2];
        else:
            self.array=vals.copy();
        self.rowno,self.colno=self.array.shape;
        #Make row and column dictionaries        
        if rows:
            self.rows=[x for x in rows];
        else:
            self.rows=[x for x in range(self.rowno)];
        self.drows={self.rows[i]:i for i in range(self.rowno)};
        if not self.rowno==len(self.rows):
            raise Exception('Provided rows and array shape do not agree.')
        if cols:
            self.cols=[x for x in cols];
        else:
            self.cols=[x for x in range(self.colno)];
        self.dcols={self.cols[i]:i for i in range(self.colno)}; 
        if not self.colno==len(self.cols):
            raise Exception('Provided columns and array shape do not agree.')        


    def getrow(self,row,nonZero=False):
        '''Returns <row> as array.'''
        return self.array[self.drows[row],:];
    
    def __getitem__(self,ij):
        try:
            return self.array[self.drows[ij[0]],self.dcols[ij[1]]] 
        except:
            print('item not found')
            return None    

    def __setitem__(self,ij,val):
        try:
            self.array[self.drows[ij[0]],self.dcols[ij[1]]]=val; 
        except:
            raise Exception('item not found')

def csv2matrix(filename):
    '''Converts csv table with rows as first column and columns as first row to
    a myTable.matrix incidence and returns it.'''
    df=pandas.read_csv(filename,header=0,index_col=0);
    return matrix(df.to_numpy(),rows=df.index.to_list(),cols=df.columns.to_list())