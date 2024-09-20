try:
    import pickle5 as pickle
except:
    import pickle

def saveobj(obj,filename):
    fo=open(filename,'wb');
    pickle.dump(obj, fo, pickle.HIGHEST_PROTOCOL)
    fo.close();


def loadobj(filename):
    fo=open(filename, 'rb');
    obj=pickle.load(fo);
    fo.close();
    return obj
    
def list2file(L,fpath):
    '''Writes elements of a list to file at <fpath> row by row.'''
    with open(fpath, 'w') as file:
        for item in L:
            file.write(f"{item}\n")

