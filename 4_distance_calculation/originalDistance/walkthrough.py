from MetabolicDistance import *

##1-Loading model variables.

#You need stoichiometry matrix (S, a numpy object) and lists of reactions, metabolites,
#upper and lower boundaries for reactions, and byproducts (hub metabolites
#such as atp and proton). Metabolites in the byproducts list will be ignored during distance
#analysis. 
#You can use loadVarsFromFilenames() function to load all variables in a particular order,
#provided that S is in npy format or delimited text format, and all others are as text with 
#each list element occupying a line. See help for this function for details. 

#EXAMPLE
#If S matrix is in comma-delimited text format:
S,rxns,metabolites,lb,ub,byproducts=loadVarsFromFilenames('Input/Smatrix_regular.txt',\
'Input/reactions_regular.txt','Input/metabolites_regular.txt','Input/LB_regular.txt',\
'Input/UB_regular.txt','Input/byproducts_regular.txt',delimiterInSmatrix=',');
#If S matrix is in an npy file, use the following instead:
#S,rxns,metabolites,lb,ub,byproducts=loadVarsFromFilenames('Input/Smatrix_regular.npy',\
#'Input/reactions_regular.txt','Input/metabolites_regular.txt','Input/LB_regular.txt',\
#'Input/UB_regular.txt','Input/byproducts_regular.txt',delimiterInSmatrix=None);
#NOTE: Parentheses ("(","[","]", and ")") and other special characters are not allowed in reaction names.

##2-Forming a reaction network with unidirectional reactions

#The loaded variables have all information needed. Then you can use the establishRxnNetwork() 
#function to create an instance of rxnnetwork class, which can be used to traverse the
#reaction network. The algorithm first converts all reactions to "forward" and "reverse"
#reactions if any of these directions is available based on the upper and lower boundary constraints. 
#These two directions are indicated by "f" and "r" in the end of reaction ID
#(original ID is in the rxns list above), respectively.

#EXAMPLE
Net=establishRxnNetwork(S,rxns,metabolites,lb,ub,byproducts);
#Net is the object you can use for traversing the reaction network.  


##3-Finding Distances

#To find shortest distances from a reaction of interest (ROI) to all reachable reactions in the 
#network, you can use findDistances() function. The result is a dictionary that maps the ID of each
#reaction to its distance from the ROI. Note that the reaction ID for the ROI includes the 
#directionality letter ("f" or "r") in the end. Likewise all reaction IDs in the distance
#dictionary has the directionality indicator in the end.

#NOTE: If a reaction does not exist in the distance dictionary, then it is not reachable (distance
#may be taken as inf or an arbitrary large number).

#EXAMPLE (ROI: RM04432f)
Ddistance=findDistances('RM04432f',Net);
#To get the distance of a particular reaction from ROI:
print Ddistance['RM01608f']; 
#The answer is 3


##4-More detailed analyses

#If you want to study the shortest paths found from ROI to a reachable reaction, use findPaths()
#function instead. The result is a dictionary that maps reaction ID to path found as a list of sequential
#reactions. This function also outputs a list of reactions for which loops were encountered 
#and fixed, and a list of reactions for which errors were encountered, in that order. 

#NOTE: The shortest path is not a unique solution. There may be multiple paths of the shortest length, the program
#outputs only one.

#EXAMPLE:
Dpath,rxns_wloop,errors=findPaths('RM04432f',Net);
#To get the shortest path to a particular reaction from ROI:
print Dpath['RM01608f']; 
#The answer is ['RM04432f', 'RM03045f', 'RM03158f', 'RM01608f']
#The path from RM04432f to RM01608f is as follows:
#RM04432f: [m] : etfox + ppcoa <==> etfrd + prpncoa
#RM03045f: [m] : h2o + prpncoa <==> 3hpcoa
#RM03158f: [m] : 3hpcoa + h2o --> 3hpp + coa + h
#RM01608f: [m] : 3hpp + nad <==> h + msa + nadh

#Dpath dictionary can be converted to a distance dictionary using convertPaths2Distances()
Ddistance2=convertPaths2Distances(Dpath);

