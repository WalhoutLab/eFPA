%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new_ID = getNewIndex(IDs)
%
% Benjamín J. Sánchez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newID = getNewIndex(IDs)

%Find latest index and create next one:
IDs   = regexprep(IDs,'[^(\d*)]','');
IDs   = str2double(IDs);
newID = max(IDs) + 1;
newID = num2str(newID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%