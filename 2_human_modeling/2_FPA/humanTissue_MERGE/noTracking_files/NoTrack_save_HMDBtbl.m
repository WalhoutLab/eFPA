tbl = array2table(refMat);
tbl.Properties.VariableNames = measuredTissue;
newName = allCmp_iHumanName;
reps = {'lipoic acid', 'androstenediol', '5-formyl-THF', 'inositol'};
newName(ismember(newName, reps)) = strcat(newName(ismember(newName, reps)),...
                                         '_MetaboAnalystR_name_',...
                                        allCmp_MSEAname(ismember(newName, reps)));
tbl.Properties.RowNames = newName;
writetable(tbl,'output/HMDB_metaboliteset_processed.csv','WriteRowNames',true);