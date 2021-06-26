function model = defineConstriants_min(model, infDefault,smallFluxDefault, flux_lb, rxnLabel)
% This function is to define the uptake constrainst for a native human
% model RECON2.2. It is not designed or tested for any other model.
%
% USAGE:
%
%    model = defineConstriants(model, infDefault,smallFluxDefault, FVA)
%
% INPUTS:
%    model:             input RECON2.2 model (COBRA model structure)
%    infDefault:        the default value for infinite fluxes
%    smallFluxDefault:  the default value for trace uptake fluxes
%    MFA:               the Metabolic Flux Measurement data
%
% OUTPUT:
%   model:              the constrianed model
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020

exRxns = findExcRxns(model);
% open all exchange to a small flux default
model.lb(exRxns) = -smallFluxDefault;
% define the freely avaiable inorganic media content 
media = {'HMR_9082',...
        'HMR_9150',...
        'HMR_9076',...
        'HMR_9096',...
        'HMR_9079',...
        'HMR_9047',...
        'HMR_9081',...
        'HMR_9077',...
        'EX_nh4[e]',...
        'HMR_9074',...
        'HMR_9072',...
        'HMR_9048'%,...
        %'HMR_9269',...%start here are added because of essential 
        %'HMR_9167',...
        %'EX_retinal[e]'...
        %'EX_HC00004[e]',...
        %'12dgr120_x'
        };
model.lb(ismember(model.rxns,media)) = -infDefault;% media ion set to free

% define the vitamin input
vitamins = {'HMR_9109',...
        'HMR_9083',...
        'HMR_9145',...
        'HMR_9146',...
        'HMR_9378',...
        'HMR_9144',...
        'HMR_9143',...
        'HMR_9159',...
        'EX_adpcbl[e]',...
        };
model.lb(ismember(model.rxns,vitamins)) = -0.01;%artificially set as -0.01
% set the maintaince
model = changeRxnBounds(model,'HMR_3964','l',1.07);%1.07 according to palsson
% the flux in the analysis is reported as mmol/gdw/h

% remove parentathsis
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% we use the minimal uptake rate for a metabolite when multiple conditions are provided (this will be used to calculate epsilon)
% and the constriant will be the real uptake rate when only one condition
% is provided
for i = 1:size(rxnLabel,1)
    myrxn = rxnLabel{i};
    if flux_lb(i) < 0 %could be uptaken
        if any(strcmp(model.rxns,myrxn))
            model.lb(strcmp(model.rxns,myrxn)) = flux_lb(i);
        else % add sink
            if ~any(strcmp(model.rxns,['sink_',myrxn(4:end)]))
                model = addSinkReactions(model,myrxn(5:end),flux_lb(i),1000);
            else
                model.lb(strcmp(model.rxns,['sink_',myrxn(5:end)])) = flux_lb(i);
            end
        end
    else %don't allow uptake
        if any(strcmp(model.rxns,myrxn))
            model.lb(strcmp(model.rxns,myrxn)) = 0;
        else % add demand 
            model = addDemandReaction(model,myrxn(5:end));
        end
    end
end

% additionally, the essential aa histidine and nonessential aa cys has no data
% we give a guess
AA = {'HMR_9038','HMR_9065'};
model.lb(ismember(model.rxns,AA)) = -0.01;

% remove parentathsis
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');

end
