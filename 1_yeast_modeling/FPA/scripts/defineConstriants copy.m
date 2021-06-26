function model = defineConstriants(model, infDefault,smallFluxDefault, flux_lb, rxnLabel)
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
%    MFA:               the Metabolic Flux Measurement data; In this
%                       function, the we only use MFA data to determine whether allow a
%                       nutrient or not. Any allowed nutrient is freely avaible
%
% OUTPUT:
%   model:              the constrianed model
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020

if ~any(strcmp(model.rxns,'DM_hom_L_c_'))%not modified model
    % EX_hom_L(e) cannot carry flux ==> add a celluar demand for this
    % it is a co-transporting circular met
    model = addDemandReaction(model,'hom_L[c]');
    % EX_sbt-d(e) can only be uptaken but cannot carry influx ==> change transport reversibility
    model.lb(strcmp(model.rxns,'SBTle')) = -infDefault;
end

% set infinite
model.ub(isinf(model.ub)) = infDefault;
model.lb(isinf(model.lb)) = -infDefault;
% open all exchange to a small flux default
model.lb(cellfun(@(x) ~isempty(regexp(x,'^EX_','once')),model.rxns)) = -smallFluxDefault;
model.lb(cellfun(@(x) ~isempty(regexp(x,'^sink_','once')),model.rxns)) = -smallFluxDefault;
% define the freely avaiable inorganic media content 
media = {'EX_ca2(e)',...
        'EX_cl(e)',...
        'EX_fe2(e)',...
        'EX_fe3(e)',...
        'EX_h(e)',...
        'EX_h2o(e)',...
        'EX_k(e)',...
        'EX_na1(e)',...
        'EX_nh4(e)',...
        'EX_so4(e)',...
        'EX_pi(e)',...
        'EX_o2(e)'};
model.lb(ismember(model.rxns,media)) = -infDefault;% media ion set to free

% define the vitamin input
vitamins = {'EX_btn(e)',...
        'EX_chol(e)',...
        'EX_pnto_R(e)',...
        'EX_fol(e)',...
        'EX_ncam(e)',...
        'EX_pydxn(e)',...
        'EX_ribflv(e)',...
        'EX_thm(e)',...
        'EX_adpcbl(e)',...
        };
model.lb(ismember(model.rxns,vitamins)) = -infDefault;%artificially set as -0.01
% set the maintaince
model = changeRxnBounds(model,'DM_atp_c_','l',0);%1.07 according to palsson
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
            model.lb(strcmp(model.rxns,myrxn)) = -infDefault;
        else % add sink // note applicable
%             if ~any(strcmp(model.rxns,['sink_',myrxn(4:end)]))
%                 model = addSinkReactions(model,regexprep(myrxn(4:end),'_c_$','[c]'),flux_lb(i),1000);
%             else
%                 model.lb(strcmp(model.rxns,['sink_',myrxn(4:end)])) = -infDefault;
%             end
        end
    else %don't allow uptake
        if any(strcmp(model.rxns,myrxn))
            model.lb(strcmp(model.rxns,myrxn)) = 0;
        else % add demand 
            model = addDemandReaction(model,regexprep(myrxn(4:end),'_c_$','[c]'));
        end
    end
end

% additionally, the essential aa histidine and nonessential aa cys has no data
% we give a guess
AA = {'EX_his_L_e_','EX_cys_L_e_'};
model.lb(ismember(model.rxns,AA)) = -infDefault;

% remove parentathsis
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');

end
