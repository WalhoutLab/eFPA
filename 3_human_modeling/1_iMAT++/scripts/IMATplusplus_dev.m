function [OFD, PFD, N_highFit, N_zeroFit, minLow, minTotal, minTotal_OFD, MILP, MILP_PFD, HGenes, RLNames, OpenGene, latentRxn, Nfit_latent, wasteDW] = IMATplusplus_dev(model,epsilon_f,epsilon_r,ExpCateg, modelType, speedMode,minLowTol,doMinPFD, doLatent,latentCAP,sideCAP, storeProp,SideProp,ATPm,verbose)
% Uses the iMAT++ algorithm (`Yilmaz et al., 2020`) to find the optimal flux distribution that best fit categorized gene expression data. 
%
% iMAT++ algorithm performs multi-step fitting to find a flux distribution
% that best agrees with rarely, lowly and highly expressed genes in the
% expression profile. It employs a gene-centric feature that minimizes
% false negative fitting and conflicts due to ambiguous GPR assignment. 
%
% USAGE:
%
%    OFD = IMATplusplus(model,doLatent,storeProp,SideProp,epsilon_f,epsilon_r, ATPm, ExpCatag,doMinPFD,latentCAP,modelType)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    epsilon_f:         the epsilon sequence for the forward direction of all reactions in the model. Non-applicable reactions (i.e., for a irreversible,
%                       reverse-only reaction) should have a non-zero default to avoid numeric error.
%    epsilon_r:         the epsilon sequence for the reverse direction of all reactions in the model. Non-applicable reaction (i.e., for a irreversible,
%                       forward-only reaction) should have a non-zero default to avoid numeric error.
%    ExpCateg:          expression categories used for iMAT++ fitting. 
%                       The `ExpCateg` should be a structure variable with "high", "low", "dynamic" and
%                       "zero" four fields. See the walkthrough scripts (and our Github) on how to generate
%                       them from raw expression quantification data (i.e., TPM)
% OPTIONAL INPUTS:
%    modelType:         an integer to specify the input COBRA model. 
%                       1 == dual C. elegans tissue model
%                       2 == generic C. elegans model
%                       3 == other COBRA models. 
%                       iMAT++ has some custom actions (such as calculating
%                       bacteria waste) that are specific to C. elegans
%                       models. So, user should specify the correct model
%                       type to avoid errors from the custom actions.
%    speedMode:         (1, 2, or 3) to indicate which speed mode is used.
%                       We offer three modes of iMAT++ for the different 
%                       balance between computational speed and optimization
%                       stringency. The speed mode level 1 is the original
%                       iMAT++ configuration used in our tissue modeling.
%                       The level 3 gives fastest speed, so is recommanded 
%                       for all complex models, such as RECON models. 
%                       In level 3, the low reactions (dependent on rarely 
%                       and lowly expressed genes) will be constrained by 
%                       rigid boundaries after flux minimization. Additionally, 
%                       we release the MILP stringency to gain computational 
%                       speed. In level 2, we only release MILP strigency.
%                       But in general, the three modes give similar flux
%                       predictions.
%    minLowTol:         the numerical tolerance for total flux of the 
%                       reactions dependent on lowly and rarely expressed 
%                       genes. This parameter tunes the strigency of 
%                       constraining low reaction fluxes to the minimal level.
%                       Default value (1e-5) provides very strigent constraint
%                       that pushes the total flux to the minimal level.
%                       However, if the lowly expressed and rarely expressed
%                       genes extensively conflict with highly expressed
%                       genes, we recommend to use a larger tolerance such 
%                       as the default epsilon (ideally allowing one
%                       mis-fitting).
%    doMinPFD:          whether to perform flux minimization of PFD 
%                       (after minimizing total flux for lowly and rarely expressed genes). 
%                       PFD flux minimization is required to obtain OFD, but  
%                       could be omitted if one only wants to get MILP fitting 
%                       result.
%    doLatent:          whether to perform recursive fitting of latent reactions
%    latentCAP:         the total flux cap for recursive fitting of latent
%                       reactions. The total flux will be capped at (1 +
%                       latentCAP)*OriginalTotalFlux; The default cap is
%                       0.01 (1%)
%    storeProp:         (dual-C. elegans tissue network only, use empty ('[]') for other models) the maximum allowed storage molecule uptake, as w/w percentage of bacteria uptake                  
%    SideProp:          (dual-C. elegans tissue network only, use empty ('[]') for other models) the maximum allowed side metabolites uptake, as w/w percentage of bacteria uptake                  
%    ATPm:              (C. elegans network only, use empty ('[]') for other models) the Non Growth Associated Maintenance (NGAM) used in fitting               
%    verbose:           (0 or 1) to show the MILP log or not
%
%
% OUTPUT:
%   OFD:                the Optimal Flux Distribution (OFD)
% OPTIONAL OUTPUTS:
%   PFD:                the primary flux distribution (PFD)
%   N_highFit:          the number of highly expressed genes fitted
%   N_zeroFit:          the number of rarely expressed genes fitted
%   minLow:             the minimal total flux of reactions dependent on
%                       rarely and lowly expressed genes
%   minTotal:           the minimal total flux of PFD (primary flux distribution)
%   minTotal_OFD:       the minimal total flux of OFD
%   MILP:               the MILP problem (in COBRA format) in the final
%                       flux minimization of OFD. This MILP can serve as a
%                       startpoint for any further analysis such FVA.
%   MILP_PFD:           the MILP problem (in COBRA format) right before the
%                       first Flux minimization. This MILP is readily
%                       constrained for high/zero genes' fitting , and flux
%                       minimization of low reactions. This includes all
%                       data-driven constraints, so it is a conservative 
%                       but good startpoint for further custom analysis.
%   HGenes:             the list of highly expressed genes to be fitted
%   RLNames:            the list of reactions dependent on rarely expressed genes
%   OpenGene:           the list of fitted (carrying flux) highly expressed genes
%   latentRxn:          the list of all identified latent reactions to be fitted 
%   Nfit_latent:        the (total) number of latent reactions fitted
%   wasteDW:            (C. elegans network only) the percentage of wasted bacterial biomass (DW/DW)

% Citation:
% Yilmaz, L. S., Li, X., Nanda, S., Fox, B., Schroeder, F., & Walhout, A. J. (2020). Modeling tissue‐relevant Caenorhabditis elegans metabolism at network, pathway, reaction, and metabolite levels. Molecular Systems Biology, 16(10).
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 2020
%% this is the main integration function
% apply default parameters
if (nargin < 5) % type of the input model 
    modelType = 1; % 1==>dual C.elegans; 2==>generic C. elegans; 3==>non-c.elegans
end
if (nargin < 6) 
    speedMode = 1; % by default, do normal IMAT++ (speed mode = 1)
end
if (nargin < 7)
    if speedMode ~= 3
        minLowTol = 1e-5; % default value of low reaction flux tolerance
    else
        minLowTol = 0; % not applicable to speed mode 3
    end
end
if (nargin < 8)
    doMinPFD = true;
end
if (nargin < 9)
    doLatent = true;
end
if (nargin < 10) % the flux tolerance cap of total flux in the latent fitting. By default, we use 5%
    latentCAP = 0.05;
end
if (nargin < 11) % for ihuman model only
    if modelType == 3
        sideCAP = 1.1;
    else
        sideCAP = [];
    end
end
if (nargin < 12) % for C. elegans' dual-tissue model only
    if modelType == 1
        storeProp = 0.01;
    else
        storeProp = [];
    end
end
if (nargin < 13) % for C. elegans' dual-tissue model only
    if modelType == 1
        SideProp = 0.02;
    else
        SideProp = [];
    end
end
if (nargin < 14) % for C. elegans' model only
    if modelType ~= 3
        ATPm = 10;
    else
        ATPm = [];
    end
end
if (nargin < 15) 
    verbose = 0; % by default, don't show the MILP details
end

% set global constant 
bacMW=966.28583751; % only will be used for C. elegans model
if speedMode > 1
    relMipGapTol = 0.001; % we release the MIPgap to 0.1% 
    % this released MipGap only applies to latent flux fitting step.
    % The optimization stringency of PFD is still kept.
    % users can release the MipGap for PFD manually if needed
else
    relMipGapTol = 1e-12;
end
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for the current solver! Please use Gurobi for best performance!\n')
end
%% mapping the gene categories to reactions
fprintf('Start flux fitting... \n');
tic()
fprintf('Processing expression data... \n');
worm = model;
% process gene expression data
expression_gene=struct;
expression_gene.value = [3*ones(length(ExpCateg.high),1);2*ones(length(ExpCateg.dynamic),1);1*ones(length(ExpCateg.low),1);zeros(length(ExpCateg.zero),1)];
expression_gene.gene = [ExpCateg.high;ExpCateg.dynamic;ExpCateg.low;ExpCateg.zero];
[expressionRxns] = mapExpressionToReactions_xl(worm, expression_gene);
RHNames = worm.rxns(expressionRxns == 3);
RLNames = worm.rxns(expressionRxns == 0); % for the first step integration
% only genes associated with high reactions are high genes! (not the input
% high category genes!)

HGenes = {};
for i = 1: length(ExpCateg.high)
    myrxns = worm.rxns(any(worm.rxnGeneMat(:,strcmp(ExpCateg.high(i),worm.genes)),2),1);
    if any(ismember(myrxns,RHNames))
        HGenes = [HGenes;ExpCateg.high(i)];
    end
end

%% begin to do iMAT++ integration
%% Step1: do MILP integration of highly expressed genes and rarely expressed genes
toc()
fprintf('fitting high genes with MILP... \n');
tic()
% setup model specifc constraints for C. elegans mdoel
if modelType == 1
    if storeProp>=0 && SideProp>=0 && ATPm>=0
        % adjust the ATPm according to the bacterial uptake. This requires 
        % manual tuning
        worm = changeRxnBounds(worm,'RCC0005_I',ATPm,'l');
        worm = changeRxnBounds(worm,'RCC0005_X',ATPm,'l');
        % change the side proportion
        worm.S(end-1, strcmp('EXC0050_L',worm.rxns)) = storeProp*bacMW*0.01;
        worm.S(end, strcmp('EXC0050_L',worm.rxns)) = SideProp*bacMW*0.01;
    else
        error('Please check your input parameters for side proportion, storage proportion, and ATPm!');
    end
elseif modelType == 2
    worm = changeRxnBounds(worm,'RCC0005',ATPm,'l');
elseif modelType == 3
    fprintf('Doing integration for ihuman model, no ATPm imposed... \n');
    % worm.S(end, strcmp('EX_majorNutr',worm.rxns)) = SideProp;
end
% perform the gene-centric MILP fitting of highly expressed genes and reactions dependent on rarely expressed genes
[~, ~,~,solution, MILProblem] = ExpressionIntegrationByMILP(worm, RHNames, RLNames,HGenes, epsilon_f, epsilon_r,[],[],verbose);
toc()
% convert the objective value in a solution to the model constraints
MILProblem =  solution2constraint(MILProblem,solution);
NfitInd = length(MILProblem.b);%save the index of this constraint for later use

%% Step2: do minimization of low reactions
fprintf('Minimizing low flux... \n');
tic()
% add new variables (absolute flux proxy) for flux minimization
MILProblem = addAbsFluxVariables(MILProblem, worm);

% set the initial solution (as the last step solution)
% only recommend to use when numerical issue is severe
% MILProblem.x0 = [solution.full;abs(solution.full(1:length(worm.rxns)))];

% minimize low flux (reactions dependent on lowly and rarely expressed
% genes)
% please note that, to have minimal number of variables in the MILP problem, 
% we do the flux minimization by the absolute flux proxy variables we just
% generated for total flux minimization

% create a new objective function
% Creating c vector (objective function)
c_v = zeros(size(worm.S,2),1);
c_yh1 = zeros(length(RHNames),1);
c_yl = zeros(length(RLNames),1);
c_yh2 = zeros(length(RHNames),1);
c_minFluxLowRxns = zeros(length(worm.rxns),1);
c_minFluxLowRxns(expressionRxns == 1 | expressionRxns == 0) = 1; % set the coeffi of "low" and "rare" reactions as 1
c_gene = zeros(length(HGenes),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
MILProblem.c = c;
% set the objective sense as minimize 
MILProblem.osense = 1;
% sovle the MILP problem
solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
% we use customized solver interface to fine-tune the solver parameters for gurobi.
% However, user is able to replace `solveCobraMILP_XL` with original
% `solveCobraMILP` function in COBRAtoolbox. If using the original
% function, please specify the optimal, feasible and integer tolerance!

if solution.stat == 0 && solverOK % when failed to solve (infeasible), we start to tune solver parameter 
    % NOTE: SPECIFIC TO GUROBI SOLVER!
    fprintf('initial solving failed! Start the auto-tune...#1\n')
    gurobiParameters = struct();
    gurobiParameters.Presolve = 0;
    solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
    if solution.stat == 0
        fprintf('initial solving failed! Start the auto-tune...#2\n')
        gurobiParameters.NumericFocus = 3;
        solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0
            error('MILP solving failed! Please inspect the reason!');
        end
    end
end        
if solution.stat ~= 1
    if solution.stat == 3
        % pass with warning
        warning('Low fluxes are not fully eliminated! (solver time out! Change timeLimit if needed!)');
    else % unknown error
        error('MILP solving failed! Please inspect the reason!');
    end
end
toc()
fprintf('Minimizing low flux completed! \n');
minLow = solution.obj;
%% step3: minimize total flux as needed
if speedMode == 3
    numTol = 1e-8; % too small tolerance causes numerical instability; too large lose the flux minimization's effect
    % constrain by current flux plus a numeric tolerance 
    MILProblem.ub(expressionRxns == 1 | expressionRxns == 0) = abs(solution.full(expressionRxns == 1 | expressionRxns == 0)) + numTol;
    MILProblem.lb(expressionRxns == 1 | expressionRxns == 0) = -abs(solution.full(expressionRxns == 1 | expressionRxns == 0)) - numTol;   
    % all absolute zero fluxes stay zero
    MILProblem.ub((expressionRxns == 1 | expressionRxns == 0) & (solution.full(1:length(worm.rxns)) == 0)) = 0;
    MILProblem.lb((expressionRxns == 1 | expressionRxns == 0) & (solution.full(1:length(worm.rxns)) == 0)) = 0;
    
    % then remove redundant variables (the zero reaction integer variables)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = ones(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = zeros(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    lowInd = boolean([c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns]);
    MILProblem.A(:,lowInd) = [];% remove variables
    NfitZero = sum(abs(solution.full(expressionRxns == 0)) < 1e-8); % default int tol of iMAT++
    MILProblem.b(NfitInd) = MILProblem.b(NfitInd) - NfitZero; % update the Nfit constraints
    MILProblem.lb(lowInd) = [];% update boundary
    MILProblem.ub(lowInd) = [];% update boundary 
    MILProblem.vartype(lowInd) = [];% update variable type
    RLNames_ori = RLNames;
    RLNames = [];% update RLNames (should be empty because variables are deleted)
else
    % set a minFlux tolerance 
    tol = minLowTol; % it could be solver tolarence or some larger number to increase the numeric stability and allow flexible fitting of zero and low reactions
    solution.obj = solution.obj + tol;
    MILProblem = solution2constraint(MILProblem,solution);
end
% MILProblem.x0 = solution.full;
MILP_PFD = MILProblem;
MILP = MILProblem; % write the MILP output. If latent is done, this MILP will be updated.
if doMinPFD
    %% 12052020 add major nutrient uptake constriants - for human integration
    if modelType == 3 && sideCAP ~= -1 % run the dynamic constriant for unlimited bacUp
        MILProblemTemp = MILProblem;
        MILProblemTemp.c = zeros(size(MILProblem.A,2),1);
        MILProblemTemp.c(find(strcmp(worm.rxns,'EX_sideNutr'))) = 1;
        MILProblemTemp.osense = 1;         
        solution = solveCobraMILP_XL(MILProblemTemp, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
            fprintf('initial solving failed! Start the auto-tune...#1\n')
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            solution = solveCobraMILP_XL(MILProblemTemp,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
            if solution.stat == 0
                fprintf('initial solving failed! Start the auto-tune...#2\n')
                gurobiParameters.NumericFocus = 3;
                solution = solveCobraMILP_XL(MILProblemTemp,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
                if solution.stat == 0
                    error('MILP solving failed! Please inspect the reason!');
                end
            end
        end  
        if solution.stat ~= 1
            if solution.stat == 3
                % pass with warning
                warning('the minimal side nutrient demand is not fully minimized! (solver time out! Change timeLimit if needed!)');
            else % unknown error
                error('MILP solving failed! Please inspect the reason!');
            end
        end
        fprintf('the minimal side nutrient demand is %f \n',solution.obj);
        % impose boundary constraint
        MILProblem.ub(find(strcmp(worm.rxns,'EX_sideNutr'))) = sideCAP * solution.obj;
        worm.ub(find(strcmp(worm.rxns,'EX_sideNutr'))) = sideCAP * solution.obj;
    end
    
    % minimize total flux
    % again, to reduce computational burden, we use the same proxy
    % variables used in the low flux minimization
    % create a new objective function
    % Creating c vector (objective function)
    c_v = zeros(size(worm.S,2),1);
    c_yh1 = zeros(length(RHNames),1);
    c_yl = zeros(length(RLNames),1);
    c_yh2 = zeros(length(RHNames),1);
    c_minFluxLowRxns = ones(length(worm.rxns),1);
    c_gene = zeros(length(HGenes),1);
    c = [c_v;c_yh1;c_yl;c_yh2;c_gene;c_minFluxLowRxns];
    MILProblem.c = c;
    fprintf('Minimizing total flux... \n');
    tic()
    solution = solveCobraMILP_XL(MILProblem, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
    if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
        fprintf('initial solving failed! Start the auto-tune...#1\n')
        gurobiParameters = struct();
        gurobiParameters.Presolve = 0;
        solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0
            fprintf('initial solving failed! Start the auto-tune...#2\n')
            gurobiParameters.NumericFocus = 3;
            solution = solveCobraMILP_XL(MILProblem,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
            if solution.stat == 0
                error('MILP solving failed! Please inspect the reason!');
            end
        end
    end  
    if solution.stat ~= 1
        if solution.stat == 3
            % pass with warning
            warning('Total flux is not fully minimized! (solver time out! Change timeLimit if needed!)');
        else % unknown error
            error('MILP solving failed! Please inspect the reason!');
        end
    end
    
    if modelType == 3 && sideCAP ~= -1 % run the dynamic constriant for unlimited bacUp
        % maximazing NADPH production
        solution.obj = solution.obj * 1.1;
        MILProblem = solution2constraint(MILProblem,solution);
        
        totalNADPH = cellfun(@(x) ~isempty(regexp(x,['^NADPH \[(\w|\s)*\]$'],'once')),model.metNames);
        S_nadph_prod = model.S(totalNADPH,:);
        S_nadph_prod(S_nadph_prod<0) = 0;
        S_nadph_prod = sum(S_nadph_prod);
        MILProblemTemp = MILProblem;
        MILProblemTemp.c = zeros(size(MILProblem.A,2),1);
        MILProblemTemp.c(1:length(S_nadph_prod)) = S_nadph_prod;
        MILProblemTemp.osense = -1;         
        solution = solveCobraMILP_XL(MILProblemTemp, 'timeLimit', 300, 'logFile', 'MILPlog', 'printLevel', verbose);
        if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter #NOTE: SPECIFIC TO GUROBI SOLVER!%
            fprintf('initial solving failed! Start the auto-tune...#1\n')
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            solution = solveCobraMILP_XL(MILProblemTemp,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
            if solution.stat == 0
                fprintf('initial solving failed! Start the auto-tune...#2\n')
                gurobiParameters.NumericFocus = 3;
                solution = solveCobraMILP_XL(MILProblemTemp,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose);
                if solution.stat == 0
                    error('MILP solving failed! Please inspect the reason!');
                end
            end
        end  
        if solution.stat ~= 1
            if solution.stat == 3
                % pass with warning
                warning('the minimal side nutrient demand is not fully minimized! (solver time out! Change timeLimit if needed!)');
            else % unknown error
                error('MILP solving failed! Please inspect the reason!');
            end
        end
        fprintf('the maximum NADPH flux is %f \n',solution.obj);
       
    end
    
    minTotal = solution.obj;
%     solution.obj = solution.obj*(1+latentCAP); % 5% minFlux constraint!
%     MILProblem2 = solution2constraint(MILProblem,solution);
%     MILProblem2.x0 = solution.full;  
    %% 12052020 add major nutrient uptake constriants - for human integration
    % when flux minimization is done, remove the nutrient constraint
%     if modelType == 3 && sideCAP ~= -1 % run the dynamic constriant for unlimited bacUp
%         % impose boundary constraint
%         MILProblem2.ub(find(strcmp(worm.rxns,'EX_sideNutr'))) = 1000;
%         worm.ub(find(strcmp(worm.rxns,'EX_sideNutr'))) = 1000;
%     end
%     MILP2 = MILProblem2; 
    % this is a special error-handling section. In rare cases, the highly
    % expressed genes are strongly conflicting with rarely expressed genes.
    % so, our minimization of low reaction fluxes will essentially kills
    % the fitting of most if not all high genes. Finally the minimal total flux
    % could be a very small number, even zero. We raise an error for such 
    % invalid fitting.
    if solution.obj <= minLowTol % minLowTol represents a small flux tolerance 
        error('All the highly expressed genes conflict with at least one rarely expressed gene. The flux fitting is not valid!')
    end
else
    minTotal = 0;
end
fprintf('Primary Flux Fitting completed! \n');
toc()
% make output of primary flux distribution and some optional outputs
FluxDistribution = solution.full(1:length(worm.rxns));
PFD = solution.full(1:length(worm.rxns));
if speedMode == 3 % the low fitting is not controlled by binary variable
    nameL = worm.rxns(ismember(worm.rxns,RLNames_ori));
    fluxL = solution.full(ismember(worm.rxns,RLNames_ori));
    ClosedLReaction = nameL(abs(fluxL)<1e-9);
else
    nameL = worm.rxns(ismember(worm.rxns,RLNames));
    yL = boolean(solution.int((length(RHNames)+1):(length(RHNames)+length(RLNames))));
    ClosedLReaction = nameL(yL);
end
OpenGene = HGenes(boolean(solution.int((2*length(RHNames)+length(RLNames))+1:end)));
N_highFit = length(OpenGene);
N_zeroFit = length(ClosedLReaction);
if doLatent
    %% step4. make the latent rxns fitting
    [FluxDistribution,latentRxn,Nfit_latent,minTotal_OFD,MILP] = fitLatentFluxes(MILP2, worm,PFD, HGenes,epsilon_f,epsilon_r,latentCAP,verbose,relMipGapTol,sideCAP);
    OFD = FluxDistribution;
else
    OFD = [];
    latentRxn = [];
    Nfit_latent = [];
    minTotal_OFD = [];
end
%% provide some major evaluation of the flux distribution (only for C. elegans model)
if modelType == 1
    fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050_L'))); % bac
    bacWaste = ...
        FluxDistribution(strcmp(worm.rxns,{'EXC0051_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0052_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0053_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0054_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0055_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0056_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0057_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0058_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0059_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0060_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0063_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0065_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0072_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC_L[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0142_L'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC_L[e]')});
    wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050_L')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC_L[e]')}));
    fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
elseif modelType == 2
    fprintf('the bacterial uptake is: %f\n',FluxDistribution(strcmp(worm.rxns,'EXC0050'))); % bac
    bacWaste = ...
        FluxDistribution(strcmp(worm.rxns,{'EXC0051'})) * MolMass(model.metFormulas{strcmp(model.mets,'protein_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0052'})) * MolMass(model.metFormulas{strcmp(model.mets,'rna_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0053'})) * MolMass(model.metFormulas{strcmp(model.mets,'dna_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0054'})) * MolMass(model.metFormulas{strcmp(model.mets,'peptido_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0055'})) * MolMass(model.metFormulas{strcmp(model.mets,'lps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0056'})) * MolMass(model.metFormulas{strcmp(model.mets,'lipa_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0057'})) * MolMass(model.metFormulas{strcmp(model.mets,'outerlps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0058'})) * MolMass(model.metFormulas{strcmp(model.mets,'pe_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0059'})) * MolMass(model.metFormulas{strcmp(model.mets,'pg_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0060'})) * MolMass(model.metFormulas{strcmp(model.mets,'ps_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0063'})) * MolMass(model.metFormulas{strcmp(model.mets,'clpn_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0065'})) * MolMass(model.metFormulas{strcmp(model.mets,'glycogen_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0072'})) * MolMass(model.metFormulas{strcmp(model.mets,'soluble_BAC[e]')})+...
        FluxDistribution(strcmp(worm.rxns,{'EXC0142'})) * MolMass(model.metFormulas{strcmp(model.mets,'phospholipid_BAC[e]')});
    wasteDW = bacWaste / (-FluxDistribution(strcmp(worm.rxns,'EXC0050')) * MolMass(model.metFormulas{strcmp(model.mets,'BAC[e]')}));
    fprintf('the total waste of bulk bacteria is: %f\n',wasteDW);
else
    wasteDW = [];
    fprintf('Please inspect the flux distribution manually for a non-C. elegans model\n');
end
%% fix the output for speed mode level 3
if speedMode == 3
    RLNames = RLNames_ori;% recover RLNames list for output
end
end
