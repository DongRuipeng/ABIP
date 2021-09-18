clear;
clc; 
close all;

Probname = {'25FV47', '80BAU3B', 'ADLITTLE', 'AFIRO', 'AGG', 'AGG2', 'AGG3', ...
    'BANDM', 'BEACONFD', 'BLEND', 'BNL1', 'BNL2', 'BOEING1', 'BOEING2', ...
    'BORE3D', 'BRANDY', 'CAPRI', 'CRE_A', 'CRE_B', 'CRE_C', 'CRE_D', ...
    'CYCLE', 'CZPROB', 'D2Q06C', 'D6CUBE', 'DEGEN2', 'DEGEN3', 'DFL001', ...
    'E226', 'ETAMACRO', 'FFFFF800', 'FINNIS', 'FIT1D', 'FIT1P', 'FIT2D', ...
    'FIT2P', 'FORPLAN', 'GANGES', 'GFRD_PNC', 'GREENBEA', 'GREENBEB', 'GROW7', ...
    'GROW15', 'GROW22', 'ISRAEL', 'KB2', 'KEN_7', 'KEN_11', 'KEN_13', ...
    'KEN_18', 'LOTFI', 'MAROS', 'MAROS_R7', 'MODSZK1', 'NESM', 'OSA_07', ...
    'OSA_14', 'OSA_30', 'OSA_60', 'PDS_02', 'PDS_06', 'PDS_10', 'PDS_20', ...
    'PEROLD', 'PILOT', 'PILOT4', 'PILOT87', 'PILOT_JA', 'PILOT_WE', 'PILOTNOV', ...
    'QAP8', 'QAP12', 'QAP15', 'RECIPE', 'SC50A', 'SC50B', 'SC105', ...
    'SC205', 'SCAGR7', 'SCAGR25', 'SCFXM1', 'SCFXM2', 'SCFXM3', 'SCORPION', ...
    'SCRS8', 'SCSD1', 'SCSD6', 'SCSD8', 'SCTAP1', 'SCTAP2', 'SCTAP3', ...
    'SEBA', 'SHARE1B', 'SHARE2B', 'SHELL', 'SHIP04L', 'SHIP04S', 'SHIP08L', ...
    'SHIP08S', 'SHIP12L', 'SHIP12S', 'SIERRA', 'STAIR', 'STANDATA', 'STANDGUB', ...
    'STANDMPS', 'STOCFOR1', 'STOCFOR2', 'STOCFOR3', 'TRUSS', 'TUFF', 'VTP_BASE', ...
    'WOOD1P', 'WOODW'};

nprob = length(Probname);

Problist = [1:nprob];

abip = 1; 
abip_time = zeros(nprob, 1); 
abip_ipm_iter = zeros(nprob, 1); 
abip_admm_iter = zeros(nprob, 1);  

scs = 1; 
scs_time = zeros(nprob, 1); 
scs_admm_iter = zeros(nprob, 1);

for di = 1:length(Problist) 
    
    probID = Problist(di);
    name = Probname{probID};
    load(strcat('./netlib/feasible/', Probname{Problist(di)},'.mat'));
    A = Problem.A; 
    b = Problem.b; 
    c = Problem.aux.c; 
    lbounds = Problem.aux.lo; 
    ubounds = Problem.aux.hi; 
    [m, n] = size(A);
    sp = nnz(A)/(m*n);
    
    if scs || abip
       [A,b,c,info] = presolve(A,b,c,lbounds,ubounds);
       
       if scs
          K = struct('l',size(A, 2));
          data.A = sparse(A'); 
          data.c = -full(b); 
          data.b = full(c);
          params_scs = struct('eps', 1e-3, 'alpha', 1.8, 'max_iters', 1000000, 'verbose', 0);
            
          tic; 
          [~, y, ~, info_scs] = scs_direct(data, K, params_scs); 
          time_scs = toc; 
          [~, objp_scs] = postsolve(y, info);
          
          scs_time(di) = time_scs; 
          scs_admm_iter(di) = info_scs.iter;
       end
       
       if abip
          data.A = sparse(A); 
          data.c = full(c); 
          data.b = full(b);
          params_abip = struct('max_admm_iters', 1000000, 'alpha', 1.8, 'adaptive_lookback', 10, 'sparsity_ratio', sp, 'verbose', 0);
            
          tic; 
          [x, ~, ~, info_abip] = abip_indirect(data, params_abip); 
          time_abip = toc; 
          [~, objp_abip] = postsolve(x, info);
          
          abip_time(di) = time_abip; 
          abip_ipm_iter(di) = info_abip.ipm_iter; 
          abip_admm_iter(di) = info_abip.admm_iter; 
       end        
    end
    
    if abip
        fprintf('%10s & %5d & %5d & %3.2e & %3.2e & %3.2e & %3.2e & %5d & %5d & %3.2e\\\\ \\hline \n', ...
            name, m, n, objp_abip, info_abip.resPri, info_abip.resDual, info_abip.relGap, info_abip.ipm_iter, info_abip.admm_iter, time_abip);
    end
    
    if scs
        fprintf('%10s & %5d & %5d & %3.2e & %3.2e & %3.2e & %3.2e & - & %5d & %3.2e\\\\ \\hline \n', ...
            name, m, n, objp_scs, info_scs.resDual, info_scs.resPri, info_scs.relGap, info_scs.iter, time_scs);
    end
end