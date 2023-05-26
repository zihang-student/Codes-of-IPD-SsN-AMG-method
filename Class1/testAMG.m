


load('wrongAMG.mat');

amg_options.maxit=50;

[zeta,itpcg,respcg,info] = Hybrid_AMG(prob_data,amg_options);
