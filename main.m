clear all;
close all;
clc;


%% Select a sampling type =================================================
param.typ = 'VISTA'; % VISTA
% param.typ = 'VRS'; % Variable density random
% param.typ = 'UIS'; % Uniform interleaved


%% Select appropriate parameter values ====================================
param.p     = 96; % Number of phase encoding steps
param.t     = 32; % Number of frames
param.R     = 8;  % Net acceleration rate
param.alph  = 0.28;      % 0<alph<1 controls sampling density; 0: uniform density, 1: maximally non-uniform density
param.sig   = param.p/5; % Std of the Gaussian envelope for sampling density
param.sd    = 10; % Seed to generate random numbers; a fixed seed should reproduce VISTA


%% Probably you don't need to mess with these paramters ===================
% If unsure, leave them empty and the default value will be employed.
param.nIter= []; % Number of iterations for VISTA (defualt: 120)
param.ss   = []; % Step-size for gradient descent. Default value: 0.2; 
param.tf   = []; % Step-size in time direction wrt to phase-encoding direction; use zero for constant temporal resolution. Default value: 0.0
param.s    = []; % Exponent of the potenital energy term. Default value 1.4
param.g    = []; % Every gth iteration is relocated on a Cartesian grid. Default value: floor(param.nIter/6)
param.uni  = []; % At param.uni iteration, reset to equivalent uniform sampling. Default value: floor(param.nIter/2)
param.W    = []; % Scaling of time dimension; frames are "W" units apart. Default value: max(param.R/6,1)
param.sz   = []; % Display size of samples. Default value: 3.5
param.dsp  = 5; % Display frequency (verbosity), every dsp-th iteration will be displayed. Default value: 1
param.fs   = 1; % Does time average has to be fully sampled, 0 for no, 1 for yes. Only works with VISTA. Default value: 1
param.fl   = []; % Start checking fully sampledness at fl^th iteration. Default value: floor(param.nIter*5/6)

                  
                 
%% Check parameters values
param = checkParam(param);


%% Call VISTA to compute the 2D sampling pattern
samp = VISTA(param);


%% Plot the results
plotSamp(samp, param);


%% Save sampling
% save(['samp_', param.typ, '_', num2str(param.p),'x',num2str(param.t) '_R', num2str(param.R)], 'samp');
