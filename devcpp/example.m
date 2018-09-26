% example running script
clear all;
clc;

load('TestData.mat');
lb = [-1 -1 0 0 -1 -1 0 0];
ub = [+1 +1 0 0 +1 +1 0 0];

mex -I"O:\Array Benchmark\devcpp" DE_CSM.cpp DESolver.cpp
[out1, out2] = DE_CSM(frequencies, c, real(CSM), imag(CSM), lb, ub, mic_info);

%%
clear all; clc;
mex -I"O:\Array Benchmark\devcpp" DE_CSM.cpp DESolver.cpp
CsmReader;

f = 500;
lb = [-2 2 0  0];
ub = [2 2 4 2];
tic;
[out1, out2] = DE_CSM(f, csound, double(squeeze(cpreal(1,:,:))), double(squeeze(cpimag(1,:,:))), lb, ub, double(mic_info.'));
toc

%% Energy fun check
format long;
Energy_fun_C(ub,cpreal,cpimag,1,500,csound,mic_info(1,:),mic_info(2,:),mic_info(3,:))
format short;