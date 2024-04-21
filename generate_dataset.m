%% Generate dataset 
clear all;
clc
Ry = 10;
Ru = .1;
Rww = 10;
Rvv = 10;
N = 100; 
df = .1;
% Set the model parameter :
% Ry : weight state matrices
% Ru : weight consig matrices
% Rww : weight of process error covariance 
% Rvv : weight of measurement noise covariance 
% df : discount factor 
% N : number of sample-data
param.Ry = Ry;
param.Ru = Ru;
param.Rww = Rww;
param.Rvv = Rvv;
param.N = N;
param.df = df; 
%%
plant = sysmdl_unstable(N,df)