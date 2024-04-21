clear all;clc
disp('Case 1 : cart-pole');
disp('Case 2 : batch-distillation'); 
disp('Case 3 : unstable system'); 



n = input('Enter a number that define your plant: ');


disp('====================Parameter KalmanNet==========================')
nh = input('jumlah hidden unit : ');
bs = input('jumlah batch size : '); 
dl = input('nilai dropout : ');
disp('=====================Parameter Gangguan===========================')
Rww = input('Covariance gangguan : ');
Rvv = input('Covariance gangguan : ');
Q = input('Covariance gangguan : ');
R = input('Covariance gangguan : ');
df = input('Faktor Diskon : ');
N = input('Jumlah trayektori : ');
params.df = df;
params.N = N;
params.Q = Q;
params.R = R;
params.Rww = Rww;
params.Rvv = Rvv;
params.nh = nh; 
params.bs = bs;
params.dl = dl;

switch n 
    case 1
        % disp('1st case study : cart-pole system');
        [results,resultConv,params] = caseStudy1(params);
        aa = input('Save the results as (''ex:case1.mat'') : ');
        save aa results resultConv params
    case 2
        % disp('2nd case study : batch distillation')
        [results,resultConv,params] = caseStudy2(params);
        aa = input('Save the results as (''ex:case2.mat'') : ');
        save aa results resultConv params
    case 3
        % disp('3rd case study : unstable system')
        [results,resultConv,params] = caseStudy3(params);
        aa = input('Save the results as (''ex:case3.mat'') : ');
        save aa results resultConv params
end