function [hasil] = convAnalysis(dataset_cp,N,results)
%% Kajian Konvergensi
T = N;
x = dataset_cp.x;
y = dataset_cp.y;
x_hat = dataset_cp.x_hat;
K_lqg = results.K_lqg;
K_kn_vi = results.K_kn_vi;
K_kf_vi = results.K_kf_vi;
K_kn_lqr = results.K_kn_lqr;
%% Build F and Fw via Toeplitz Matrix 
A = results.sys.A;
B = results.sys.B;
nrow = size(A,2);
rmat{1} = (zeros(nrow,nrow));
rmat{2} = (ones(nrow,nrow));
temp = 1; 
for i = 3:N
rmat{i} = (A^temp);
temp=temp+1;
end 
rmatMat = cell2mat(rmat)';
% r = [zeros(nrow,N);ones(nrow,N)];
% for i = size(r,1)+1:N
%     r(:,i) = A^(i-1+1)*B;
% end
c = zeros(nrow,N);
Fw = toeplitz(c,rmatMat);
%%
oCell{1} = ones(nrow,nrow);
for i = 2:N
    oCell{i} = A^i;
end
oMat = cell2mat(oCell)'; 
%%
frmat{1} = (zeros(nrow,size(B,2)));
frmat{2} = B;
(ones(nrow,nrow));
ftemp = 1; 
for i = 3:N
frmat{i} = (A^ftemp)*B;
ftemp=ftemp+1;
end 
frmatMat = cell2mat(frmat)';
fc = zeros(nrow,N);
Fu = toeplitz(fc,frmatMat);
F = [oMat Fu];
%%
u_kn_vi = results.u_kn_vi;
u_kf_vi = results.u_kf_vi;
u_kn_lqr = results.u_kn_lqr;

x_hat_kn_vi = results.x_hat_kn_vi;
x_hat_kn_lqr = results.x_hat_kn_lqr;
x_hat_kf_vi = results.x_hat_kf_vi;

n = size(x,2);
m = size(y,2);
sig = .3;
scov = sqrt(m*((n+m)*T+n)*log(9/sig));
x_nw = dataset_cp.x_nw;

deltax_kn_vi = x_nw-x_hat_kn_vi; 
deltax_kn_lqr = x_nw-x_hat_kn_lqr;
deltax_kf_vi = x_nw-x_hat_kf_vi;


Rww = results.Rww;
Qw = kron(eye(size(Rww,1)),Rww);

sigmaU1 = cov(u_kn_vi);
sigmaU2 = cov(u_kf_vi);
sigmaU3 = cov(u_kn_lqr);

a = 24*norm(Fw,2)*norm(Qw,2)^.5;
b1 = rank(deltax_kn_vi);
b2 = rank(deltax_kf_vi);
b3 = rank(deltax_kn_lqr);

GammaU1 = 4*a*sqrt(b1*n*(T+1));
GammaU2 = 4*a*sqrt(b2*n*(T+1));
GammaU3 = 4*a*sqrt(b3*n*(T+1));


GammaX1 = (norm(F,2)*GammaU1)+((16*norm(Fw,2)*norm(sigmaU1,2)^.5*norm(Qw,2))*(norm(u_kn_vi,2)+GammaU1));
GammaX2 = (norm(F,2)*GammaU2)+((16*norm(Fw,2)*norm(sigmaU2,2)^.5*norm(Qw,2))*(norm(u_kf_vi,2)+GammaU2));
GammaX3 = (norm(F,2)*GammaU3)+((16*norm(Fw,2)*norm(sigmaU3,2)^.5*norm(Qw,2))*(norm(u_kn_lqr,2)+GammaU3));

GammaK1 = (GammaU1+GammaX1)*scov;
GammaK2 = (GammaU2+GammaX2)*scov;
GammaK3 = (GammaU3+GammaX3)*scov;

%%
[~, S1, ~] = svd(dataset_cp.x-x_hat);
[U, S, V] = svd(x_hat);


max_singular_value = max(diag(S1));
min_singular_value = min(diag(S));
kappaX = max_singular_value/min_singular_value;
tempNilai = 1/(min_singular_value*(1-kappaX));


deltaK_kn_vi = norm(K_lqg-K_kn_vi,2);
deltaK_kf_vi = norm(K_lqg-K_kf_vi,2);
deltaK_kn_lqr = norm(K_lqg-K_kn_lqr,2);

m = size(u_kn_vi,2);

rho=0.2;
deltaKacu_kn_vi = abs(tempNilai)*(deltaK_kn_vi/sqrt(n)+(sqrt(m))*rho');
deltaKacu_kf_vi = abs(tempNilai)*(deltaK_kf_vi/sqrt(n)+(sqrt(m))*rho');
deltaKacu_kn_lqr = abs(tempNilai)*(deltaK_kn_lqr/sqrt(n)+(sqrt(m))*rho');
%% 
hasil.kappaX = kappaX;
hasil.deltaKacu_kn_vi = deltaKacu_kn_vi; 
hasil.deltaKacu_kf_vi = deltaKacu_kf_vi; 
hasil.deltaKacu_kn_lqr = deltaKacu_kn_lqr; 
hasil.deltaK_kn_vi = deltaK_kn_vi;
hasil.deltaK_kf_vi = deltaK_kf_vi;
hasil.deltaK_kn_lqr = deltaK_kn_lqr;
end