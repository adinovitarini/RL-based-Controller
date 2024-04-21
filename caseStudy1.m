function [results resultConv params] = caseStudy1(params)
N = params.N ;
Rww = params.Rww;
Rvv = params.Rvv;
Q = params.Q;
R = params.R; 
df = params.df; 
nh = params.nh;
bs = params.bs;
dl = params.dl;
cartpole = sysmdl_cartpole(N,df);
% cartpole = sysmdl_unstable(N,df);
init = .1*ones(size(cartpole.sys.A,1),N);
dataset_cp = GenerateSeq(cartpole.sys,N,Rww,Rvv,Q,R);
dataset_cp_ts = GenerateSeq(cartpole.sys,N,Rww,Rvv,Q,R);
delta_x_cp = [dataset_cp.x(:,1:end-1);dataset_cp.y];
target_cp = [dataset_cp.x_nw(:,1:end-1);dataset_cp.y_nw];
u_cp = dataset_cp.u(:,1:end-1);
%%
tic;
% (delta_x,target,u,A,B,C,N)
[KG_cp,net_cp] = KalmanNet(delta_x_cp,target_cp,u_cp,cartpole.sys.A,cartpole.sys.B,cartpole.sys.C,N,nh,bs,dl);
time_kn_cp = toc;
%%  Implement KG for state estimation process (u:LQR)
delta_y_cp = dataset_cp.y-dataset_cp.y_nw;
x_hat_net_nw_cp = (KG_cp*delta_y_cp')+dataset_cp.x_nw(:,1:end-1);
x_hat_net_cp = (KG_cp*delta_y_cp')+(dataset_cp.x(:,1:end-1)-dataset_cp.x_nw(:,1:end-1));
y_hat_net_cp = cartpole.sys.C*x_hat_net_cp;
%% Four scenario Test 
cp = cartpole.sys;
%%  First scenario : KF-VI   
[kf_vi] = combine_kf_vi(cp.A,cp.B,cp.C,Q,R,N,dataset_cp.L,dataset_cp.y);
P_kf_vi = kf_vi.P;
u_kf_vi = kf_vi.u; 
x_hat_kf_vi = kf_vi.x_hat; 
y_kf_vi = kf_vi.y;
J_kf_vi = kf_vi.J;
K_kf_vi = kf_vi.K;
L_kf_vi = kf_vi.L;
%%  Second scenario : KalmanNet-LQR   
[kn_lqr] = combine_kn_lqr(cp.A,cp.B,cp.C,Q,R,N,x_hat_net_cp);
P_kn_lqr = kn_lqr.P;
u_kn_lqr = kn_lqr.u; 
x_hat_kn_lqr = kn_lqr.x_hat; 
y_kn_lqr = kn_lqr.y;
J_kn_lqr = kn_lqr.J;
K_kn_lqr = kn_lqr.K;
%%  Third scenario : KalmanNet-VI
[kn_vi] = combine_kn_vi(cp.A,cp.B,cp.C,Q,R,N,x_hat_net_cp);
P_kn_vi = kn_vi.P;
u_kn_vi = kn_vi.u; 
x_hat_kn_vi = kn_vi.x_hat; 
y_kn_vi = kn_vi.y;
J_kn_vi = kn_vi.J;
K_kn_vi = kn_vi.K;
%% Fourth scenario : LQG konvensional 
tic
x_lqg = dataset_cp.x(:,1:end-1);
y_lqg = dataset_cp.y;
L_lqg = dataset_cp.L;
[K_lqg,~,~] = lqr(cp,Q,R);
for i = 1:N
    u_lqg(i) = K_lqg*x_lqg(:,i);
    x_lqg(:,i+1) = (cp.A-cp.B*K_lqg)*x_lqg(:,i)+L_lqg(i)*(y_lqg(i)-cp.C*x_lqg(:,i));
    y_lqg(i) = cp.C*x_lqg(:,i);
end
l_time = toc;
J_lqg = value_func(x_lqg,u_lqg,Q,R,N);
%%
results.u_kn_vi = u_kn_vi;
results.u_kn_lqr = u_kn_lqr;
results.u_kf_vi = u_kf_vi;
results.y_kn_vi = y_kn_vi;
results.y_kn_lqr = y_kn_lqr;
results.y_kf_vi = y_kf_vi;
results.x_hat_kn_vi = x_hat_kn_vi;
results.x_hat_kn_lqr = x_hat_kn_lqr;
results.x_hat_kf_vi = x_hat_kf_vi;
results.u_lqg = u_lqg;
results.y_lqg = y_lqg;
results.x_hat_lqg = x_lqg;
results.sys = cp;
results.Rww = Rww;
results.Rvv = Rvv; 
results.Q = Q;
results.R = R;
results.K_kn_vi = K_kn_vi;
results.K_kn_lqr = K_kn_lqr;
results.K_kf_vi = K_kf_vi;
results.K_lqg = K_lqg;
results.kn_vi = kn_vi;
results.kf_vi = kf_vi; 
results.kn_lqr = kn_lqr; 

%%
[resultConv] = convAnalysis(dataset_cp,N,results);
