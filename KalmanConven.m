function [LL,S] = KalmanConven(sys,N,w,v,Rww,Rvv)
%% Initialisasi
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
x_hat_kf = .1*ones(size(A,1),N);
u = zeros(size(B,2),1);
%%  Init
% Rww = 0.1*ones(size(A,1),size(A,1)); %Q=Rww
% Rvv = 0.0001; %R=Rvv
P = rand(size(A,1),size(A,1),N); %Matriks covariance state estimate
%% Iteration
for i = 1:N
    x_hat_kf(:,i+1) = A*x_hat_kf(:,i)+B*u + w(:,i);       
    y(:,i) = C*x_hat_kf(:,i)+v(:,i);
    %Prior Error Covariance
    P(:,:,i+1) = A*P(:,:,i)*A'+Rww;
    P1_new(:,:,i) = A*P(:,:,i)+P(:,:,i)*A'-P(:,:,i)*C'*inv(Rvv)*C*P(:,:,i)+Rww;
    %Measurement Update
    S(:,:,i) = C*P1_new(:,:,i)*C'+Rvv;
    %Calculate Kalman Gain 
    LL(:,:,i) = P(:,:,i)*C'*inv(S(:,:,i));
    x_hat_kf(:,i+1) = x_hat_kf(:,i)+LL(:,:,i)*(y(:,i)-C*x_hat_kf(:,i));
end