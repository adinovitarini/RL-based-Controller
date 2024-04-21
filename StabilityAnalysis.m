%% Stability Analysis
function Phi = StabilityAnalysis(sys,K_lstm,N)
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
% K_lstm = result.drqn.k;
%%
for i = 1:N
    Lambda(:,:,i) = A-B*K_lstm(i,:);
    norm_Lambda (i) = norm(Lambda(:,:,i)); 
end
temp(:,:,1) = Lambda(:,:,1);
for i = 2:N
    temp(:,:,i) = Lambda(:,:,i)*temp(:,:,i-1);
end
%%
Phi.Lambda = Lambda;
Phi.norm_Lambda = norm_Lambda;
Phi.phi = temp(:,:,end);
Phi.norm_phi = norm(Phi.phi);