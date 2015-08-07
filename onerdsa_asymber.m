% Prashanth L.A., Jul. 2015
%
% An RDSA variant of the 1SPSA code. The primary difference is in the generation of perturbation r.v.s. 
% In this case, the latter are sampled from an asymmetric Bernoulli distribution.
% 
% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% epsilon ->
% numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations
% replications -> number of independent simulations
% theta_0 -> initial point 
%
function [w x y z] = onerdsa_asymber(p, sigma, type, epsilon, numSimulations, replications, theta_0)
% the following are chosen by standard guidelines
alpha =1; % exponent for stepsize
gamma =.101; % exponent for perturbation constant
a=1; % initial step-size value
c=1.9; % initial perturbation constant
A=50;

errtheta=0;
lossfinal=0;            %variable for cum. loss values
theta_lo=-2.048*ones(p,1);   %lower bounds on theta  
theta_hi=2.047*ones(p,1);    %upper bounds on theta 

delta = zeros(p,1);
lossfinaleval='loss_myexample';  %choice of loss function for final perf. evaluations (noise-free)                            % evaluation
loss='loss_myexample_noise';     %loss function used in algorithm operations (with noise)
rand('seed',61415927)
randn('seed',6111113)

% set optimal theta based on type
Ltheta0 = feval(lossfinaleval, p, theta_0, type);
if(Ltheta0==0) 
    Ltheta0=1; 
end

thetaStar = getOptima(p, type);

lossesAllReplications = zeros(1,p);
nmseAllReplications = zeros(1,p);
mseTheta0=(theta_0-thetaStar)'*(theta_0-thetaStar);
% outer loop for replications
for i=1:replications
  theta=theta_0;
  for k=0:numSimulations/2-1
    ak = a/(k+1+A)^alpha;
    ck = c/(k+1)^gamma;
    % Generate asymmetric Bernoulli perturbations
    unifrands = unifrnd(0,1,p,1);
    for j=1:p
        if unifrands(j,1) < ((1+epsilon)/(2+epsilon))
            delta(j,1) = -1;
        else
            delta(j,1) = 1+epsilon;
        end
    end
    thetaplus = theta + ck*delta;
    thetaminus = theta - ck*delta;
    yplus=feval(loss, p, thetaplus, sigma, type);
    yminus=feval(loss, p, thetaminus, sigma, type);
    ghat = (1/(1+epsilon))*((yplus - yminus)/(2*ck))*delta;
    theta=theta-ak*ghat;
    % Two lines below invoke component-wise constraints
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
  end
  lossvalue=feval(lossfinaleval, p, theta, type);
  lossfinal=lossfinal+lossvalue;
  errtheta=errtheta+(theta-thetaStar)'*(theta-thetaStar);
  nmseAllReplications(1, i) = (theta-thetaStar)'*(theta-thetaStar)/mseTheta0;
  lossesAllReplications(1, i) = lossvalue/Ltheta0;
end

% Display normalized loss values
%sprintf('Normalized loss: %5.4f', lossfinal/replications/Ltheta0)

% Display results: MSE normalized
%sprintf('Normalised MSE: %5.4f, Std dev: %5.4f',errtheta/replications/mseTheta0, std(nmseAllReplications))
w=lossfinal/replications/Ltheta0;
x=std(lossesAllReplications);
y=errtheta/replications/mseTheta0;
z=std(nmseAllReplications);
disp(mat2str(theta,4));
%disp(quantile(nmseAllReplications, [0.95 0.05]));
