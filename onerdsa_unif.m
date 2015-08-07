% Prashanth L.A., Jul. 2015
%
% An RDSA variant of the 1SPSA code from J.C. Spall. The primary difference is in the generation of perturbation r.v.s. 
% In this case, the latter are sampled from an asymmetric Bernoulli
% distribution. Moreover, the implementation is made very basic, with many
% of the original checks, e.g. blocking, being removed.
%
% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations
% replications -> number of independent simulations
% theta_0 -> initial point
%
function onerdsa_unif(p, sigma, type, numSimulations, replications, theta_0)
alpha =1;
gamma =.101;
a=1;
c=1.9;      				%chosen by standard guidelines
A=50;
errtheta=0;
lossfinal=0;            %variable for cum. loss values
theta_lo=-2.048*ones(p,1);   %lower bounds on theta  
theta_hi=2.047*ones(p,1);    %upper bounds on theta 

lossfinaleval='loss_myexample';  %choice of loss function for final perf. evaluations (noise-free)                            % evaluation
loss='loss_myexample_noise';     %loss function used in algorithm operations (with noise)
rand('seed',61415927)
randn('seed',6111113)

Ltheta0 = feval(lossfinaleval, p, theta_0, type);
if(Ltheta0==0) 
    Ltheta0=1; 
end

% get optimal theta based on type
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
    % Generate uniform [-1,1] perturbations
    delta = unifrnd(-1,1,p,1);
    thetaplus = theta + ck*delta;
    thetaminus = theta - ck*delta;
    yplus=feval(loss, p, thetaplus, sigma, type);
    yminus=feval(loss, p, thetaminus, sigma, type);
    ghat = 3*((yplus - yminus)/(2*ck))*delta;
    theta=theta-ak*ghat;
    % Two lines below invoke component-wise constraints
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
  end
  lossvalue=feval(lossfinaleval, p, theta, type);
  lossfinal=lossfinal+lossvalue;
  errtheta=errtheta+(theta-thetaStar)'*(theta-thetaStar); 
  lossesAllReplications(1, i) = lossvalue/Ltheta0;
  nmseAllReplications(1, i) = (theta-thetaStar)'*(theta-thetaStar)/mseTheta0;
end

% Display results: Normalized loss and mean square error, with sample standard deviations
str = sprintf('Normalized loss: %e +- %e, Normalised MSE: %e +- %e', lossfinal/replications/Ltheta0, std(lossesAllReplications), errtheta/replications/mseTheta0, std(nmseAllReplications));
disp(str);
%disp(mat2str(theta,4));
