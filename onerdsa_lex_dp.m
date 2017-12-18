% Prashanth L.A., Jul. 2015
%
% An RDSA variant of the 1SPSA code. The primary difference is in the generation of perturbation r.v.s. 
% In this case, the latter are sampled from an asymmetric Bernoulli distribution.
% Deterministic Perturbations are generated from lexicographic sequence.

% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% epsilon ->
% numSimulation -> this is the simulation budget that impacts the number of iterations
% replications -> number of independent simulations
% theta_0 -> initial point 
%
function [w x y z] = onerdsa_lex_dp(p, sigma, type, numSimulations, replications, theta_0)
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
  for k=0:(numSimulations/(2*3^p))-1
    ak = a/(k+1+A)^alpha;
    ck = c/(k+1)^gamma;
    
    % Generating lexicograpic sequence  
    delta = zeros(3^p,p);
    for t = 1:p  
       temp = [-1*ones(2*3^(p-t),1); 2*ones(3^(p-t),1)]; 
       delta(:,t) = repmat(temp,3^(t-1),1);
    end
    ghat = 0;
    
    for j = 1:3^p
        thetaplus = theta + ck*delta(j,:)';
        thetaminus = theta - ck*delta(j,:)';
        yplus=feval(loss, p, thetaplus, sigma, type);
        yminus=feval(loss, p, thetaminus, sigma, type); 
        ghat = ghat + delta(j,:)'*((yplus - yminus)/(2*ck));
    end
    
    ghat = ghat/(2*3^p);
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

disp(['Number of iterations of outer for loop is : ',num2str(k+1)]);


% Display normalized loss values
%sprintf('Normalized loss: %5.4f', lossfinal/replications/Ltheta0)

% Display results: MSE normalized
str = sprintf('Normalised MSE: %10.9f, Std dev: %10.9f',errtheta/replications/mseTheta0, std(nmseAllReplications));
disp(str);
str1 = sprintf('Std Error: %10.9f',std(nmseAllReplications)/sqrt(p));
disp(str1);
w=lossfinal/replications/Ltheta0;
x=std(lossesAllReplications);
y=errtheta/replications/mseTheta0;
z=std(nmseAllReplications);
disp(mat2str(theta,4));
%disp(quantile(nmseAllReplications, [0.95 0.05]));
