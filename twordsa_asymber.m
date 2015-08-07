% Prashanth L.A., Jul. 2015
%
% An RDSA variant of the 2SPSA code. The primary difference is in the generation of perturbation r.v.s. 
% In this case, the latter are sampled from an asymmetric Bernoulli
% distribution, i.e., each component of the perturbation takes value -1
% with probability (1+epsilon)/(2+epsilon) and (1+epsilon) with probability
% 1/(2+epsilon)
%
% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% epsilon -> distribution parameter for asymmetric Bernoulli perturbations.
% numSimulation -> this is the simulation budget that impacts the number of 2RDSA-AsymBer iterations
% replications -> number of independent simulations
% theta_0 -> initial point for 1RDSA-AsymBer (If N=0, then 2RDSA-AsymBer starts at theta_0)
%
function [w x y z] = twordsa_asymber(p, sigma, type, epsilon, numSimulations,replications, theta_0)
%value of numerator in a_k sequence for all iterations of 1RDSA-AsymBer 
%and first N-measurement-based iterations in the initialization of 2RDSA-AsymBer
a1=1;
%value of numerator in a_k in 2RDSA-AsymBer part
a2=10;
A1=50;           %stability constant for 1RDSA-AsymBer
A2=0;           %stability constant for 2RDSA-AsymBer
c1=1.9;         %numerator in c_k for 1RDSA-AsymBer 
c2=2*1.9;         %numerator in c_k for 2RDSA-AsymBer
alpha1=1;      %a_k decay rate for 1RDSA-AsymBer 
alpha2=0.6;         %a_k decay rate for 2RDSA-AsymBer
gamma1=.101;      %c_k decay rate for 1RDSA-AsymBer
gamma2=.1666701;        %c_k decay rate for 2RDSA-AsymBer
epsilon1=0.0001;
N=0.2*numSimulations;				%no. of function meas. for 1RDSA-AsymBer initialization	
loss='loss_myexample_noise';	  %loss function for use in algorithm (usually with noise)
lossfinaleval='loss_myexample'; %loss function for "true" evaluation of algorithm (no noise)

% Loss value for the initial point
Ltheta0 = feval(lossfinaleval, p, theta_0, type);
if(Ltheta0==0) 
    Ltheta0=1; 
end
    
% the optima is problem-dependent. For quadratic losses,
% thetaStar(i,1)=-0.9091 for all i=1,..,p, For fourth-order loss,
% thetaStar=0
thetaStar = getOptima(p, type);

rand('seed',31415297)
randn('seed',3111113)                       
delta = zeros(p,1);
%
%the loop 1:replications below is for doing multiple replications for use in averaging to 
%evaluate the relative performance.  
%
%the first loop in 2RDSA-AsymBer below uses the standard 1RDSA-AsymBer form to initialize 2RDSA-AsymBer
%
%the second loop does 2RDSA-AsymBer following guidelines in Spall ASP (Chap. 7 of ISSO))
%
%lines below initialize various recursions for the gradient/Hess. averaging
%and for final error reporting based on the average of the solutions for 
%"replications" replications.
%
meanHbar=0;
errtheta=0;
losstheta=0;				%cum. sum of loss values
theta_lo=-2.048*ones(p,1);   %lower bounds on theta  
theta_hi=2.047*ones(p,1);    %upper bounds on theta 

nmseAllReplications = zeros(1,p);
lossesAllReplications = zeros(1,p);
mseTheta0=(theta_0-thetaStar)'*(theta_0-thetaStar);

%DUMMY STATEMENT FOR SETTING DIMENSIONS OF Hhat (AVOIDS OCCASIONAL
%ERROR MESSAGES)
Hhat=eye(p);
for j=1:replications
%INITIALIZATION OF PARAMETER AND HESSIAN ESTIMATES
  theta=theta_0;
  Hbar=500*eye(p);
  %INITIAL N ITERATIONS OF 1RDSA-AsymBer PRIOR TO 2RDSA-AsymBer ITERATIONS
  for k=1:N/2    %use of N-avg is to account for avg used in setting lossold 
    a_k=a1/(k+A1)^alpha1;
    c_k=c1/k^gamma1; 
    unifrands = unifrnd(0,1,p,1);
    for l=1:p
        if unifrands(l,1) < ((1+epsilon1)/(2+epsilon1))
        delta(l,1) = -1;
        else
            delta(l,1) = 1+epsilon1;
        end
    end
    thetaplus=theta+c_k*delta;
    thetaminus=theta-c_k*delta;  
    yplus=feval(loss, p, thetaplus, sigma, type);
    yminus=feval(loss, p, thetaminus, sigma, type);
    ghat=(yplus-yminus)./(2*c_k*delta);   
%   theta update
    theta=theta-a_k*ghat;
    % Two lines below invoke constraints
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
  end
%
% START 2RDSA-AsymBer ITERATIONS FOLLOWING INITIALIZATION
%  
for k=1:(numSimulations-N)/3
    a_k=a2/(k+A2)^alpha2;
    c_k=c2/k^gamma2;
% GENERATION OF AVERAGED GRADIENT AND HESSIAN (NO AVERAGING IF gH_avg=1)                     
      unifrands = unifrnd(0,1,p,1);
        for l=1:p
            if unifrands(l,1) < ((1+epsilon)/(2+epsilon))
            delta(l,1) = -1;
            else
                delta(l,1) = 1+epsilon;
            end
        end
      thetaplus=theta+c_k*delta;
      thetaminus=theta-c_k*delta;  
      yplus=feval(loss, p, thetaplus, sigma, type);
      yminus=feval(loss, p, thetaminus, sigma, type);      
      ghat=(1/(1+epsilon))*(yplus-yminus)/(2*c_k)*delta;  
% GENERATE THE HESSIAN UPDATE
      M_n = delta*delta'*(1/(2*(1+epsilon)^2));
      beta = (1+epsilon)*(1+ (1+epsilon)^3)/(2+epsilon);
      kappa = beta * (1 - (1+epsilon)^2/beta);
      for idx=1:p
          M_n(idx,idx) = 1/kappa*((2*(1+epsilon)^2)*M_n(idx,idx) - (1+epsilon));
      end      
% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS. PER ITERATION      
      y=feval(loss, p, theta, sigma, type);      
      Hhat = ((yplus+yminus-2*y)/(c_k^2))*M_n;
      Hhat=.5*(Hhat+Hhat');    
    Hbar=(k/(k+1))*Hbar+Hhat/(k+1);          
%   THE THETA UPDATE (FORM BELOW USES GAUSSIAN ELIMINATION TO AVOID DIRECT 
%   COMPUTATION OF HESSIAN INVERSE)
    Hbarbar=sqrtm(Hbar*Hbar+.000001*eye(p)/k);    
    % The main update step
    theta=theta-a_k*(Hbarbar\ghat);
    % Two lines below invoke constraints
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
 end
meanHbar=meanHbar+Hbar;
errtheta=errtheta+(theta-thetaStar)'*(theta-thetaStar); 
losstheta=losstheta+feval(lossfinaleval, p, theta, type);
lossesAllReplications(1,j) = feval(lossfinaleval, p, theta, type);
nmseAllReplications(1, j) = (theta-thetaStar)'*(theta-thetaStar)/mseTheta0;
end

% Return Normalized loss and mean square error, with sample standard
% deviations
w=losstheta/replications/Ltheta0;
x=std(lossesAllReplications);
y=errtheta/replications/mseTheta0;
z=std(nmseAllReplications);
%disp(mat2str(theta,4));
