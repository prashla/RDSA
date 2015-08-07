% J.C. Spall, Aug. 1998
% twospsaconstrained
% Code for evaluation of second-order SPSA (2SPSA) versus first-order 
% SPSA (1SPSA, as in Chap. 7 of ISSO).  Code is for comparative evaluation purposes; hence,  
% it includes much that is not required for a basic implementation. Further, it is 
% in no way "optimized" for efficiency or generality; this is strictly research
% code for the purpose of getting a basic idea of how 2SPSA works.

% Code includes the capability for initializing 2SPSA by running 1SPSA for 
% N measurements.  This code allows for averaging of the 
% SP gradients and Hessian estimates at EACH iteration after the initial (N) 
% measurements where only 1SPSA is used for estimating theta.  We use "theta"
% for the 2SPSA recursion.  Code allows for checking for simple constraint 
% violation (componentwise constraints).
%
% UPDATE MAR. 2006: Feedback and weighting versions of this code are available from the author;
% this can generally provide enhanced performance. (Reference: Spall, J. C. (2006), â€œFeedback and Weighting Mechanisms 
% for Improving Jacobian (Hessian) Estimates in the Adaptive Simultaneous Perturbation Algorithm,â€? Proceedings of the 
% American Control Conference, 14-16 June 2006, Minneapolis, MN, paper ThB09.1 in CD-ROM.)
%
% Prashanth L.A., Jul. 2015
% The code below is from J.C. Spall's ISSO book website, except that it
% does not perform any "blocking" test for the parameter.
% 
% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations
% replications -> number of independent simulations
% theta_0 -> initial point for 1SPSA (If N=0, then 2SPSA starts at theta_0)

function twospsa(p, sigma, type, numSimulations, replications, theta_0)
%value of numerator in a_k sequence for all iterations of 1SPSA 
%and first N-measurement-based iterations in the initialization of 2SPSA
a1=1;
%value of numerator in a_k in 2SPSA part
a2=10;
A1=50;           %stability constant for 1SPSA
A2=0;           %stability constant for 2SPSA
c1=1.9;         %numerator in c_k for 1SPSA 
c2=2*1.9;         %numerator in c_k for 2SPSA
ctilda=2*1.9;     %numerator in ctilda_k for 2SPSA;
alpha1=1;      %a_k decay rate for 1SPSA 
alpha2=0.6;         %a_k decay rate for 2SPSA
gamma1=.101;      %c_k decay rate for 1SPSA
gamma2=.1666701;        %c_k decay rate for 2SPSA
N=0.2*numSimulations;				%no. of function meas. for 1SPSA initialization	
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

%
%the loop 1:replications below is for doing multiple replications for use in averaging to 
%evaluate the relative performance.  
%
%the first loop in 2SPSA below uses the standard 1SPSA form to initialize 2SPSA
%
%the second loop does 2SPSA following guidelines in Spall ASP (Chap. 7 of ISSO))
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


lossesAllReplications = zeros(1,p);
nmseAllReplications = zeros(1,p);
mseTheta0=(theta_0-thetaStar)'*(theta_0-thetaStar);
%DUMMY STATEMENT FOR SETTING DIMENSIONS OF Hhat (AVOIDS OCCASIONAL
%ERROR MESSAGES)
Hhat=eye(p);
for j=1:replications
%INITIALIZATION OF PARAMETER AND HESSIAN ESTIMATES
  theta=theta_0;
  Hbar=500*eye(p);
  %INITIAL N ITERATIONS OF 1SPSA PRIOR TO 2SPSA ITERATIONS
  for k=1:N/2    %use of N-avg is to account for avg used in setting lossold 
    a_k=a1/(k+A1)^alpha1;
    c_k=c1/k^gamma1; 
    delta=2*round(rand(p,1))-1;
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
% START 2SPSA ITERATIONS FOLLOWING INITIALIZATION
%  
for k=1:(numSimulations-N)/4
%for k=1:(numSimulations-N)/4
    a_k=a2/(k+A2)^alpha2;
    c_k=c2/k^gamma2;
    ctilda_k=ctilda/k^gamma2;
% GENERATION OF AVERAGED GRADIENT AND HESSIAN (NO AVERAGING IF gH_avg=1)                 
      delta=2*round(rand(p,1))-1;
      thetaplus=theta+c_k*delta;
      thetaminus=theta-c_k*delta;  
      yplus=feval(loss, p, thetaplus, sigma, type);
      yminus=feval(loss, p, thetaminus, sigma, type);
      ghat=(yplus-yminus)./(2*c_k*delta);  
% GENERATE THE HESSIAN UPDATE
      deltatilda=2*round(rand(p,1))-1;
      thetaplustilda=thetaplus+ctilda_k*deltatilda;
      thetaminustilda=thetaminus+ctilda_k*deltatilda;
% LOSS FUNCTION CALLS      
      yplustilda=feval(loss, p, thetaplustilda, sigma, type);
      yminustilda=feval(loss, p, thetaminustilda, sigma, type);
      ghatplus=(yplustilda-yplus)./(ctilda_k*deltatilda);
      ghatminus=(yminustilda-yminus)./(ctilda_k*deltatilda);
% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS. PER ITERATION
      deltaghat=ghatplus-ghatminus;
      for i=1:p
        Hhat(:,i)=deltaghat(i)./(2*c_k*delta);
      end
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

% Display results: normalized loss and normalized mean square error, 
% both with sample standard deviation
str = sprintf('Normalized loss: %e +- %e, Normalised MSE: %e +- %e', std(lossesAllReplications), losstheta/replications/Ltheta0, errtheta/replications/mseTheta0, std(nmseAllReplications));
disp(str);
%disp(mat2str(theta,4));
