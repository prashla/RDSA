% D. Sai koti reddy,Feb. 2016
%
% An RDSA variant of the enhanced 2SPSA code. The primary difference is in the generation of perturbation r.v.s. 
% In this case, the latter are sampled from an uniform distribution.
% Code for evaluation of enhanced second-order RDSA_unif, which includes non-identical 
% weighting for the Hessian inputs and feedback.  

% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss
% numSimulation -> this is the simulation budget that impacts the number of 2RDSA-Unif iterations
% replications -> number of independent simulations
% theta_0 -> initial point for 1RDSA-Unif (If N=0, then 2RDSA-Unif starts at theta_0)
%
function y = enhanced_twordsa_unif(p, sigma, type, numSimulations, replications, theta_0)
%value of numerator in a_k sequence for all iterations of 1RDSA-Unif 
%and first N-measurement-based iterations in the initialization of 2RDSA-Unif
a1=1;
%value of numerator in a_k in 2RDSA-Unif part
a2=10;
A1=50;           %stability constant for 1RDSA-Unif
A2=0;           %stability constant for 2RDSA-Unif
c1=1.9;         %numerator in c_k for 1RDSA-Unif 
c2=2*1.9;         %numerator in c_k for 2RDSA-Unif
alpha1=1;      %a_k decay rate for 1RDSA-Unif 
alpha2=0.6;         %a_k decay rate for 2RDSA-Unif
gamma1=.101;      %c_k decay rate for 1RDSA-Unif
gamma2=.1666701;        %c_k decay rate for 2RDSA-Unif
%w=0;           %denominator part of weight sequence
        %numerator in weight sequence for Hessian averaging
d=.501;         %decay rate in weight sequence 
N=0.2*numSimulations;				%no. of function meas. for 1RDSA-Unif initialization	
loss='loss_myexample_noise';	  %loss function for use in algorithm (usually with noise)
lossfinaleval='loss_myexample'; %loss function for "true" evaluation of algorithm (no noise)



% a1=1;
% %value of numerator in a_k in 2SPSA part
% a2=1;
% A1=50;           %stability constant for 1SPSA
% A2=50;           %stability constant for 2SPSA
% c1=1.9;         %numerator in c_k for 1SPSA 
% c2=0.01;         %numerator in c_k for 2SPSA
% alpha1=1;      %a_k decay rate for 1SPSA 
% alpha2=1;         %a_k decay rate for 2SPSA
% gamma1=.101;      %c_k decay rate for 1SPSA
% gamma2=.166667;        %c_k decay rate for 2SPSA
% w=1/10;
% d=0.501;
% N=0.2*numSimulations;				%no. of function meas. for 1SPSA initialization	
% loss='loss_myexample_noise';	  %loss function for use in algorithm (usually with noise)
% lossfinaleval='loss_myexample'; %loss function for "true" evaluation of algorithm (no noise)


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
%the first loop in enhanced 2RDSA-Unif below uses the standard 1RDSA-Unif form to initialize enhanced 2RDSA-Unif
%
%the second loop does 2RDSA-Unif following guidelines in Spall ASP (Chap. 7 of ISSO))
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
  B=triu(ones(p,p))/p;        %matrix used as part of 4th-order loss function
   Hbar=1.05*2*B'*B;
%    Hbar=500*eye(p);
   Hbarbar=Hbar;
   w=0; 
  
  %INITIAL N ITERATIONS OF 1RDSA-Unif PRIOR TO 2RDSA-Unif ITERATIONS
  for k=1:N/2    %use of N-avg is to account for avg used in setting lossold 
    a_k=a1/(k+A1)^alpha1;
    c_k=c1/k^gamma1; 
    % Generate uniform U[-1,1] distributed perturbations
    delta=unifrnd(-1,1,p,1);
    thetaplus=theta+c_k*delta;
    thetaminus=theta-c_k*delta;  
    yplus=feval(loss, p, thetaplus, sigma, type);
    yminus=feval(loss, p, thetaminus, sigma, type);
    ghat = 3*((yplus - yminus)/(2*c_k))*delta;   
%   theta update
    theta=theta-a_k*ghat;
    % Two lines below invoke constraints
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
  end
  % START enhanced 2RDSA-Unif ITERATIONS FOLLOWING INITIALIZATION
%  
for k=1:(numSimulations-N)/3
    a_k=a2/(k+A2)^alpha2;
    c_k=c2/(k)^gamma2;
%     w_k=1/(k+1);
    
%     w_k=w/(k+1)^d;
    w=w+(c_k)^(4);
    w_k=(c_k)^(4)/w;
    
%     a_k=a2/(k+1+A2)^alpha2;
%     c_k=c2/(k+1)^gamma2;
%     w_k=w/(k+1)^d;
    
% GENERATION OF AVERAGED GRADIENT AND HESSIAN (NO AVERAGING IF gH_avg=1)                 
      delta=unifrnd(-1,1,p,1);
      thetaplus=theta+c_k*delta;
      thetaminus=theta-c_k*delta;  
      yplus=feval(loss, p, thetaplus, sigma, type);
      yminus=feval(loss, p, thetaminus, sigma, type);
      ghat=3*(yplus-yminus)/(2*c_k)*delta;  
% GENERATE THE HESSIAN UPDATE
      M_n = delta*delta';
      for idx=1:p
          M_n(idx,idx) = 5/2*(M_n(idx,idx) - 1/3);
      end      
% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS. PER ITERATION      
      y=feval(loss, p, theta, sigma, type);      
      Hhat = 9/2*((yplus+yminus-2*y)/(c_k^2))*M_n;
%   if (k == 1)
%       Psi_k = zeros(p,p);
%   else    
     %Hbardig=Hbar;
     Hbaroffdig=Hbarbar;
     %M_ndig=M_n;
     M_noffdig=M_n;
     for idx=1:p
          M_noffdig(idx,idx) = 0;
          Hbaroffdig(idx,idx) = 0;
     end 
     M_ndig=M_n-M_noffdig;
     Hbardig=Hbarbar-Hbaroffdig;
     Psi_k=(9/2)*(M_ndig*(delta'*Hbaroffdig*delta)+M_noffdig*(delta'*Hbardig*delta));
     Psi_k=.5*Psi_k+.5*Psi_k';
%   end
     Hhatinput=.5*(Hhat+Hhat')-Psi_k;    
   % Hbar=(k/(k+1))*Hbar+Hhat/(k+1); 
     Hbar=(1-w_k)*Hbar+w_k*Hhatinput;         
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
str = sprintf('Normalized loss: %e +- %e, Normalised MSE: %e +- %e', losstheta/replications/Ltheta0, std(lossesAllReplications)/(replications^.5), errtheta/replications/mseTheta0, std(nmseAllReplications)/(replications^.5));
disp(str);
%disp(mat2str(theta,4));  