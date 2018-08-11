%J.C. Spall, Jan. 2000
%This code implements SPSA with constraints for theta to lie in 
%a specified hypercube (i.e., component-wise constraints).  Allows for multiple replications   
%for purposes of statistical evaluation based on knowledge of true (noise-free) loss value
%(set replications=1 if user only wants one run).

% Prashanth L.A., Jul. 2015
% The code below is from J.C. Spall's ISSO book website, except that it
% does not perform any "blocking" test for the parameter.
% 
% Parameters:
% p -> dimension of the problem
% sigma -> noise parameter. Noise is (p+1)-dimensional Gaussian with variance sigma^2
% type -> 1 for quadratic, 2 for fourth-order loss, 3 = Powell singular function, 4 = Rosenbrock function, 5 = Rastrigin function
% numSimulation -> this is the simulation budget that impacts the number of 2SPSA iterations
% replications -> number of independent simulations
% theta_0 -> initial point
%
function [all] = onespsa(p, sigma, type, numSimulations, replications, theta_0)
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
rand('seed',31415297)
randn('seed',3111113)

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
  % inner loop runs 1SPSA for numIterations (each iter=2 function evals)
  for k=0:numSimulations/2-1
    ak = a/(k+1+A)^alpha;
    ck = c/(k+1)^gamma;
    delta = 2*round(rand(p,1))-1;
    thetaplus = theta + ck*delta;
    thetaminus = theta - ck*delta;
    yplus=feval(loss, p, thetaplus, sigma, type);
    yminus=feval(loss, p, thetaminus, sigma, type);
    ghat = (yplus - yminus)./(2*ck*delta);
    theta=theta-ak*ghat;
    % Project theta onto a bounded set, component-wise
    theta=min(theta,theta_hi);
    theta=max(theta,theta_lo);
  end
  lossvalue=feval(lossfinaleval, p, theta, type);
  lossfinal=lossfinal+lossvalue;
  errtheta=errtheta+(theta-thetaStar)'*(theta-thetaStar); 
  lossesAllReplications(1, i) = lossvalue/Ltheta0;
  nmseAllReplications(1, i) = (theta-thetaStar)'*(theta-thetaStar)/mseTheta0;
end
disp(['Number of iterations of outer for loop is : ',num2str(k+1)]);

% Display results: normalized loss and mean square error
str = sprintf('Normalized loss: %5.4f +- %5.4f, Normalised MSE: %5.4f +- %5.4f', lossfinal/replications/Ltheta0, std(lossesAllReplications), errtheta/replications/mseTheta0, std(nmseAllReplications));
disp(str);
str = sprintf('Normalised MSE: %10.9f, Std dev: %10.9f',errtheta/replications/mseTheta0, std(nmseAllReplications));
disp(str);
str1 = sprintf('Std Error: %10.9f',std(nmseAllReplications)/sqrt(replications));
disp(str1);
str2 = sprintf('%5.4f (%5.4f) #%d',errtheta/replications/mseTheta0,std(nmseAllReplications)/sqrt(replications),k+1);
disp(str2);
all = [errtheta/replications/mseTheta0,std(nmseAllReplications)/sqrt(replications),k+1];

disp(mat2str(theta,4));
