function thetaStar = getOptima(p, type)
% type 1 = quadratic, 2 = fourth order function from Spall's 2000 paper
% 3 = Powell singular function, 4 = Rosenbrock function, 5 = Rastrigin function

% Create upper triangular matrix with entries 1/p
tmpMat = zeros(1,p);
tmpMat(1,1) = 1;
Amat=1/p*toeplitz(tmpMat,ones(1,p));
% based on type, return the appropriate loss - quadratic if 1, 4th order
% otherwise
if type == 1
    b = ones(p, 1);    
    thetaStar = -inv(Amat + Amat')*b;
elseif type == 2 
    thetaStar = zeros(p,1);
elseif type == 3 
    thetaStar = zeros(p,1);
elseif type == 4 
    thetaStar = ones(p,1);
elseif type == 5 
    thetaStar = zeros(p,1);
end

