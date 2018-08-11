function y = loss_myexample(p, theta, type)
% type 1 = quadratic, 2 = fourth order function from Spall's 2000 paper,
% 3 = Powell singular function, 4 = Rosenbrock function, 5 = Rastrigin function

% Create upper triangular matrix with entries 1/p
tmpMat = zeros(1,p);
tmpMat(1,1) = 1;
Amat=1/p*toeplitz(tmpMat,ones(1,p));

% based on type, return the appropriate loss - quadratic if 1, 4th order
% otherwise
if type == 1
    b = ones(p, 1);
    y=theta'*Amat*theta + theta'* b;
elseif type == 2 % Spall 2000
    y=theta'*(Amat'*Amat)*theta + 0.1* sum((Amat*theta).^3) + 0.01 * sum((Amat*theta).^4);
elseif type == 3 % Powell singular function (n = 10)
%     h1(x,ξx) = −1 − (xi−1 +10xi)2 +5(xi+1 −xi+2)2 +(xi −2xi+1)4 +10(xi−1 −xi+2)4 +ξx, 
    temp = 0;
    for i = 2:p-2
       temp = temp + (theta(i-1) + 10.*theta(i)).^2 + 5.*(theta(i+1) - theta(i+2)).^2 + (theta(i) - 2.*theta(i+1)).^4 + 10.*(theta(i-1) - theta(i+2)).^4;
    end
    powell = -1 - temp;
    y = -powell;
 elseif type == 4 % Rosenbrock function (n = 10)
%     h2(x,ξx)=−1−  100(xi+1 −x2i)2 +(xi −1)2 +ξx,
    temp = 0;
    for i = 1:p-1
       temp = temp + 100.*(theta(i+1) - theta(i).^2).^2 + (theta(i) - 1).^2;
    end
    powell = -1 - temp;
    y = -powell;
  elseif type == 5 % Rastrigin function (n = 10)
%     h2(x,ξx)=h4(x,ξx)=−  x2i −10cos(2πxi) −10n−1+ξx,
    temp = 0;
    for i = 1:p
       temp = temp + theta(i).^2 - 10.*cos(2*pi*theta(i));
    end
    powell = -temp -10*p -1;
    y = -powell;
end
