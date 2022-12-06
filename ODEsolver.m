function output = ODEsolver(alpha,beta,a,b,n,r,sigmaP2,k,Slater,Wlater)
% solver for boundary problem 
% u'' + P(x)*u' + Q(x)*u = R(x)     (alpha<x<beta)
% r1*u'(alpha)  + r2*u(alpha)   = r3
% r4*u'(beta)   + r5*u(beta)    = r6
% r can be zero

%plug in parametres
h = (beta-alpha)/n;
M = zeros(n+1);
D = zeros(n+1,1);
X = zeros(n+1,1);

x = alpha;
%for alpha/x1
X(1) = alpha;
M(1,1) = 1;
D(1) = a;

%for x2 to xn
for i = 2:1:n
    x = x + h;
    X(i) = x;
    M(i,i-1) = k*sigmaP2/2 + r*x*h*k/2;
    M(i,i) = -k*sigmaP2 - k*r*h^2 - h^2;
    M(i,i+1) = k*sigmaP2/2 - r*x*h*k/2;
end
D(2:n) = -h^2.*laterValue(Slater,Wlater,X(2:n));

%beta/xn+1
X(n+1) = beta;
M(n+1,n+1) = 1;
D(n+1) = b;

%solve linear equations
U = M\D;
output = [X U];
end
