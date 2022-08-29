% Margin Analysis of exp(x)
% u'' + P(x)*u' + Q(x)*u = R(x)     (alpha<x<beta)
% r1*u'(alpha)  + r2*u(alpha)   = r3
% r4*u'(beta)   + r5*u(beta)    = r6
% r can be zero

%initiate
close all;
clear;
clc;
clf;

MaxMargin = zeros(3,2);
for i = 1:1:4
    n = 10^i;
    Solution = ODEsolver(0,1,0,1,1,1,0,exp(1),n,@(x)(1),@(x)(1),@(x)(3*exp(x)));
    Margin = exp(Solution(:,1)) - Solution(:,2);
    MaxMargin(i,1) = -i;
    MaxMargin(i,2) = log(max(Margin));
end
plot(MaxMargin(:,1),MaxMargin(:,2))
axis([-5,0,-5,0])