function Solution1 = FBP1Dsolver(dh,alpha,K,b,v,n,r,tau,sigmaP2,k,Slater,Wlater)
% Solver for freeeeee boundary problem:
% u'' + P(x)*u' + Q(x)*u = R(x)     (alpha<x<beta)
% u(alpha)  = a ,alpha is the lower boundary
% u(beta)   = 0 ,beta is the upper boundary
% u'(beta)  = v
% dh is the difference of b, used to calculate quasi-derivative

derivative1 = 0;
while(abs(v-derivative1) > dh)

    Solution0 = ODEsolver(alpha-dh,K,(K-alpha+dh),b,n,r,sigmaP2,k,Slater,Wlater);
    Nodes0 = Solution0(1:2,1:2);
    derivative0 = (Nodes0(2,2)-Nodes0(1,2))/(Nodes0(2,1)-Nodes0(1,1));
    
    Solution1 = ODEsolver(alpha,K,(K-alpha),b,n,r,sigmaP2,k,Slater,Wlater);
    Nodes1 = Solution1(1:2,1:2);
    derivative1 = (Nodes1(2,2)-Nodes1(1,2))/(Nodes1(2,1)-Nodes1(1,1));
    
    Solution2 = ODEsolver(alpha+dh,K,(K-alpha-dh),b,n,r,sigmaP2,k,Slater,Wlater);
    Nodes2 = Solution2(1:2,1:2);
    derivative2 = (Nodes2(2,2)-Nodes2(1,2))/(Nodes2(2,1)-Nodes2(1,1));
    
    Db = (derivative2 - derivative0)/(2*dh); %find d(u'(beta))/d(beta)
    alpha1 = alpha + (v - derivative1)/(Db);

    if(0<alpha1&&alpha1<K)
        alpha = alpha1;
    elseif(alpha1>K)
        alpha = (alpha+K)/2;
    elseif(alpha1<0)
        alpha = alpha/2;
    end
    dh = min(dh,1/Nodes2(2,2));
    plot(Solution1(:,1),Solution1(:,2))
end

