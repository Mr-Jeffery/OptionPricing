function AmericanPut()
% calculate the optimal exercise price for American put option
% arrays start with capitial letter
% variables and parametres start with lowercase
        
    taumax = 30;% max time before expire
    n = 100;% node number
    m = 10;% iteration time
    r = 0.05;% risk-free interest rate
    q = 0.05;% continuous dividends rate
    sigma = 0.25;% volatility of the asset
    s = 130;% exercise price
    
    h = taumax/n;% universial spacing
    Tau = 0:h:taumax;% set up the grid

 
    % b0 is the exercise gain at tau=0
    if (q~=0&&r<q)
       b0 = s*(r/q);
    else
        b0 = s;
    end
    B=@(x)b0-x.*(x>0);
    Bpoint = zeros(m+1,n+1);
    Bpoint(1,:) = B(Tau);% initial guess of the boundary
    Bpoint(1,1) = b0;
    % plot(Tau,Bpoint(1,:));% optional, used for ploting the initial guess grid

    %start interation
    for j=2:m+1
        Bpoint(j,:) = s.*exp(Tau.*(q-r)).*N( Tau, Bpoint(j-1,:))./D(Tau,Bpoint(j-1,:)) ;
        plot(Tau,Bpoint(j,:))% optional, used for ploting the grid obtained in current iteration
    end
    

    %functions for iteration
    %all of them can take a whole array as input
    function N = N(Tau,Bt)        
        N = zeros(1,n+1);
        N(1) = 1;
        for i=2:n+1
            tau = Tau(i);
            bt = Bt(i);
            func = exp(r.*Tau(1,1:i)).*normcdf(dminus(tau-Tau(1,1:i),bt./Bt(1,1:i)));
            N(i) = normcdf(dminus(tau,Bt(i)/s)) + r*trapz(Tau(1:i),func);
        end

        function dminus = dminus(tau,z)
        dminus = zeros(1,length(tau));
            for k=1:length(tau)
                if(tau(k)~=0)
                dminus(k) = (log(z(k))+(r-q-0.5*sigma^2).*tau(k))./(sigma.*sqrt(tau(k)));
                end
            end
        end
    end
    
    
    function D = D(Tau,Bt)
        D = zeros(1,n+1);
        D(1) = 1;
        for i=2:n+1
            tau = Tau(i);
            bt = Bt(i);
            func = exp(q.*Tau(1,1:i)).*normcdf(dplus(tau-Tau(1,1:i),bt./Bt(1,1:i)));
            D(i) = normcdf(dplus(tau,Bt(i)/s)) + q*trapz(Tau(1:i),func);
        end
                
        function dplus = dplus(tau,z)
            dplus = zeros(1,length(tau));
            for k=1:length(tau)
                if(tau(k)~=0)
                    dplus(k) = (log(z(k))+(r-q+0.5*sigma^2).*tau(k))./(sigma.*sqrt(tau(k)));
                end
            end
        end

    end
    
    

end