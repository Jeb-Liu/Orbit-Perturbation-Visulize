function [P_n,DP_n] = Pn(n, m, sinphi, cosphi)

    P_n(0+1) = 1;
    P_n(1+1) = sqrt(3)*sinphi;
    if n>=2
        for count=2:n
            P_n(count+1) = sqrt(4*n^2-1)/n*sinphi*P_n(count+1 -1) - (n-1)/n*sqrt((2*n+1)/(2*n-3))*P_n(count+1 -2);
        end
    end
    
    DP_n(0+1) = 0;
    DP_n(1+1) = sqrt(3)*cosphi;
    if n>=2
        for count=2:n
            DP_n(count+1) = sqrt(4*n^2-1)/n*sinphi*DP_n(count+1 -1) ...
                           + sqrt(4*n^2-1)/n*cosphi*P_n(count+1 -1) ...
                           - (n-1)/n*sqrt((2*n+1)/(2*n-3))*DP_n(count+1 -2);
        end
    end
    
end

