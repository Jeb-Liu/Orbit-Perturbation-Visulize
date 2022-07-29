function [P_nm,DP_nm] = Pnm(n, m, sinphi, cosphi)
    %===============================C-1-9================================%
    P_nm = zeros(n);
    P_nm(1,1) = sqrt(3) * cosphi;
    % 对角线 n>=2 %
    for count = 2:n
        P_nm(count, count) = sqrt((2*n+1)/(2*n)) * cosphi * P_nm(count-1, count-1);
    end
    % 对角线下一格 n>=1 %
    for count = 1:n-1
        P_nm(count+1, count) = sqrt(2*n+3) * sinphi * P_nm(count, count);
    end
    % 列 %
    for mm = 1:n-1
        for nn = mm+2:n
            P_nm(nn,mm) = sqrt((4*nn^2-1)/(nn^2-mm^2)) * sinphi * P_nm(nn-1,mm) ...
                         - sqrt( (2*nn+1)*((nn-1)^2-mm^2)/(2*nn-3)/(nn^2-mm^2) )*P_nm(nn-2,mm);
        end
    end
    
    %===============================C-1-15===============================%
    DP_nm = zeros(n);
    DP_nm(1,1) = -sqrt(3) * sinphi;
    % 对角线 n>=2 %
    for count = 2:n
        DP_nm(count, count) = sqrt((2*n+1)/2/n) * ( cosphi*DP_nm(count-1, count-1) - sinphi*P_nm(count-1, count-1) );
    end
    % 对角线下一格 n>=1 %
    for count = 1:n-1
        DP_nm(count+1, count) = sqrt(2*n+3) * ( cosphi*P_nm(count, count) + sinphi*DP_nm(count,count) );
    end
    % 列 %
    for mm = 1:n-1
        for nn = mm+2:n
            DP_nm(nn,mm) = sqrt((4*nn^2-1)/(nn^2-mm^2)) * ( cosphi*P_nm(nn-1,mm)+sinphi*DP_nm(nn-1,mm) ) ...
                         - sqrt( (2*nn+1)*((nn-1)^2-mm^2)/(2*nn-3)/(nn^2-mm^2) )*DP_nm(nn-2,mm);
        end
    end
    
end

