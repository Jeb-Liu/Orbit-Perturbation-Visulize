function [Fu] = C_1_10(x, y, z)
    %
    % gravity for non-sphere earth
    %
    global C S
    mu = 3.986004415e14;
    ae = 6.378136300e6;
    
    r = sqrt(x^2 + y^2 + z^2);
    sinphi = z/r;
    cosphi = sqrt(x^2 + y^2)/r;
    sinlambda = y/sqrt(x^2 + y^2);
    coslambda = x/sqrt(x^2 + y^2);
    
    n=12;m=n;
    
    [P_n,DP_n] = Pn(n,m,sinphi,cosphi);
    [P_nm,DP_nm] = Pnm(n,m,sinphi,cosphi);
    [T_nm,DT_nm] = Tnm(n,m,sinlambda,coslambda,r,ae,C,S);
    
    % ==========================C-1-11========================== %
    drdrr = [x/r , y/r , z/r].';
    dphidrr = [-sinphi*coslambda/r, -sinphi*sinlambda/r, cosphi/r].';
    dlambdadrr = [-y/(x^2 + y^2), x/(x^2 + y^2), 0].';
    
    % ==========================C-1-12========================== %
    ans1=0;ans2=0;ans3=0;ans4=0;ans5=0;
    % C-1-12 dvdr %
    for nn=2:n
        ans1 = ans1 + (nn+1)*(ae/r)^nn*C(nn+1,0+1)*P_n(nn+1);
    end
    for nn=2:n
        for mm=1:nn
            ans2 = ans2 + (nn+1)*P_nm(nn,mm)*T_nm(nn+1,mm+1);
        end
    end
    dvdr = -mu/r^2 * ( ans1+ans2 );
    
    % C-1-12 dvdphi %
    for nn=2:n
        ans3 = ans3 + (nn+1)*(ae/r)^nn*C(nn+1,0+1)*DP_n(nn+1);
    end
    for nn=2:n
        for mm=1:nn
            ans4 = ans4 + DP_nm(nn,mm)*T_nm(nn+1,mm+1);
        end
    end
    dvdphi = mu/r * ( ans3+ans4 );
    
    % C-1-12 dvdlambda %
    for nn=2:n
        for mm=1:nn
            ans5 = ans5 + P_nm(nn,mm)*DT_nm(nn+1,mm+1);
        end
    end
    dvdlambda = mu/r * ans5;
    
    Fu = dvdr*drdrr + dvdphi*dphidrr + dvdlambda*dlambdadrr;
    Fu(1) = +Fu(1) -mu/r^3*x;
    Fu(2) = +Fu(2) -mu/r^3*y;
    Fu(3) = +Fu(3) -mu/r^3*z;
end

