function [T,DT] = Tnm(n,m,sinlambda,coslambda,r,ae,C,S)
    
    %==========================C-1-7===============================%
    sin_m(0+1) = 0;
    sin_m(1+1) = sinlambda;
    for count=2:m
        sin_m(count+1) = 2*coslambda*sin_m(count+1 -1)-sin_m(count+1 -2);
    end

    cos_m(0+1) = 1;
    cos_m(1+1) = coslambda;
    for count=2:m
        cos_m(count+1) = 2*coslambda*cos_m(count+1 -1)-cos_m(count+1 -2);
    end
    
    %===============================C-1-13===============================%
    %%从0开始所以+1
    T = zeros(n+1);
    DT = zeros(n+1);

    for nn=2:n
        for mm=1:nn
            T(nn+1,mm+1)=(ae/r)^nn*( C(nn+1,mm+1)*cos_m(mm)+S(nn+1,mm+1)*sin_m(mm) );
            DT(nn+1,mm+1)=mm*(ae/r)^nn*( C(nn+1,mm+1)*cos_m(mm)-S(nn+1,mm+1)*sin_m(mm) );
        end
    end
    
end