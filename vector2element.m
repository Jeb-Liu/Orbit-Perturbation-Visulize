% clear;clc;
% 
% sigma_start = [6493.491894^2 7.7253696509E-004 0.0958123859 -3.0051431274E+000 -1.002638139 3.3319386917E-001];
% % syms symE;
% %     E=vpasolve(sigma_start(6)==symE-sigma_start(2)*sin(symE));
% %     f=2*atan(sqrt((1+sigma_start(2)^2)/(1-sigma_start(2)^2))*tan(E/2));
% %     f=double(f);
% site_out = [6493.491894^2/1000 7.7253696509E-004 0.0958123859 -3.0051431274E+000 -1.002638139 3.3319386917E-001];
% [r1,v1,f1] = element2vector_v0(sigma_start);
% [r2,v2,f2] = element2vector(sigma_start);
% ele1 = vector2element1(r1,v1);
% ele2 = vector2element1(r2,v2);
% error1 = ele1 - sigma_start;
% error2 = ele2 - sigma_start;

function [sigma_curr,f] = vector2element(r_curr,v_curr)
    mu = 3.986004415E14;
    r = norm(r_curr);
    v = norm(v_curr);
    vr = dot(r_curr,v_curr)/r;
    hh = cross(r_curr,v_curr);
    h = norm(hh);
    %倾角i
    ii = acos(hh(3)/h);
    nn = cross([0 0 1],hh);
    n = norm(nn);
    %RAAN
    if n~=0
    W = acos(nn(1)/n);
    if nn(2)<0
    W = 2*pi - W;
    end
    else
    W=0;
    end
    
    ee = 1/mu*((v^2-mu/r)*r_curr-r*vr*v_curr);
    e = norm(ee);
    
    %近地点幅角w
    if n~=0
        if e>1E-10
            w = acos(dot(nn,ee)/n/e);
            if ee(3)<0
        w = 2*pi - w;
            end
        else
            w=0;
        end
    else
        w=0;
    end
    
%     w = atan(ee(3)/( (ee(1)*cos(W)+ee(2)*sin(W)*sin(ii)) ));
    
    %真近点角
    if e>1E-10
        f = acos(dot(ee,r_curr)/e/r);
        if vr<0
            f = 2*pi - f;
        end
    else
        cp =cross(nn,r_curr);
        if cp(3)>=0
            f = acos(dot(nn,r_curr)/e/r);
        else
            f = 2*pi - acos(dot(nn,r_curr)/e/r);
        end
    end
    
    %半长轴a
    a = h^2/(mu*(1-e^2));
    
    %平近点角
    E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));
    M = E-e*sin(E);
    
    sigma_curr = [a e ii W w M];
end

