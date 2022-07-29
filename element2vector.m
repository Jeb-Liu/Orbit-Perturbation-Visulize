function [r_curr,v_curr,f] = element2vector(sigma_curr)
    mu = 3.986004415E14;
    a = sigma_curr(1);
    e = sigma_curr(2);
    ii = sigma_curr(3);
    W = sigma_curr(4);
    w = sigma_curr(5);
    M = sigma_curr(6);
    
    %求真近点角f和地心距r
    syms symE;
    E=vpasolve(M==symE-e*sin(symE));
    f=2*atan(sqrt((1+e^2)/(1-e^2))*tan(E/2));
    f=double(f);
    r=a*(1-e^2)/(1+e*cos(f));
    r_p=a*(1-e^2);
    h=sqrt(a*mu*(1-e^2));
    
    %坐标系变换
    rp = (h^2/mu) * (1/(e*cos(f))) * (cos(f) * [1;0;0] + sin(f) * [0;1;0]);
    vp = (mu/h) * (-sin(f)*[1;0;0] + (e+cos(f))*[0;1;0]);
    
    R3_W = [cos(W) sin(W) 0;
           -sin(W) cos(W) 0;
              0      0    1];
    R1_i = [1    0      0;
            0  cos(ii) sin(ii);
            0 -sin(ii) cos(ii)];
    
    R3_w = [cos(w) sin(w) 0;
           -sin(w) cos(w) 0;
              0      0    1];
          
    Q = R3_W'*R1_i'*R3_w';

    %轨道坐标系位置
    r_o=[r*cos(f) r*sin(f) 0];

    %赤道惯性坐标系-轨道坐标系 方向余弦阵
    ie=[cos(W)*cos(w)-sin(W)*cos(ii)*sin(w) sin(W)*cos(w)+cos(W)*cos(ii)*sin(w) sin(ii)*sin(w)];
    ip=[-cos(W)*sin(w)-sin(W)*cos(ii)*cos(w) -sin(W)*sin(w)+cos(W)*cos(ii)*cos(w) sin(ii)*cos(w)];
    
    
    %赤道惯性坐标系位置速度
    r_curr=r*cos(f).*ie + r*sin(f).*ip;
    v_curr=-mu/h*sin(f).*ie+mu/h*(e+cos(f)).*ip;
    
    r_curr2=(Q*rp)';
    v_curr2=(Q*vp)';
 
end

