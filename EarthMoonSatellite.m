function dr=EarthMoonSatellite(t,r_r)
%=====太阳位置J2000赤道�?=====
global JD_T ecliptic_fix;
AU=1.49597870700e11;
L = deg2rad(280.46645 + 36000.76983*(JD_T+t) + 0.0003032*(JD_T+t)^2);
l = deg2rad(357.52910 + 35999.05030*(JD_T+t) - 0.0001559*(JD_T+t)^2 - 0.00000048*(JD_T+t)^3);
dO = ( deg2rad(1.914600 - 0.004817*(JD_T+t) - 0.000014*(JD_T+t)^2) )*sin(l) ...
    + ( deg2rad(0.019993 - 0.000101*(JD_T+t)) )*sin(2*l) + deg2rad(0.000290)*sin(3*l);
O = L + dO;
f = l + dO;
e = 0.016708617 - 0.000042037*(JD_T+t) - 0.0000001236*(JD_T+t)^2;
R = 1.000001018*(1-e^2)/(1+e*cos(f))*AU;
O_2000 = O - deg2rad(0.01397)*(2020-2000);
Sun_rr_ecl = [R*cos(O_2000) R*sin(O_2000) 0];
Sun_rr_Hel = (ecliptic_fix.' * Sun_rr_ecl.').';

%=====各距�?=====
global A_M C_D
%r_r(1)==moonX
%r_r(2)==moonY
%r_r(3)==moonZ
%r_r(7)==X
%r_r(8)==Y
%r_r(9)==Z
rSun2Satellite=[r_r(7)-Sun_rr_Hel(1),r_r(8)-Sun_rr_Hel(2),r_r(9)-Sun_rr_Hel(3)];
rMoon2Satellite=[r_r(7)-r_r(1),r_r(8)-r_r(2),r_r(9)-r_r(3)];
rSun2Moon=[r_r(1)-Sun_rr_Hel(1),r_r(2)-Sun_rr_Hel(2),r_r(3)-Sun_rr_Hel(3)];

absrSun2Moon=norm(rSun2Moon);
absrSun2Satellite=norm(rSun2Satellite);
absrMoon2Earth=sqrt(r_r(1)^2 + r_r(2)^2 + r_r(3)^2);
absrSatellite2Earth=sqrt(r_r(7)^2 + r_r(8)^2 + r_r(9)^2);
absrMoon2Satellite=norm(rMoon2Satellite);
dr=zeros(12,1);

muSun = 1.32712440041E20;
muEarth = 3.986004415E14;
muMoon = 4.902785E12;
ae = 6.378136300e6;

J2 = 1.0826355e-3;
J3 = -2.5324105e-6;
J4 = -1.6198976e-6;
AJ2 = 3/2*J2*(ae/absrSatellite2Earth)^2;
AJ3 = 5/2*J3*(ae/absrSatellite2Earth)^3;
AJ4 = 5/8*J4*(ae/absrSatellite2Earth)^4;

%Fu = C_1_10( r_r(7), r_r(8), r_r(9) );

rho=H_P_atomsphere(r_r,absrSatellite2Earth,Sun_rr_Hel);
% rho=0;


dr(1) = r_r(4);
dr(4) = - muEarth/absrMoon2Earth^3*r_r(1) ...
        - muSun/absrSun2Moon^3 * rSun2Moon(1) ;

dr(2) = r_r(5);
dr(5) = - muEarth/absrMoon2Earth^3*r_r(2) ...
        - muSun/absrSun2Moon^3 * rSun2Moon(2) ;

dr(3) = r_r(6);
dr(6) = - muEarth/absrMoon2Earth^3*r_r(3) ...
        - muSun/absrSun2Moon^3 * rSun2Moon(3) ;
%{======================================================
dr(7) = r_r(10);
dr(10) = -muEarth/absrSatellite2Earth^3*r_r(7)...
    *( 1+AJ2*(1-5*r_r(9)^2/absrSatellite2Earth^2)...
    +AJ3*(3*r_r(9)/absrSatellite2Earth-7*r_r(9)^3/absrSatellite2Earth^3)...
    -AJ4*(3-42*r_r(9)^2/absrSatellite2Earth^2+63*r_r(9)^4/absrSatellite2Earth^4) )...
         - 1/2*rho*A_M*C_D*r_r(10)*r_r(10) ;
%              - muMoon/absrMoon2Satellite^3 * rMoon2Satellite(1) ...
%          - muSun/absrSun2Satellite^3 * rSun2Satellite(1) ;


         

dr(8) = r_r(11);
dr(11) = -muEarth/absrSatellite2Earth^3*r_r(8)...
    *( 1+AJ2*(1-5*r_r(9)^2/absrSatellite2Earth^2)...
    +AJ3*(3*r_r(9)/absrSatellite2Earth-7*r_r(9)^3/absrSatellite2Earth^3)...
    -AJ4*(3-42*r_r(9)^2/absrSatellite2Earth^2+63*r_r(9)^4/absrSatellite2Earth^4) )...
        - 1/2*rho*A_M*C_D*r_r(11)*r_r(11)  ;
%              - muMoon/absrMoon2Satellite^3 * rMoon2Satellite(2) ...
%          - muSun/absrSun2Satellite^3 * rSun2Satellite(2) ;


         

dr(9) = r_r(12);
dr(12) = -muEarth/absrSatellite2Earth^3*r_r(9)...
    *( 1+AJ2*(3-5*r_r(9)^2/absrSatellite2Earth^2)...
    +AJ3*(6*r_r(9)/absrSatellite2Earth-7*r_r(9)^3/absrSatellite2Earth^3-3/5*(absrSatellite2Earth/r_r(9)) )...
    -AJ4*(15-70*r_r(9)^2/absrSatellite2Earth^2+63*r_r(9)^4/absrSatellite2Earth^4) )...
        - 1/2*rho*A_M*C_D*r_r(12)*r_r(12)  ;
%              - muMoon/absrMoon2Satellite^3 * rMoon2Satellite(3) ...
%          - muSun/absrSun2Satellite^3 * rSun2Satellite(3) ;


         

end