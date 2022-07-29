clear;clc;

%北斗时T0=2006年1月1日0点0分0秒
%北斗时T0 JD=2453736.5
%start 北斗时 Week 782 T:432000.0000
%start JD=2453736.5+782*7*86400+432000
%end   北斗时 Week 782 T:0.0000
%GE0:3074.66277 42164000

%=====历书数据=====
alc_data=csvread('alc.csv');   %同步轨道
%alc_data=csvread('alc11.csv');%高轨道
alc_data(:,9) = alc_data(:,1)*7*86400 + alc_data(:,2) - (782*86400*7+432000);
alc_length = length(alc_data);
%=====EGM2008=====
global C S
C = csvread('EGM2008_C.csv');
S = csvread('EGM2008_S.csv');
%=====H_P_model=====
global H_P
H_P = csvread('H-P-atomsphere-model.csv');
%===============================变量=========================================
global JD_T ecliptic_fix;
mu = 3.986004415E14;
JD_T = 2453736.5+782*7*86400+432000;

%=====地轴旋转矩阵ecliptic to ECI=====
earth_axial_tilt=deg2rad(-23.44);
ecliptic_fix=[1 0 0;0 cos(earth_axial_tilt) -sin(earth_axial_tilt);0 sin(earth_axial_tilt) cos(earth_axial_tilt)];
%=====ECI to ECEF=====
%ECI2ECEF=[cos(2*pi/86184*T) -sin(2*pi/86184*T) 0;sin(2*pi/86184*T) cos(2*pi/86184*T) 0;0 0 1];
     
%=================================初始状态=====================================
%=====计算参数=====
step = 60;
drawtimes = 60;
day = 1;
t_end = 86400*day;
%=====卫星初始=====
global A_M C_D
A_M=200/50;%17.4/50
C_D=2.0;
%=====初始轨道根数=====
%sigma_start=[6508140 1e-5 deg2rad(77.6) deg2rad(45) deg2rad(45) 0];
sigma_start=[8058997.305 9.019997558e-4 1.937007381 deg2rad(122) pi/2 0];%回归 太阳同步 冻结
%sigma_start=[7878140 0.01 deg2rad(30) 0 0 0];
%sigma_start=[26562560,0.01,deg2rad(89.9),-pi/2,-pi/2,0.333];
%历书数据
%sigma_start=[alc_data(1,3)^2 alc_data(1,4) alc_data(1,5) alc_data(1,6) alc_data(1,7) alc_data(1,8)];
sigma_end=[alc_data(alc_length,3)^2 alc_data(alc_length,4) alc_data(alc_length,5) alc_data(alc_length,6) alc_data(alc_length,7) alc_data(alc_length,8)];

[end_r,end_v,end_f]=element2vector(sigma_end);
[orbit_r,orbit_v,orbit_f]=element2vector(sigma_start);
% orbit_i = deg2rad(30);
% norm_orbit_v = 3074.66277;
% orbit_r = [42164000 1E-8 1E-8];
% orbit_v = [1E-8 norm_orbit_v*cos(orbit_i) norm_orbit_v*sin(orbit_i)];

%=====月球初始位置=====
sigma_moon_start=moonJ2000(JD_T);
[j2000ecl_moon_r,j2000ecl_moon_v,j2000ecl_moon_f]=element2vector(sigma_moon_start);
moon_r = (ecliptic_fix.' * j2000ecl_moon_r.').';
moon_v = (ecliptic_fix.' * j2000ecl_moon_v.').';
% moon_i = deg2rad(5.145);
% norm_moon_v = 1075.913124;
% moon_r = [363228900 1E-8 1E-8];
% moon_v = [1E-8 norm_moon_v*cos(moon_i) norm_moon_v*sin(moon_i)];

%==============================计算================================================
%=====微分方程=====
options=odeset('RelTol',1E-13, 'AbsTol', [1E-15,1E-15,1E-15,1E-15,1E-15,1E-15,1E-15,1E-15,1E-15,1E-15,1E-15,1E-15]);
[T,R]=ode45('EarthMoonSatellite',0:step:t_end, [moon_r,moon_v,orbit_r,orbit_v],options);
%[T,R]=ode45('threeBody',0:step:t_end, [moon_r,moon_v,orbit_r,orbit_v],options);

%=====卫星位置速度=====
curr_rr=[R(:,7),R(:,8),R(:,9)];
curr_vv=[R(:,10),R(:,11),R(:,12)];
curr_r=sum(abs(curr_rr).^2,2).^(1/2);
curr_v=sum(abs(curr_vv).^2,2).^(1/2);
fix_rr=(ecliptic_fix*curr_rr(1:length(T),:).').';
ecef_rr = zeros(max(size(T)),3);
for count=1:length(T)
    ecef_rr(count,:)=(ecliptic_fix*[cos(2*pi/86184*T(count)) -sin(2*pi/86184*T(count)) 0;sin(2*pi/86184*T(count)) cos(2*pi/86184*T(count)) 0;0 0 1].' ...
             *curr_rr(count,:).').';
end

%=====月球位置速度=====
moon_rr=[R(:,1),R(:,2),R(:,3)];
moon_vv=[R(:,4),R(:,5),R(:,6)];
fix_moon_rr=(ecliptic_fix*moon_rr(1:length(T),:).').';

%=====误差=====
error_r = norm( curr_rr(length(curr_rr),:) )- norm(end_r);
error_v = norm( curr_vv(length(curr_vv),:) )- norm(end_v);

%=====卫星轨道根数=====
sigma_t = zeros(max(size(T)),6);
for count1=1:length(curr_rr)
    sigma_t(count1,:) = vector2element(curr_rr(count1,:),curr_vv(count1,:));
end

%=====卫星轨道根数变化=====
d_sigma_t=zeros(max(size(T)),6);
d_sigma_t(1)=0;
for d=1:length(T)-1
    d_sigma_t(d+1,:)=sigma_t(d+1,:)-sigma_t(d,:);
    while abs(d_sigma_t(d+1,5))>(pi-1)
        if d_sigma_t(d+1,5)>0
            d_sigma_t(d+1,5) = d_sigma_t(d+1,5) - pi;
        else
            d_sigma_t(d+1,5) = d_sigma_t(d+1,5) + pi;
        end
    end
    
end

for d=1:length(T)
    if d_sigma_t(d:4)>pi
        d_sigma_t(d:4)=0;
    end
end


%===============================绘图=======================================
%figure(1)
figure(1);
%=====sigma=====
subplot(5,2,1);
plot(T/86400,sigma_t(:,1));
xlabel('day');title('a [m]');hold on;

subplot(5,2,3);
plot(T/86400,sigma_t(:,2));
title('e [-]');xlabel('day');hold on;

subplot(5,2,5);
plot(T/86400,sigma_t(:,3));
title('i [deg]');xlabel('day');hold on;

subplot(5,2,7);
plot(T/86400,sigma_t(:,4));
title('{\Omega}[deg]');xlabel('day');hold on;

subplot(5,2,9);
plot(T/86400,sigma_t(:,5));
title('w [s]');xlabel('day');hold on;

%=====d_sigma=====
subplot(5,2,2);
plot(T/86400,d_sigma_t(:,1));
xlabel('day');title('{\Delta}a [m]');

subplot(5,2,4);
plot(T/86400,d_sigma_t(:,2));
title('{\Delta}e [-]');xlabel('day');

subplot(5,2,6);
plot(T/86400,d_sigma_t(:,3));
title('{\Delta}i [deg]');xlabel('day');

subplot(5,2,8);
plot(T/86400,d_sigma_t(:,4));
title('{\Delta}{\Omega}[deg]');xlabel('day');

subplot(5,2,10);
plot(T/86400,d_sigma_t(:,5));
title('{\Delta}w [s]');xlabel('day');

%figure(2)
figure(2);

%=======坐标轴J2000地心黄道=====
xaxis=[0 0 0;10000000 0 0];
yaxis=[0 0 0;0 10000000 0];
zaxis=[0 0 0;0 0 10000000];
axiss=[0,0,0;10000000,0,0;0,10000000,0;0,0,10000000];
scatter3(axiss(:,1), axiss(:,2), axiss(:,3),'filled');
hold on;
plot3(xaxis(:,1), xaxis(:,2), xaxis(:,3));
hold on;
plot3(yaxis(:,1), yaxis(:,2), yaxis(:,3));
hold on;
plot3(zaxis(:,1), zaxis(:,2), zaxis(:,3));
hold on;

% %=======黄道面=====
% [s_x,s_y]=meshgrid(-400000000:40000000:400000000);
% s_z=s_x-s_x;
% s=surf(s_x,s_y,s_z,'FaceAlpha',0.1);
% s.EdgeColor = 'none';

xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');
box on;
grid on;
axis([-1000000,1000000,-1000000,1000000,-1000000,1000000]);
axis equal;
pause(1);
hold on;

%画倾斜的地球
[earth_x,earth_y,earth_z]=sphere;
%surf(6371000*earth_x, 6371000*earth_y, 6371000*earth_z);
coor = [earth_x(:),earth_y(:),earth_z(:)];
newcoor = (ecliptic_fix*coor.').';
new_earth_x = reshape(newcoor(:,1), size(earth_x));
new_earth_y = reshape(newcoor(:,2), size(earth_y));
new_earth_z = reshape(newcoor(:,3), size(earth_z));
surf(6371000*new_earth_x, 6371000*new_earth_y, 6371000*new_earth_z);

line1=animatedline;
line2=animatedline;
%=====一次性plot轨道=====%
plot3(fix_rr(:,1), fix_rr(:,2), fix_rr(:,3));%惯性系
%plot3(fix_moon_rr(:,1), fix_moon_rr(:,2), fix_moon_rr(:,3));
%plot3(ecef_rr(:,1), ecef_rr(:,2), ecef_rr(:,3));%地心地固

%=====轨道动图=====%
 for count2=1:step*drawtimes:length(T)
    
     %=====ECI ECEF MOON=====%
     %addpoints(line1,fix_rr(count2,1),fix_rr(count2,2),fix_rr(count2,3));
     %addpoints(line1,ecef_rr(count2,1),ecef_rr(count2,2),ecef_rr(count2,3));
     %addpoints(line2,fix_moon_rr(count2,1),fix_moon_rr(count2,2),fix_moon_rr(count2,3));
     drawnow;
     
     
%      %月球
%      [moon_x,moon_y,moon_z]=sphere;
%      moon_move=surf(1737100*moon_x + fix_moon_rr(count2,1) ...
%                   , 1737100*moon_y + fix_moon_rr(count2,2) ...
%                   , 1737100*moon_z + fix_moon_rr(count2,3));
     
 end

%=====PLOT alc数据=====%
% figure(1);
% subplot(5,2,1)
% scatter(alc_data(:,9)/86400,alc_data(:,3).^2);
% subplot(5,2,3)
% scatter(alc_data(:,9)/86400,alc_data(:,4));
% subplot(5,2,5)
% scatter(alc_data(:,9)/86400,alc_data(:,5));
% subplot(5,2,7)
% scatter(alc_data(:,9)/86400,alc_data(:,6)+2*pi);
% subplot(5,2,9)
% scatter(alc_data(:,9)/86400,alc_data(:,7)+pi);

%=====星下点=====%


% ecif_rr 0~2pi
ecif_rr = zeros(max(size(T)),3);
for count=1:length(T)
    ecif_rr(count,:)=([cos(2*pi/86184*T(count)) -sin(2*pi/86184*T(count)) 0;sin(2*pi/86184*T(count)) cos(2*pi/86184*T(count)) 0;0 0 1].' ...
             *curr_rr(count,:).').';
end
phi=zeros(max(size(T)),1);
lambda=zeros(max(size(T)),1);
for count3=1:length(T)
    if ecif_rr(count3,2)>=0
        phi(count3) = asin(ecif_rr(count3,3)/sqrt(ecif_rr(count3,1)^2+ecif_rr(count3,2)^2+ecif_rr(count3,3)^2));
        lambda(count3) = acos(ecif_rr(count3,1)/sqrt(ecif_rr(count3,1)^2+ecif_rr(count3,2)^2));
    else
        phi(count3) = asin(ecif_rr(count3,3)/sqrt(ecif_rr(count3,1)^2+ecif_rr(count3,2)^2+ecif_rr(count3,3)^2));
        lambda(count3) = 2*pi - acos(ecif_rr(count3,1)/sqrt(ecif_rr(count3,1)^2+ecif_rr(count3,2)^2));
    end
    
    %0~2pi to -pi~pi
    if lambda(count3)>=pi
        lambda(count3)=lambda(count3)-2*pi;
    end
    
end

%分段
count5=1;
loop(1)=0;
    for count4=1:length(T)
        if count4==length(T)
        else
            if abs(lambda(count4+1)-lambda(count4))>pi
            loop(1+count5)=count4;
            count5=count5+1;
            end
        end
    end
    
%
figure(3);
axis([-180,180,-90,90]);
grid on;
geoshow('landareas.shp','FaceColor','green');
hold on;
set(gca,'DataAspectRatio',[1 1 1]);
if norm(loop)==0
    plot(rad2deg(lambda(:)),rad2deg(phi(:)));
else
    for loops=1:length(loop)-1
    plot(rad2deg(lambda(loop(loops)+1:loop(loops+1))),rad2deg(phi(loop(loops)+1:loop(loops+1))),'b-');
    hold on;
    end
    plot(rad2deg(lambda(loop(loops+1)+2:length(T))),rad2deg(phi(loop(loops+1)+2:length(T))),'b-');
end

figure(4);
plot(T/86400,(sqrt(fix_rr(:,1).^2+fix_rr(:,2).^2+fix_rr(:,3).^2)-6371000)/1000);
title('height');xlabel('day');ylabel('km');


