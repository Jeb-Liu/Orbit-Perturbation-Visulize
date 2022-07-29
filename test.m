clear;clc;

global C S
C = readmatrix('EGM2008_C.csv');
S = readmatrix('EGM2008_S.csv');

x = 1;
y = 6371000;
z = 0;
  
Fu = C_1_10(x, y, z);

g = norm(Fu);



