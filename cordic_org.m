% This program computes cos(theta) & sin(theta) using CORDIC algorithm, where 'theta'
% is in degrees & moreover here, 'theta' should lie between 0 & 90 degrees
% for this program to work!!!
% References: Volder, J.E., 1959. The CORDIC trigonometric computing technique. IRE Transactions on electronic computers, (3), pp.330-334.
% Written by Vineel Kumar Veludandi

clear all
clc
clear all
max_iter = 20; % maximum number of iterations allowed in CORDIC algorithm/PRECISION
theta = 40; % in degrees

%--------------------------------------------------------------------------
%                      LOOKUP TABLE / PRE-COMPUTATIONS
angle_deg=zeros(1,max_iter);
itr_num = 0:max_iter-1;
for i1=1:max_iter
    angle_deg(i1)=(180/pi)*atan(1/2^(i1-1));
end
multipliers = ones(1,max_iter+1);
for i1=1:max_iter
    multipliers(i1+1)=multipliers(i1)*1/sqrt(1+2^(-2*(i1-1)));
end
multipliers(1)=[];
%--------------------------------------------------------------------------
% initial coordinates of X & Y
x0=1; y0=0; 
for i2=1:max_iter
if theta>=0
    d=1;
    temp_var1 = x0 - y0*d*2^(-(i2-1));
    temp_var2 = y0 + x0*d*2^(-(i2-1));
    x0=temp_var1;
    y0=temp_var2;
    theta = theta - d*angle_deg(i2);
elseif theta<0
    d=-1;
    temp_var1 = x0 - y0*d*2^(-(i2-1));
    temp_var2 = y0 + x0*d*2^(-(i2-1)); 
    x0=temp_var1;
    y0=temp_var2;
    theta = theta - d*angle_deg(i2);
end
end
cosine_theta = x0*multipliers(max_iter);
sine_theta = y0*multipliers(max_iter);
fprintf('cosine(theta) is %f and sin(theta) is %f',cosine_theta,sine_theta)
