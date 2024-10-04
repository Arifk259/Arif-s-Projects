%Simulating Orbital Trajectory of ISS
%ENGPHYS 3NM4 - Final Project

%% Simple Two-Body Orbit - Euler's Method
close;
clear;
clc;
format long

G = 6.67408e-11;
me = 5.9722e24;

mu = G*me/1e9; %km^3/s^2

%Initial Values
x0 = 287.6240;
y0 = -6718.3;
z0 = 14.5677;

vx0 = 4.7718;
vy0 = 0.2202;
vz0 = 6.0350;

r = @(x,y,z) sqrt(x^2+y^2+z^2);
r0 = r(x0,y0,z0);

%Orbital Paramters
a = 6720.32608;
T = (2*pi*a^(3/2))/sqrt(mu); %period
h = 1; %time step
t = [0:h:5*T];
n = length(t);



%Differential Equation 
fx = @(vx) vx;
gx = @(x,y,z) (-mu*x)/r(x,y,z)^3;

fy = @(vy) vy;
gy = @(x,y,z) (-mu*y)/r(x,y,z)^3;

fz = @(vz) vz;
gz = @(x,y,z) (-mu*z)/r(x,y,z)^3;


%Euler's Method
x_eu = zeros(1,n);
vx_eu = zeros(1,n);
y_eu = zeros(1,n);
vy_eu = zeros(1,n);
z_eu = zeros(1,n);
vz_eu = zeros(1,n);

x_eu(1) = x0;
vx_eu(1) = vx0;
y_eu(1) = y0;
vy_eu(1) = vy0;
z_eu(1) = z0;
vz_eu(1) = vz0;

for i = 1:n-1
    vx_eu(i+1) = vx_eu(i) + gx(x_eu(i),y_eu(i),z_eu(i))*h;
    x_eu(i+1) = x_eu(i) + fx(vx_eu(i))*h;

    vy_eu(i+1) = vy_eu(i) + gy(x_eu(i),y_eu(i),z_eu(i))*h;
    y_eu(i+1) = y_eu(i) + fy(vy_eu(i))*h;

    vz_eu(i+1) = vz_eu(i) + gz(x_eu(i),y_eu(i),z_eu(i))*h;
    z_eu(i+1) = z_eu(i) + fz(vz_eu(i))*h;
end


figure(1)
eutwo = plot(x_eu,y_eu);
xlabel('x (km)');
ylabel('y (km)');
title('2D Orbital Trajectory of ISS Using Euler Method');

figure(2)
euthree = plot3(x_eu,y_eu,z_eu);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('3D Orbital Trajectory of ISS Using Euler Method');

%Euler-Cromer Method
x_ec = zeros(1,n);
vx_ec = zeros(1,n);
y_ec = zeros(1,n);
vy_ec = zeros(1,n);
z_ec = zeros(1,n);
vz_ec = zeros(1,n);

x_ec(1) = x0;
vx_ec(1) = vx0;
y_ec(1) = y0;
vy_ec(1) = vy0;
z_ec(1) = z0;
vz_ec(1) = vz0;

for i = 1:n-1
    vx_ec(i+1) = vx_ec(i) + gx(x_ec(i),y_ec(i),z_ec(i))*h;
    x_ec(i+1) = x_ec(i) + fx(vx_ec(i+1))*h;

    vy_ec(i+1) = vy_ec(i) + gy(x_ec(i),y_ec(i),z_ec(i))*h;
    y_ec(i+1) = y_ec(i) + fy(vy_ec(i+1))*h;

    vz_ec(i+1) = vz_ec(i) + gz(x_ec(i),y_ec(i),z_ec(i))*h;
    z_ec(i+1) = z_ec(i) + fz(vz_ec(i+1))*h;
end

figure(3)
ectwo = plot(x_ec,y_ec);
xlabel('x (km)');
ylabel('y (km)');
title('2D Orbital Trajectory of ISS Using Euler-Cromer Method');

figure(4)
ecthree = plot3(x_ec,y_ec,z_ec);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('3D Orbital Trajectory of ISS Using Euler-Cromer Method');


%Runge Kutta Method
x_rk = zeros(1,n);
vx_rk = zeros(1,n);
y_rk = zeros(1,n);
vy_rk = zeros(1,n);
z_rk = zeros(1,n);
vz_rk = zeros(1,n);

x_rk(1) = x0;
vx_rk(1) = vx0;
y_rk(1) = y0;
vy_rk(1) = vy0;
z_rk(1) = z0;
vz_rk(1) = vz0;

for i = 1:n-1
    k1x = fx(vx_rk(i));
    k1vx = gx(x_rk(i), y_rk(i), z_rk(i));
    k1y = fy(vy_rk(i));
    k1vy = gy(x_rk(i), y_rk(i), z_rk(i));
    k1z = fz(vz_rk(i));
    k1vz = gz(x_rk(i), y_rk(i), z_rk(i));

    k2x = fx(vx_rk(i) + 1/2*k1vx*h);
    k2vx = gx(x_rk(i) + 1/2*k1x*h, y_rk(i) + 1/2*k1y*h, z_rk(i) + 1/2*k1z*h);
    k2y = fy(vy_rk(i) + 1/2*k1vy*h);
    k2vy = gy(x_rk(i) + 1/2*k1x*h, y_rk(i) + 1/2*k1y*h, z_rk(i) + 1/2*k1z*h);
    k2z = fz(vz_rk(i) + 1/2*k1vz*h);
    k2vz = gz(x_rk(i) + 1/2*k1x*h, y_rk(i) + 1/2*k1y*h, z_rk(i) + 1/2*k1z*h);

    k3x = fx(vx_rk(i) + 1/2*k2vx*h);
    k3vx = gx(x_rk(i) + 1/2*k2x*h, y_rk(i) + 1/2*k2y*h, z_rk(i) + 1/2*k2z*h);
    k3y = fy(vy_rk(i) + 1/2*k2vy*h);
    k3vy = gy(x_rk(i) + 1/2*k2x*h, y_rk(i) + 1/2*k2y*h, z_rk(i) + 1/2*k2z*h);
    k3z = fz(vz_rk(i) + 1/2*k2vz*h);
    k3vz = gz(x_rk(i) + 1/2*k2x*h, y_rk(i) + 1/2*k2y*h, z_rk(i) + 1/2*k2z*h);

    k4x = fx(vx_rk(i) + k3vx*h);
    k4vx = gx(x_rk(i) + k3x*h, y_rk(i) + k3y*h, z_rk(i) + k3z*h);
    k4y = fy(vy_rk(i) + k3vy*h);
    k4vy = gy(x_rk(i) + k3x*h, y_rk(i) + k3y*h, z_rk(i) + k3z*h);
    k4z = fz(vz_rk(i) + k3vz*h);
    k4vz = gz(x_rk(i) + k3x*h, y_rk(i) + k3y*h, z_rk(i) + k3z*h);

    x_rk(i+1) = x_rk(i) + 1/6*(k1x + 2*k2x + 2*k3x + k4x)*h;
    vx_rk(i+1) = vx_rk(i) + 1/6*(k1vx + 2*k2vx + 2*k3vx + k4vx)*h;
    y_rk(i+1) = y_rk(i) + 1/6*(k1y + 2*k2y + 2*k3y + k4y)*h;
    vy_rk(i+1) = vy_rk(i) + 1/6*(k1vy + 2*k2vy + 2*k3vy + k4vy)*h;
    z_rk(i+1) = z_rk(i) + 1/6*(k1z + 2*k2z + 2*k3z + k4z)*h;
    vz_rk(i+1) = vz_rk(i) + 1/6*(k1vz + 2*k2vz + 2*k3vz + k4vz)*h;
end


figure(5)
rktwo = plot(x_rk, y_rk);
xlabel('x (km)');
ylabel('y (km)');
title('2D Orbital Trajectory of ISS Using RK4 Method');

figure(6)
rkthree = plot3(x_rk, y_rk, z_rk);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('3D Orbital Trajectory of ISS Using RK4 Method');

%Comparison of Methods
figure(7)
eutwocomp = plot(x_eu,y_eu);
hold on
ectwocomp = plot(x_ec,y_ec);
hold on
rktwocomp = plot(x_rk,y_rk);
legend('Euler Method','Euler-Cromer Method','4th Order RK');
title('Comparison of 2D Trajectories of ISS');
xlabel('x (km)');
ylabel('y (km)');

figure(8)
euthreecomp = plot3(x_eu,y_eu,z_eu);
hold on
ecthreecomp = plot3(x_ec,y_ec,z_ec);
hold on
rkthreecomp = plot3(x_rk,y_rk,z_rk);
legend('Euler Method','Euler-Cromer Method','4th Order RK');
title('Comparison of 3D Trajectories of ISS');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');



%% Natural Cubic Spline
close;
clear;
clc;
format long
syms X

altitude = [120;125;130;135;140;145;150;155;160;165;170;175;180;185;190;195;200
     210;220;230;240;250;260;270;280;290;300;310;320;330;340;350;360;370
     380;390;400;420;440;460;480;500;520;540;560;580;600;620;640
     660;680;700;720;740;760;780;800];

th_den = 1e12*[1.72e-11; 1.04e-11; 6.79e-12; 4.66e-12; 3.32e-12; 2.44e-12; 1.83e-12; 1.41e-12
     1.10e-12; 8.65e-13; 6.91e-13; 5.58e-13; 4.54e-13; 3.72e-13; 3.07e-13; 2.55e-13
     2.13e-13; 1.51e-13; 1.09e-13; 7.94e-14; 5.87e-14; 4.39e-14; 3.32e-14; 2.53e-14
     1.95e-14; 1.51e-14; 1.18e-14; 9.22e-15; 7.28e-15; 5.77e-15; 4.60e-15; 3.68e-15
     2.95e-15; 2.38e-15; 1.92e-15; 1.56e-15; 1.26e-15; 8.40e-16; 5.63e-15; 3.80e-16
     2.59e-16; 1.78e-16; 1.23e-16; 8.60e-17; 6.07e-17; 4.34e-17; 3.15e-17; 2.32e-17
     1.74e-17; 1.33e-17; 1.04e-17; 8.29e-18; 6.74e-18; 5.60e-18; 4.73e-18; 4.06e-18
     3.54e-18];


n = length(altitude);

for i = 1:n-1
    h(i) = altitude(i+1)-altitude(i);
end
a = th_den;

for i = 2:n-1
    alpha(i) = (3/(h(i-1)*h(i)))*(a(i+1)*h(i-1) - a(i)*(altitude(i+1)-altitude(i-1)) + a(i-1)*h(i));
end
alpha;

l = zeros(n,1);
z = zeros(n,1);
mu = zeros(n,1);
l(1) = 1;

for i = 2:n-1
    l(i) = 2*(altitude(i+1)-altitude(i-1)) - h(i-1)*mu(i-1);
    mu(i) = h(i)/l(i);
    z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i);
end

l(n) = 1;
z(n) = 0;
c = zeros(n,1);

for i = n-1:-1:1
    c(i) = z(i) - mu(i)*c(i+1);
    b(i) = (a(i+1)-a(i))/h(i) - h(i)*(c(i+1) + 2*c(i))/3;
    d(i) = (c(i+1) - c(i))/(3*h(i));
end

%symbolic matrix definition
s = sym('s', [3,1]);

for i = 1:n-1
    s(i) = a(i) + b(i)*(X-altitude(i)) + c(i)*(X-altitude(i))^2 + d(i)*(X-altitude(i))^3;
end

%Plotting
for i = 1:n-1
    xlin = linspace(altitude(i), altitude(i+1));
    splot = eval(subs(s(i), xlin));
    plot(xlin, splot);
    hold on
end

plot(altitude,th_den,'o', xlin, splot)
xlabel('Altitude (km)');
ylabel('Atmospheric Density (kg/km^3)');
title('Atmospheric Density of the Thermosphere at Different Altitudes');


%% Adding Perturbations
clc;
G = 6.67408e-11;
me = 5.9722e24;
re = 6378;

mu = G*me/1e9; %km^3/s^2

%Initial Values
x0 = 287.6240;
y0 = -6718.3;
z0 = 14.5677;

vx0 = 4.7718;
vy0 = 0.2202;
vz0 = 6.0350;

r = @(x,y,z) sqrt(x^2+y^2+z^2);
r0 = r(x0,y0,z0);

%Orbital Paramters
a = 6720.32608;
T = (2*pi*a^(3/2))/sqrt(mu); %period
h = 1; %time step
t = [0:h:20*T]; %20 Periods
n = length(t);

%Drag Parameters
ms = 471727; %mass of satellite in kg
A = 2040.50/1e6; %area in km^2
Cd = 2;
rho = den(x0,y0,z0,altitude,s);


%Runge Kutta Vectors
x_rk = zeros(1,n);
vx_rk = zeros(1,n);
y_rk = zeros(1,n);
vy_rk = zeros(1,n);
z_rk = zeros(1,n);
vz_rk = zeros(1,n);

x_rk(1) = x0;
vx_rk(1) = vx0;
y_rk(1) = y0;
vy_rk(1) = vy0;
z_rk(1) = z0;
vz_rk(1) = vz0;

%Differential Equation 
fx = @(vx) vx;
gx = @(x,y,z,vx,vy,vz) (-mu*x)/r(x,y,z)^3 - (1/2*rho*Cd*A*sqrt(vx^2+vy^2+vz^2)*vx)/ms;

fy = @(vy) vy;
gy = @(x,y,z,vx,vy,vz) (-mu*y)/r(x,y,z)^3 - (1/2*rho*Cd*A*sqrt(vx^2+vy^2+vz^2)*vy)/ms;

fz = @(vz) vz;
gz = @(x,y,z,vx,vy,vz) (-mu*z)/r(x,y,z)^3 - (1/2*rho*Cd*A*sqrt(vx^2+vy^2+vz^2)*vz)/ms;


for i = 1:n-1
    rho = den(x_rk(i),y_rk(i),z_rk(i),altitude,s)
    k1x = fx(vx_rk(i));
    k1vx = gx(x_rk(i), y_rk(i), z_rk(i), vx_rk(i), vy_rk(i), vz_rk(i));
    k1y = fy(vy_rk(i));
    k1vy = gy(x_rk(i), y_rk(i), z_rk(i), vx_rk(i), vy_rk(i), vz_rk(i));
    k1z = fz(vz_rk(i));
    k1vz = gz(x_rk(i), y_rk(i), z_rk(i), vx_rk(i), vy_rk(i), vz_rk(i));

    k2x = fx(vx_rk(i) + 1/2*k1vx*h);
    k2vx = gx(x_rk(i) + 1/2*k1x*h, y_rk(i) + 1/2*k1y*h, z_rk(i) + 1/2*k1z*h, vx_rk(i) + 1/2*k1vx*h, vy_rk(i) + 1/2*k1vy*h, vz_rk(i) + 1/2*k1vz*h);
    k2y = fy(vy_rk(i) + 1/2*k1vy*h);
    k2vy = gy(x_rk(i) + 1/2*k1x*h, y_rk(i) + 1/2*k1y*h, z_rk(i) + 1/2*k1z*h, vx_rk(i) + 1/2*k1vx*h, vy_rk(i) + 1/2*k1vy*h, vz_rk(i) + 1/2*k1vz*h);
    k2z = fz(vz_rk(i) + 1/2*k1vz*h);
    k2vz = gz(x_rk(i) + 1/2*k1x*h, y_rk(i) + 1/2*k1y*h, z_rk(i) + 1/2*k1z*h, vx_rk(i) + 1/2*k1vx*h, vy_rk(i) + 1/2*k1vy*h, vz_rk(i) + 1/2*k1vz*h);

    k3x = fx(vx_rk(i) + 1/2*k2vx*h);
    k3vx = gx(x_rk(i) + 1/2*k2x*h, y_rk(i) + 1/2*k2y*h, z_rk(i) + 1/2*k2z*h, vx_rk(i) + 1/2*k2vx*h, vy_rk(i) + 1/2*k2vy*h, vz_rk(i) + 1/2*k2vz*h);
    k3y = fy(vy_rk(i) + 1/2*k2vy*h);
    k3vy = gy(x_rk(i) + 1/2*k2x*h, y_rk(i) + 1/2*k2y*h, z_rk(i) + 1/2*k2z*h, vx_rk(i) + 1/2*k2vx*h, vy_rk(i) + 1/2*k2vy*h, vz_rk(i) + 1/2*k2vz*h);
    k3z = fz(vz_rk(i) + 1/2*k2vz*h);
    k3vz = gz(x_rk(i) + 1/2*k2x*h, y_rk(i) + 1/2*k2y*h, z_rk(i) + 1/2*k2z*h, vx_rk(i) + 1/2*k2vx*h, vy_rk(i) + 1/2*k2vy*h, vz_rk(i) + 1/2*k2vz*h);

    k4x = fx(vx_rk(i) + k3vx*h);
    k4vx = gx(x_rk(i) + k3x*h, y_rk(i) + k3y*h, z_rk(i) + k3z*h, vx_rk(i) + k3vx*h, vy_rk(i) + k3vy*h, vz_rk(i) + k3vz*h);
    k4y = fy(vy_rk(i) + k3vy*h);
    k4vy = gy(x_rk(i) + k3x*h, y_rk(i) + k3y*h, z_rk(i) + k3z*h, vx_rk(i) + k3vx*h, vy_rk(i) + k3vy*h, vz_rk(i) + k3vz*h);
    k4z = fz(vz_rk(i) + k3vz*h);
    k4vz = gz(x_rk(i) + k3x*h, y_rk(i) + k3y*h, z_rk(i) + k3z*h, vx_rk(i) + k3vx*h, vy_rk(i) + k3vy*h, vz_rk(i) + k3vz*h);

    x_rk(i+1) = x_rk(i) + 1/6*(k1x + 2*k2x + 2*k3x + k4x)*h;
    vx_rk(i+1) = vx_rk(i) + 1/6*(k1vx + 2*k2vx + 2*k3vx + k4vx)*h;
    y_rk(i+1) = y_rk(i) + 1/6*(k1y + 2*k2y + 2*k3y + k4y)*h;
    vy_rk(i+1) = vy_rk(i) + 1/6*(k1vy + 2*k2vy + 2*k3vy + k4vy)*h;
    z_rk(i+1) = z_rk(i) + 1/6*(k1z + 2*k2z + 2*k3z + k4z)*h;
    vz_rk(i+1) = vz_rk(i) + 1/6*(k1vz + 2*k2vz + 2*k3vz + k4vz)*h;
end



%% Plotting
clc;

%3D and 2D radii of orbit
rvals = sqrt(x_rk.^2 + y_rk.^2 + z_rk.^2)-re;
rtwo = sqrt(x_rk.^2 + y_rk.^2)-re;

%Index of every 5 periods
tperiod = [0:h:5*T];
T1 = find(t==floor(5*T));
T2 = find(t==floor(10*T));
T3 = find(t==floor(15*T));
T4 = find(t==floor(20*T));

%Initial positions for every 5 periods
r1 = sqrt(x_rk(1)^2 + y_rk(1)^2 + z_rk(1)^2)*1000;
r2 = sqrt(x_rk(T1)^2 + y_rk(T1)^2 + z_rk(T1)^2)*1000;
r3 = sqrt(x_rk(T2)^2 + y_rk(T2)^2 + z_rk(T2)^2)*1000;
r4 = sqrt(x_rk(T3)^2 + y_rk(T3)^2 + z_rk(T3)^2)*1000;
r5 = sqrt(x_rk(T4)^2 + y_rk(T4)^2 + z_rk(T4)^2)*1000;

delr2 = r2-r1;
delr3 = r3-r2;
delr4 = r4-r3;
delr5 = r5-r4;
Tvals = [5,10,15,20];
delrvals = [delr2, delr3, delr4, delr5];


%2D Trajectory
figure(1)
traj = plot(x_rk, y_rk);
xlabel('x (km)');
ylabel('y (km)');
title('2D Trajectory of ISS Satellite');

%3D Trajectory
figure(2)
threeD_traj = plot3(x_rk,y_rk,z_rk);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('3D Trajectory of ISS Satellite')

%Period Plots
figure(3)
p1 = plot(tperiod,rvals(1:T1));
hold on
p2 = plot(tperiod, rvals(T1:T2));
hold on
p3 = plot(tperiod, rvals(T2:T3));
hold on
p4 = plot(tperiod, rvals(T3:T4));
legend('Period 1-5', 'Period 6-10', 'Period 11-15', 'Period 16-20');
xlabel('t (s)');
ylabel('altitude (km)');
title("Decay of Satellite's Altitude from Atmospheric Drag");

figure(4)
rplot = plot(tperiod, rtwo(1:T1));
hold on
rplot2 = plot(tperiod,rtwo(T1:T2));
hold on
rplot3 = plot(tperiod,rtwo(T2:T3));
hold on
rplot4 = plot(tperiod,rtwo(T3:T4));
legend('Period 1-5', 'Period 6-10', 'Period 11-15', 'Period 16-20');
xlabel('t (s)');
ylabel('altitude (km)');
title("Decay of Satellite's XY Position");

figure(5)
relplot = plot(Tvals,delrvals);
xlim([0,20]);
title('Change in Radial Distance of ISS');
ylabel('change in radial distance (m)');
xlabel('Periods');



%% Atmoshpheric Drag Function (from Thermosphere) for Interpolation

function density = den(x0,y0,z0,alt,sp)
    syms X
    r = sqrt(x0^2+y0^2+z0^2);
    re = 6378;
    h = r-re;
    j = 0;
    n = length(alt);
    for i = 1:n
        if alt(i) > h
            j = i-1;
            break
        else
            j = n-1;
        end
    end
    density = eval(subs(sp(j), h));
end




