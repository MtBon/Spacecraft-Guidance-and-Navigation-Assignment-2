% Spacecraft Guidance and Navigation (2023/2024)
% Assigment #2
% Exercise #3
% Author: Matteo Bono
%
%% Point 1

clc; clearvars; close all

addpath('sgp4');
addpath('functions');
addpath('kernels');
addpath('..\mice\');
% Load spice kernels
cspice_furnsh('assignment02.tm');

matlab_graphics

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Satellite 1
Sat1.name = 'Mango';
Sat1.ID = 36599;

% Mean Initial State for Satellite Mango
Sat1.r0_mean = [4622.232026629; 5399.3369588058; -0.0212138165769957];
Sat1.v0_mean = [0.812221125483763; -0.721512914578826; 7.42665302729053];
Sat1.x0_mean = [Sat1.r0_mean; Sat1.v0_mean];



% Covariance Matrix
 P0_ref = [   [5.6e-7  3.5e-7  -7.1e-8 
           3.5e-7  9.7e-7   7.6e-8      
          -7.1e-8  7.6e-8   8.1e-8]            zeros(3)

           zeros(3)        diag([2.8e-11 , 2.7e-11 , 9.6e-12])];



et_ref = cspice_str2et('2010-08-12T05:27:29.114');   %ET

 % Open window
 t0 = '2010-08-12T05:30:00.000';
 et0 = cspice_str2et(t0);

 % Close window
 tf = '2010-08-12T6:30:00.000';
 etf = cspice_str2et(tf);

 %Time window with 5s step
 et = et0 : 5 : etf;

Station2.name = 'SVALBARD';
Station2.lat = 78.229772;     %[Deg]
Station2.lon = 15.407786;   %[Deg]
Station2.alt = 458;         %[m]
Station2.sigma_az = 125 * 1e-3;    %[Deg]
Station2.sigma_el = 125 * 1e-3;    %[Deg]
Station2.sigma_rho = 0.01;  %[Km]
Station2.min_el = 5;       %[Deg]

% Parameters of FFRF
FFRF.sigma_az = 1;      % [Deg]
FFRF.sigma_el = 1;      % [Deg]
FFRF.sigma_rho = 1 * 1e-5; %[Km]


%% Point 1a

% First state at t0
[Sat1.xx(1,:), ~, ~ ] = keplerian_propagator( et_ref,Sat1.x0_mean, et0 , 'Earth');


for j = 2:length(et)

%State 
[ Sat1.xx(j,:), ~,~] = keplerian_propagator(et(j-1), Sat1.xx(j-1,:), et(j) , 'Earth');

end


% Compute antenna angles, satellite range wrt Station 2
[Station2.sat1.rho, Station2.sat1.azimuth, Station2.sat1.elevation] = pointing(Station2.name,Sat1.xx(:,1:3)',Sat1.xx(:,4:6)',et);

%% Point 1b

%TLE for Mango
Sat1.longstr1 = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
Sat1.longstr2 = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';

% TLE for Tango
Sat2.longstr1 = '1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996';
Sat2.longstr2 = '2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55';

sat1rec = twoline2rv(Sat1.longstr1, Sat1.longstr2, typerun,'e', opsmode, whichconst);
sat2rec = twoline2rv(Sat2.longstr1, Sat2.longstr2, typerun,'e', opsmode, whichconst);

% Mango
[year,mon,day,hr,min,sec] = invjday(sat1rec.jdsatepoch, sat1rec.jdsatepochf);
sat1_tle_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat1_epoch_et = cspice_str2et(sat1_tle_epoch_str);

% Tango
[year,mon,day,hr,min,sec] = invjday(sat2rec.jdsatepoch, sat2rec.jdsatepochf);
sat2_tle_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat2_epoch_et = cspice_str2et(sat2_tle_epoch_str);


% Set nutation corrections parameters 
ddpsi = -0.073296 * arcsec2rad; %  [rad]
ddeps = -0.009373 * arcsec2rad; %  [rad]

% Initialize
Sat1.reci = zeros(3,length(et));
Sat1.veci = zeros(3,length(et));
Sat2.reci = zeros(3,length(et));
Sat2.veci = zeros(3,length(et));

% Loop over epochs
for i = 1:length(et)

    % SGP4 propagation for both satellites
    tsince1 = (et(i) - sat1_epoch_et)/60.0;   % time from TLE epoch of Mango in minutes
    tsince2 = (et(i) - sat2_epoch_et)/60.0;   % time from TLE epoch of Tango in minutes
    [~,Sat1.rteme,Sat1.vteme] = sgp4(sat1rec,  tsince1);
    [~,Sat2.rteme,Sat2.vteme] = sgp4(sat2rec,  tsince2);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(et(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [Sat1.reci(:,i), Sat1.veci(:,i)] = teme2eci(Sat1.rteme ,Sat1.vteme, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
    [Sat2.reci(:,i), Sat2.veci(:,i)] = teme2eci(Sat2.rteme ,Sat2.vteme, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);


end

GM = cspice_bodvrd('Earth','GM',1);



% Compute antenna angles, satellite range and range-rate wrt Station 2
[Station2.sat1.rho, Station2.sat1.azimuth, Station2.sat1.elevation] = pointing(Station2.name,Sat1.reci,Sat1.veci,et);

% Station 2 with nosie
i_visibility_noise_station2 = Station2.sat1.elevation*cspice_dpr()  > Station2.min_el ;
Station2.sat1.visibility_noise = et(i_visibility_noise_station2);
t = visibility(Station2.sat1.visibility_noise,Station2.name,Sat1.name,5);

% Model for noise
% Station 2 [Km - Deg - Deg]
Station2.sat1.rho_noise = mvnrnd(Station2.sat1.rho(i_visibility_noise_station2),ones(1,length(Station2.sat1.visibility_noise))*Station2.sigma_rho^2);
Station2.sat1.azimuth_noise = mvnrnd(Station2.sat1.azimuth(i_visibility_noise_station2)*cspice_dpr(),ones(1,length(Station2.sat1.visibility_noise))*Station2.sigma_az^2);
Station2.sat1.elevation_noise = mvnrnd(Station2.sat1.elevation(i_visibility_noise_station2)*cspice_dpr(),ones(1,length(Station2.sat1.visibility_noise))*Station2.sigma_el^2);

y_meas = [Station2.sat1.rho_noise;Station2.sat1.azimuth_noise;Station2.sat1.elevation_noise]; %[Km - Deg - Deg]

% Mean motion 
Sat1.n = sqrt(GM/norm(Sat1.reci(:,1))^3);


%% Point 1c

% Unscented Transform
alpha = 0.1;
n = 6;
lambda = alpha^2 * n;
beta = 2;


% Initialize Sigma points vector
 x_sigma_sat1 = zeros(6,13);
 
% Initialize vector containing last propagation of sigma points 
xx_sigma_sat1 = zeros(6,13);

sigma_mean_sat1 = zeros(6,2);
sigma_cov_sat1 = zeros(6,6,2);

% Weights for Mean
W_m = zeros(1,2*n + 1);
W_m(1) =  1 - n/(alpha^2*n);
W_m(2:end) = 1/(2 * alpha^2 * n);

% Weights for Covariance
W_c = zeros(1,2*n + 1);
W_c(1) = (2 - alpha^2 + beta) - n/(alpha^2 * n);
W_c(2:end) = 1/(2 * alpha^2 * n);


% Initialize first iteration values
 sigma_cov_sat1(:,:,1) = P0_ref*1e4;
 sigma_mean_sat1(:,1) = Sat1.x0_mean;
 
% First iteration to reach t0 from t_ref
P_sq1 = sqrtm(lambda * sigma_cov_sat1(:,:,1)) ;


    for i = 1 : 2 * n + 1

        % Compute Sigma Points
         if(i==1)

            x_sigma_sat1(:,i) = sigma_mean_sat1(:,1) ;
 

         elseif(i>1 && i<8)

            x_sigma_sat1(:,i) = sigma_mean_sat1(:,1) + P_sq1(:,i-1);
    
         else
    
            x_sigma_sat1(:,i) = sigma_mean_sat1(:,1) - P_sq1(:,i-n-1);
         end

        % Propagate Sigma points
        [xx_sigma_sat1(:,i), ~, ~ ] = keplerian_propagator(et_ref,x_sigma_sat1(:,i), et0, 'Earth');

        % Sample Mean
        sigma_mean_sat1(:,2) = sigma_mean_sat1(:,2) + W_m(i) * xx_sigma_sat1(:,i); 
        

    end

   
    for i = 1: 2*n+1

    sigma_cov_sat1(:,:,2) = sigma_cov_sat1(:,:,2) + W_c(i) * ((xx_sigma_sat1(:,i) - sigma_mean_sat1(:,2)) * (xx_sigma_sat1(:,i) - sigma_mean_sat1(:,2))');
    
    end

    % Update Variance matrices
    P_sq1 = sqrtm(lambda * sigma_cov_sat1(:,:,2)) ;


P0 = sigma_cov_sat1(:,:,end);

%Noise Matrix
R = diag([Station2.sigma_rho^2,Station2.sigma_az^2 , Station2.sigma_el^2]);

% Full vector of times starting from t0
time = [et0 ,Station2.sat1.visibility_noise];

% Unscented Kalman Filter
[Sat1.xk_up,Sat1.Pk_up,Sat1.sigma_pos,Sat1.sigma_vel] = UKF(time,P0,R,[Sat1.reci(:,1);Sat1.veci(:,1)],y_meas,Station2,Sat1.n,0);


State = [Sat1.reci(:,i_visibility_noise_station2); Sat1.veci(:,i_visibility_noise_station2)];

Sat1.err_pos = Sat1.xk_up(1:3,2:end) - State(1:3,:);
Sat1.norm_err_pos = vecnorm(Sat1.err_pos,2);

Sat1.err_vel = Sat1.xk_up(4:6,2:end) - State(4:6,:);
Sat1.norm_err_vel = vecnorm(Sat1.err_vel,2);

%% Point 2a



% Relative State at t0 in ECI frame
Dx0_eci = [Sat2.reci(:,1) ; Sat2.veci(:,1)] - [Sat1.reci(:,1) ; Sat1.veci(:,1)];

% Relative state in LVLH Frame
[Dr0, Dv0] = ECI_2_LVLH(Dx0_eci(1:3),Dx0_eci(4:6));


%% Point 2b
Dxf = zeros(6,length(et));
% First relative position
Dxf(:,1) = [Dr0 ; Dv0];

% Perform Propagation of CW equations
for j = 2 : length(et)

[Dxf(:,j) ,~, tt] = CW_propagation(et(j-1), Dxf(:,j-1) ,et(j),Sat1.n);

end

% Retrieve FFRF measurements from Geometry
[Rho,Azimuth,Elevation] = FFRF_measurements(Dxf);

% Noise 
Azimuth = mvnrnd(Azimuth*cspice_dpr,ones(1,length(et)) * FFRF.sigma_az^2);
Elevation = mvnrnd(Elevation*cspice_dpr, ones(1,length(et)) * FFRF.sigma_el^2);
Rho = mvnrnd(Rho,ones(1,length(et)) * FFRF.sigma_rho^2);

%% Point 2c

post_visibility = et > t;
% Time vector of 20 min after the time window
interval = et(t+1):5:et(t+1) + 20*60;

% Initial covariance
DP0 = diag([0.01, 0.01, 0.1, 0.0001, 0.0001, 0.001]) ;

% Noise Matrix for FFRF system
FFRF.R = diag([FFRF.sigma_rho^2,FFRF.sigma_az^2 , FFRF.sigma_el^2]);

% retrieve measurements during the time interval
rho_meas = Rho(t+1: t + length(interval));
az_meas = Azimuth(t+1: t + length(interval));
el_meas = Elevation(t+1: t + length(interval));
FFRF.meas = [rho_meas;az_meas;el_meas];

% Retrieve state in the initial time
DeltaX0 = Dxf(:,t+1);

% Unscented Kalman Filter
[Dxx,DP,Delta_sigma_pos,Delta_sigma_vel] = UKF(interval,DP0,FFRF.R,DeltaX0,FFRF.meas,Station2,Sat1.n,1); %Flag to 1 for FFRF

% Computing errors
rel_states = Dxf(:,t+1:t+length(interval));

rel_err_pos = Dxx(1:3,:) - rel_states(1:3,:);
rel_norm_err_pos = vecnorm(rel_err_pos,2);

rel_err_vel = Dxx(4:6,:) - rel_states(4:6,:);
rel_norm_err_vel = vecnorm(rel_err_vel,2);

%% Plots

% Plot passes (azimuth and elevation) w noise
figure()
scatter(Station2.sat1.azimuth_noise, Station2.sat1.elevation_noise,'b','filled','DisplayName', Sat1.name)
hold on;
plot([-180 180],[Station2.min_el Station2.min_el],'--')
axis([-180,180,0, 50])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('','Minimum Elevation');

% Errors for UKF Mango
figure()
subplot(2,1,1)
plot(Station2.sat1.visibility_noise,Sat1.norm_err_pos,'b','LineWidth',2);
hold on;
plot(Station2.sat1.visibility_noise,Sat1.sigma_pos(2:end),'--r','LineWidth',2);
legend('Position error','3 $\sigma$','Interpreter','latex');

subplot(2,1,2)
plot(Station2.sat1.visibility_noise,Sat1.norm_err_vel,'b','LineWidth',2);
hold on;
plot(Station2.sat1.visibility_noise,Sat1.sigma_vel(2:end),'--r','LineWidth',2);
legend('Velocity error','3 $\sigma$','Interpreter','latex');

% Errors for UKF FFRF system
figure()
subplot(2,1,1)
plot(interval,rel_norm_err_pos,'b','LineWidth',2);
hold on;
plot(interval,Delta_sigma_pos,'--r','LineWidth',2);
legend('Position error','3 $\sigma$','Interpreter','latex');
title('FFRF Position');
subplot(2,1,2)
plot(interval,rel_norm_err_vel,'b','LineWidth',2);
hold on;
plot(interval,Delta_sigma_vel,'--r','LineWidth',2);
legend('Velocity error','3 $\sigma$','Interpreter','latex');
title('FFRF Velocity');



%% Functions
function [r_LVLH, v_LVLH] = ECI_2_LVLH(r,v)


% Compute Row vectors
h = cross(r,v);
i_ax = r/norm(r);
k_ax = h/norm(h);
j_ax = cross(k_ax,i_ax);

% Rotation Matrix
rotm_LVLH = [i_ax, j_ax, k_ax];

% Position in LVLH Frame
r_LVLH = rotm_LVLH * r;

% Derivative of Matrix

d_ez = 1/norm(r) * (v - dot(r,v) * r);
d_ey = zeros(3,1);
d_ex = cross(h,d_ez);

dot_rotm_LVLH = [d_ex, d_ey , d_ez]';

v_LVLH = dot_rotm_LVLH * r + rotm_LVLH * v;


end

function t = visibility(et,stationname,satname,tstep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluate the boundaries of the visibility windows given as
% input the visibility times with a fixed step
%
%   INPUT:      -et: Time of visibility in ET   [s], [1,n]
%               -tstep: Given step of et, (1 minute if not specified)
%
%   OUTPUT:     -t: Index where visibility window finish
%
%
%   AUTHOT: Matteo Bono
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if(nargin<4)
    tstep = 60;
end

fprintf('@%s for #%s: Start Visibility: %s\n',stationname,satname,cspice_et2utc(et(1),'C',3));

flag=true;
k = 2;
while(flag==true && k<length(et)+1)

    if((et(k) - et(k-1)) == tstep)
        
    else
        flag = false;
        fprintf('@%s for #%s: End of Visibility: %s\n',stationname,satname,cspice_et2utc(et(k-1),'C',3));
        fprintf('--------------------------------------------------------------------------\n');
       
       

    end
    

    if(flag == false && k<length(et))
        flag = true;
        fprintf('@%s for #%s: Start Visibility: %s\n',stationname,satname,cspice_et2utc(et(k),'C',3));
    end

    if(flag==true && k==length(et))
        %flag = false;
        fprintf('@%s for #%s: End of Visibility: %s\n',stationname,satname,cspice_et2utc(et(k),'C',3));
        fprintf('--------------------------------------------------------------------------\n');
        t=k;
    end
k = k+1;
end
end

function [xf, tt, xx] = keplerian_propagator(t0, x0, t1 , attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 


    GM = cspice_bodvrd(attractor, 'GM', 1);


options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-9; ones(3,1)*1e-12]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 t1-t0], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end

function [xf, tt, xx] = keplerian_propagator_J2(t0, x0, t1 , attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 


    GM = cspice_bodvrd(attractor, 'GM', 1);


options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-9; ones(3,1)*1e-12]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs_J2(t,x,GM), [0 t1-t0], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end

function [dxdt] = keplerian_rhs_J2(t, x, GM)
%KEPLERIAN_RHS_J2  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator with J2 perturbation.
%
%
% Author
%   Name: Matteo
%   Surname: Bono
%   University: Politecnico di Milano 
%   
%   % 
% Inputs:
%   t   : [ 1, 1] epoch 
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

rotm = cspice_pxform('J2000','ITRF93',t);
radii = cspice_bodvrd('Earth','RADII',3);
pos_ECEF = rotm * rr;
J2 = 0.0010826269;
aj2 = 3/2 * GM * J2 * pos_ECEF/norm(pos_ECEF)^3 * (radii(1)/norm(pos_ECEF))^2 .* ( 5 * (pos_ECEF(3)/norm(pos_ECEF))^2 - [1;1;3]);
rotm_ECI = cspice_pxform('ITRF93','J2000',t);
aj2 = rotm_ECI * aj2;
% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2) + aj2;

end

function [dxdt] = keplerian_rhs(~, x, GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 24/10/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2);

end

function matlab_graphics()
%{
Set graphical values for better looking plots
%}

% interpreter:

set(0, 'defaultTextInterpreter', 'tex')
set(0, 'defaultAxesTickLabelInterpreter', 'tex')

% figure properties:

%colors
set(0, 'defaultFigureColormap',turbo(256));
set(0, 'defaultFigureColor', [1; 1; 1]);
% grid
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')

% surfaces:
% transparency
set(0, 'defaultSurfaceEdgeAlpha', 0.3);

% lines:

defaultLineWidth = 1.5;

% plots
set(0,'defaultLineLineWidth', defaultLineWidth);
% stairs
set(0,'defaultStairLineWidth', defaultLineWidth); % needs a different command for no reason apparently
% ylines
% set(0,'defaultYLineLineWidth', defaultLineWidth)
% legend:
set(0, 'defaultLegendLocation','best');
set(0, 'defaultLegendFontSize',7);

% axes:
% grid 
% set(0, 'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
% set(0, 'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
% font
set(0, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');

% color
set(0, 'defaultAxesColor', 'none');
% fontSize
set(0,'defaultAxesFontSize',16);
end

function [rho, azimuth, elevation] = pointing(stationName,rr_ECI,vv_ECI,et)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to compute the measurement of Range, Azimuth and
% Elevation of the propagated state of a Satellite with respect to a given
% Ground Station
%
%
%   INPUT:  -Station Name
%           -rr_ECI : propagated position state of the Spacecraft   [3,n]
%           -vv_ECI : propagated velocity state of the Spacecraft   [3,n]
%           -et     : vector of times in seconds                    [n]
%
%
%   OUTPUT: -Rho:       range measured                              [1,n]
%           -Azimuth:   Azimuth angle in Degree                     [1,n]
%           -Elevation: Elevation angle in Degree                   [1,n]
%
% Author:   Matteo Bono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
n = length(et);
ROT_ECI2TOPO = zeros(6,6,n);
rv_station_eci = zeros(6,n);
rv_station_sat_eci = zeros(6,n);
rv_station_sat_topo = zeros(6,n);
rho = zeros(1,n);
azimuth = zeros(1,n);
elevation = zeros(1,n);

topoFrame = [stationName, '_TOPO'];

for i = 1: length(et)

% Transformation from ECI to topocentric frame
ROT_ECI2TOPO(:,:,i) = cspice_sxform('J2000', topoFrame, et(i));

% Compute spacecraft position in ECI
rv_station_eci(:,i) = cspice_spkezr(stationName, et(i), 'J2000', 'NONE', 'EARTH');

% Compute station-satellite vector in ECI
rv_station_sat_eci(:,i) = [rr_ECI(:,i); vv_ECI(:,i)] - rv_station_eci(:,i);

% Convert state into topocentric frame
rv_station_sat_topo(:,i) = ROT_ECI2TOPO(:,:,i) * rv_station_sat_eci(:,i);

% Compute range, azimuth and elevation using cspice_reclat
[rho(i), azimuth(i),  elevation(i)] = cspice_reclat(rv_station_sat_topo(1:3,i));

end
end

function [xf , xx, tt] = CW_propagation(et0, x0,etf,n)
%
% This function propagates the relative motion of two satellites with CW
% equation
options = odeset('reltol', 1e-12, 'abstol', 1e-12*ones(6,1));

% Perform integration
[tt, xx] = ode78(@(t,x) CW_eq(t,x,n),[0 etf-et0], x0, options);

% Extract state vector 
xf = xx(end,1:6);

end

function [dxdt] = CW_eq(~,x,n)

% Initialize right-hand-side
dxdt = zeros(6,1);

dxdt(1:3) = x(4:6);
dxdt(4) = 3 * n^2 * x(1) + 2 * n * x(5)^2;
dxdt(5) = -2 * n * x(4);
dxdt(6) = -n^2 * x(3);

end

function [Rho,Azimuth,Elevation] = FFRF_measurements(Dx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function retrieve the measurements of Rho , Azimuth and Elevation of
% the relative state vector between a first Spacecraft with FFRF on board
% and a second one. The measurements are acquire in the LVLH Reference
% frame solidal to the S/C N.1
%
%   INPUT: -Dx : Propagated Relative state                      [6,n]
%          -FFRF : Struct with measurements noise parameters
%
%   OUTPUT: -Rho        [1,n]
%           -Azimuth    [1,n]
%           -Elevation  [1,n]
%
% Author : Matteo Bono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length of Time span
n = size(Dx,2);

% Initializations
Azimuth = zeros(1,n);
Elevation = zeros(1,n);
Rho = zeros(1,n);

for i = 1 : n

    Azimuth(i) = atan2(Dx(2,i),Dx(1,i));
    Elevation(i) = atan2(Dx(3,i),sqrt(Dx(1,i)^2 + Dx(2,i)^2));
    Rho(i) = norm(Dx(1:3,i));

end




end

function [xk_up,Pk_up,sigma_pos,sigma_vel] = UKF(time,P0,R,x0,y_meas,station,N,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function represent the Unscented Kalman Filter.
%
%       INPUT: -time: vector of time
%              -P0: initial covariance matrix
%              -R: Diagonal Measurements Noise matrix 
%              -x0: Initial State [6,1]
%              -y_meas: Measurements acquired [3,length(time)]
%              -station: Struct of Station
%              -flag: 1 for LVLH relative propagation, 0 for ECI
%              -N: Mean motion of Satellite 1
%
%      OUTPUT: -xk_up: estimated states [6,length(time)]
%              -Pk_up: estimated covariance [6,6,length(time)]
%              -sigma_pos: 3sigma for Pos
%              -sigma_vel: 3sigma for Vel
%
% Author: Matteo Bono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametersd of UT
alpha = 0.1;
n = 6;
lambda = alpha^2 * n;
beta = 2;

% Weights for Mean
W_m = zeros(1,2*n + 1);
W_m(1) =  1 - n/(alpha^2*n);
W_m(2:end) = 1/(2 * alpha^2 * n);

% Weights for Covariance
W_c = zeros(1,2*n + 1);
W_c(1) = (2 - alpha^2 + beta) - n/(alpha^2 * n);
W_c(2:end) = 1/(2 * alpha^2 * n);

% Initialize Sigma points vector
x_sigma = zeros(n,13);
  
% Initialize vector containing last propagation of sigma points 
xx_sigma = zeros(n,13);

% Initialization vector of mean values 
xk_m = zeros(n,length(time)-1); 
yk_m = zeros(3,length(time)-1);
yyk = zeros(3,2*n+1);
xk_up = zeros(n,length(time));


% Initialization vector of covariance matrices for every step
Pk_m = zeros(n,n,length(time)-1);
Pyy = zeros(3,3,length(time)-1);
Pxy = zeros(n,3,length(time)-1);

% Initialization of measurements
rho_meas = zeros(1,2*n+1);
azimuth_meas = zeros(1,2*n+1);
elevation_meas = zeros(1,2*n+1);

sigma_pos = zeros(1,length(time));
sigma_vel = zeros(1,length(time));


% Initialize first iteration values
 Pk_up(:,:,1) = P0;
 xk_up(:,1) = x0;
 sigma_pos(:,:,1) =  3 * sqrt(max(eig(P0(1:3,1:3))));
 sigma_vel(:,:,1) =  3 * sqrt(max(eig(P0(4:6,4:6))));

% First iteration
P_sq1 = sqrtm(lambda * Pk_up(:,:,1));


for j = 2 : length(time)

    for i = 1 : 2 * 6 + 1

        % Compute Sigma Points
         if(i==1)

            x_sigma(:,i) = xk_up(:,j-1) ;

         elseif(i>1 && i<8)

            x_sigma(:,i) = xk_up(:,j-1) + P_sq1(:,i-1);
            
         else
    
            x_sigma(:,i) = xk_up(:,j-1) - P_sq1(:,i-n-1);

         end

     if(flag==0)
        % Propagate Sigma points
        [xx_sigma(:,i), ~, ~ ] = keplerian_propagator_J2(time(j-1),x_sigma(:,i), time(j), 'Earth');
        
        % Sigma points in Measurement space %[Km - Rad - Rad]
       
        [rho_meas(i), azimuth_meas(i), elevation_meas(i)] = pointing(station.name,xx_sigma(1:3,i),xx_sigma(4:6,i),time(j));
     else

         [xx_sigma(:,i) ,~, ~] = CW_propagation(time(j-1), x_sigma(:,i) ,time(j),N);

         % Retrieve FFRF measurements from Geometry
         [rho_meas(i),azimuth_meas(i),elevation_meas(i)] = FFRF_measurements(xx_sigma(:,i));

     end
        % Sample Mean
        xk_m(:,j-1) = xk_m(:,j-1) + W_m(i) * xx_sigma(:,i); 

        % [Km - Deg - Deg]
        yyk(:,i) = [rho_meas(i);  azimuth_meas(i)*cspice_dpr; elevation_meas(i)*cspice_dpr];   

       % Sample Mean of measurements
        yk_m(:,j-1) = yk_m(:,j-1) + W_m(i) * yyk(:,i); 
        
    end

   
    for i = 1: 2*6+1

    Pk_m(:,:,j-1) = Pk_m(:,:,j-1) + W_c(i) * ((xx_sigma(:,i) - xk_m(:,j-1)) * (xx_sigma(:,i) - xk_m(:,j-1))');
    
    % Measurement differences
    diff_rho = yyk(1,i) - yk_m(1,j-1);
    diff_az = angdiff(yyk(2,i)*cspice_rpd,yk_m(2,j-1)*cspice_rpd);
    diff_el = angdiff(yyk(3,i)*cspice_rpd,yk_m(3,j-1)*cspice_rpd);
    diff_meas = [diff_rho ; diff_az*cspice_dpr ; diff_el*cspice_dpr];
    
    Pyy(:,:,j-1) = Pyy(:,:,j-1) + W_c(i) * ((diff_meas) * (diff_meas'));
    Pxy(:,:,j-1) = Pxy(:,:,j-1) + W_c(i) * ((xx_sigma(:,i) - xk_m(:,j-1)) * (diff_meas'));

    end
    Pyy(:,:,j-1) = Pyy(:,:,j-1) + R;

    % Kalman Gain
    Kk = Pxy(:,:,j-1) / (Pyy(:,:,j-1));

    % Update State
    xk_up(:,j) = xk_m(:,j-1) + Kk * (y_meas(:,j-1) - yk_m(:,j-1));
    Pk_up(:,:,j) = Pk_m(:,:,j-1) - Kk * Pyy(:,:,j-1) * Kk';

    % Force simmetry
    Pk_up(:,:,j) = (Pk_up(:,:,j) + Pk_up(:,:,j)')/2;

    % Update Matrix for sigma points
    P_sq1 = sqrtm(lambda * Pk_up(:,:,j));

    % 3 Sigma for Pos and Vel
    sigma_pos(j) = 3 * sqrt(max(eig( Pk_up(1:3,1:3,j))));
    sigma_vel(j) = 3 * sqrt(max(eig( Pk_up(4:6,4:6,j)))); 

    
  
 
end

end

