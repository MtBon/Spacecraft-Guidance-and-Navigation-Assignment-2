% Spacecraft Guidance and Navigation (2023/2024)
% Assigment #2
% Exercise #2
% Author: Matteo Bono
%
%% Point 1

clc; clearvars; cspice_kclear; close all

% Replace the path below with your own installation path if not present as
% a subfolder
addpath('sgp4');
addpath('functions');
addpath('kernels');
% Load spice kernels
cspice_furnsh('assignment02.tm');

matlab_graphics

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Satellite
Sat1.name = 'Mango';
Sat1.ID = 36599;

% Mean Initial State for Satellite Mango
r0_mean_sat1 = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v0_mean_sat1 = [0.812221125483763; -0.721512914578826; 7.42665302729053];
x0_mean_sat1 = [r0_mean_sat1; v0_mean_sat1];

% Covariance Matrix
 P0 = [   [5.6e-7  3.5e-7  -7.1e-8 
           3.5e-7  9.7e-7   7.6e-8      
          -7.1e-8  7.6e-8   8.1e-8]            zeros(3)

           zeros(3)        diag([2.8e-11 , 2.7e-11 , 9.6e-12])];

 et_ref = cspice_str2et('2010-08-12T05:27:29.114');   %ET

 % Open visibility window
 t0 = '2010-08-12T05:30:00.000';
 et0 = cspice_str2et(t0);

 % Close visibility window
 tf = '2010-08-12T11:00:00.000';
 etf = cspice_str2et(tf);

 %Time window with 60s step
 et = et0 : 60 : etf;

% Stations
Station1.name = 'KOUROU';
Station1.lat = 5.25144;     %[Deg]
Station1.lon = -52.80466;   %[Deg]
Station1.alt = -14.67;      %[m]
Station1.sigma_az = 100 * 1e-3;    %[Deg]
Station1.sigma_el = 100 * 1e-3;    %[Deg]
Station1.sigma_rho = 0.01;  %[Km]
Station1.min_el = 10;       %[Deg]

Station2.name = 'SVALBARD';
Station2.lat = 78.229772;     %[Deg]
Station2.lon = 15.407786;   %[Deg]
Station2.alt = 458;         %[m]
Station2.sigma_az = 125 * 1e-3;    %[Deg]
Station2.sigma_el = 125 * 1e-3;    %[Deg]
Station2.sigma_rho = 0.01;  %[Km]
Station2.min_el = 5;       %[Deg]

%Propagate Mango state
Sat1.xx = zeros(length(et),6);
% First state at visibility window
[Sat1.xx(1,:), ~, ~] = keplerian_propagator(et_ref, x0_mean_sat1, et0 , 'Earth');

for j = 2:length(et)

%State during visibility window
[ Sat1.xx(j,:), ~,~] = keplerian_propagator(et(j-1), Sat1.xx(j-1,:), et(j) , 'Earth');

end

% Compute antenna angles, satellite range and range-rate wrt Station 1
[Station1.sat1.rho, Station1.sat1.azimuth, Station1.sat1.elevation] = pointing(Station1.name,Sat1.xx(:,1:3)',Sat1.xx(:,4:6)',et);

% Compute antenna angles, satellite range and range-rate wrt Station 2
[Station2.sat1.rho, Station2.sat1.azimuth, Station2.sat1.elevation] = pointing(Station2.name,Sat1.xx(:,1:3)',Sat1.xx(:,4:6)',et);


%TLE for Mango
longstr1 = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
longstr2 = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';

sat1rec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec] = invjday(sat1rec.jdsatepoch, sat1rec.jdsatepochf);
sat1_tle_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat1_epoch_et = cspice_str2et(sat1_tle_epoch_str);

% Set nutation corrections parameters (same date as first example -> same values)
ddpsi = -0.111690*arcsec2rad; %  [rad]
ddeps = -0.006855*arcsec2rad; %  [rad]

% Loop over epochs
reci_sat1 = zeros(3,length(et));
veci_sat1 = zeros(3,length(et));

for i = 1:length(et)

    % SGP4 propagation
    tsince = (et(i) - sat1_epoch_et)/60.0;   % time from TLE epoch in minutes
    [~,rteme_sat1,vteme_sat1] = sgp4(sat1rec,  tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(et(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_sat1(:,i), veci_sat1(:,i)] = teme2eci(rteme_sat1, vteme_sat1, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);

end

% Compute antenna angles, satellite range and range-rate wrt Station 1
[Station1.sat1.rho_noise, Station1.sat1.azimuth_noise, Station1.sat1.elevation_noise] = pointing(Station1.name,reci_sat1,veci_sat1,et);






% Compute antenna angles, satellite range and range-rate wrt Station 2
[Station2.sat1.rho_noise, Station2.sat1.azimuth_noise, Station2.sat1.elevation_noise] = pointing(Station2.name,reci_sat1,veci_sat1,et);


% 
% for i = 1:length(et)
% 
%     % Model Noise for station 1
%     Station1.sat1.rho_noise(i) = Station1.sat1.rho_noise(i) + rand(1) * Station1.sigma_rho;
%     Station1.sat1.azimuth_noise(i) = Station1.sat1.azimuth_noise(i)*cspice_dpr() + rand(1) * Station1.sigma_az;
%     Station1.sat1.elevation_noise(i) = Station1.sat1.elevation_noise(i)*cspice_dpr() + rand(1) * Station1.sigma_el;
% 
%     % Model Noise for station 2
%     Station2.sat1.rho_noise(i) = Station2.sat1.rho_noise(i) + rand(1) * Station2.sigma_rho;
%     Station2.sat1.azimuth_noise(i) = Station2.sat1.azimuth_noise(i) *cspice_dpr() + rand(1) * Station2.sigma_az;
%     Station2.sat1.elevation_noise(i) = Station2.sat1.elevation_noise(i) *cspice_dpr() + rand(1) * Station2.sigma_el;
% 
% end

% Model for noise
%Station 1
Station1.sat1.rho_noise = mvnrnd(Station1.sat1.rho,diag(ones(length(et),1)*Station1.sigma_rho));
Station1.sat1.azimuth_noise = mvnrnd(Station1.sat1.azimuth*cspice_dpr(),diag(ones(length(et),1)*Station1.sigma_az));
Station1.sat1.elevation_noise = mvnrnd(Station1.sat1.elevation*cspice_dpr(),diag(ones(length(et),1)*Station1.sigma_el));

% Station 2
Station2.sat1.rho_noise = mvnrnd(Station2.sat1.rho,diag(ones(length(et),1)*Station2.sigma_rho));
Station2.sat1.azimuth_noise = mvnrnd(Station2.sat1.azimuth*cspice_dpr(),diag(ones(length(et),1)*Station2.sigma_az));
Station2.sat1.elevation_noise = mvnrnd(Station2.sat1.elevation*cspice_dpr(),diag(ones(length(et),1)*Station2.sigma_el));

% Visibility windows
i_visibility_station1 = Station1.sat1.elevation * cspice_dpr() > Station1.min_el;
et(i_visibility_station1)

i_visibility_noise_station1 = Station1.sat1.elevation_noise  > Station1.min_el;
et(i_visibility_noise_station1)

i_visibility_station2 = Station2.sat1.elevation * cspice_dpr > Station2.min_el ;
et(i_visibility_station2)
i_visibility_noise_station2 = Station2.sat1.elevation_noise  > Station2.min_el ;
et(i_visibility_noise_station2)

%%
figure()
subplot(1,2,1)
plot(et/cspice_spd(), Station1.sat1.azimuth * cspice_dpr(), 'DisplayName', Sat1.name)
title(['@',Station1.name])
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
legend(Sat1.name)

subplot(1,2,2)
plot(et/cspice_spd(), Station2.sat1.azimuth *cspice_dpr(), 'DisplayName', Sat1.name)
title(['@',Station2.name])
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
legend(Sat1.name)

% Azimuith with noise
figure()
subplot(1,2,1)
plot(et/cspice_spd(), Station1.sat1.azimuth_noise, 'DisplayName', Sat1.name)
title(['@',Station1.name])
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
legend(Sat1.name)

subplot(1,2,2)
plot(et/cspice_spd(), Station2.sat1.azimuth_noise , 'DisplayName', Sat1.name)
title(['@',Station2.name])
xlabel('Epoch [MJD2000]')
ylabel('Azimuth [deg]')
legend(Sat1.name)

% Elevation w/o noise
figure()
subplot(1,2,1)
plot(et/cspice_spd(), Station1.sat1.elevation*cspice_dpr(),'DisplayName', Sat1.name)
title(['@',Station1.name])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend(Sat1.name)

subplot(1,2,2)
plot(et/cspice_spd(), Station2.sat1.elevation*cspice_dpr(),'DisplayName', Sat1.name)
title(['@',Station2.name])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend(Sat1.name)

% Elevation with noise
figure()
subplot(1,2,1)
plot(et/cspice_spd(), Station1.sat1.elevation_noise,'DisplayName', Sat1.name)
title(['@',Station1.name])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend(Sat1.name)

subplot(1,2,2)
plot(et/cspice_spd(), Station2.sat1.elevation_noise,'DisplayName', Sat1.name)
title(['@',Station2.name])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend(Sat1.name)

% Plot passes (azimuth and elevation) w/o noise
figure()
subplot(1,2,1)


scatter(Station1.sat1.azimuth(i_visibility_station1)*cspice_dpr(), Station1.sat1.elevation(i_visibility_station1)*cspice_dpr(),'b','filled','DisplayName', Sat1.name)
hold on;
plot([-180 180],[Station1.min_el Station1.min_el],'--')
axis([-180,180,0, 50])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('','Minimum Elevation');
subplot(1,2,2)


scatter(Station2.sat1.azimuth(i_visibility_station2)*cspice_dpr(), Station2.sat1.elevation(i_visibility_station2)*cspice_dpr(),'b','filled','DisplayName', Sat1.name)
hold on;
plot([-180 180],[Station2.min_el Station2.min_el],'--')
axis([-180,180,0, 50])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('','Minimum Elevation');

% Plot passes (azimuth and elevation) with Noise
figure()
subplot(1,2,1)


scatter(Station1.sat1.azimuth_noise(i_visibility_noise_station1), Station1.sat1.elevation_noise(i_visibility_noise_station1),'b','filled','DisplayName', Sat1.name)
hold on;
plot([-180 180],[Station1.min_el Station1.min_el],'--')
axis([-180,180,0, 50])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('','Minimum Elevation');
subplot(1,2,2)


scatter(Station2.sat1.azimuth_noise(i_visibility_noise_station2), Station2.sat1.elevation_noise(i_visibility_noise_station2),'b','filled','DisplayName', Sat1.name)
hold on;
plot([-180 180],[Station2.min_el Station2.min_el],'--')
axis([-180,180,0, 50])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('','Minimum Elevation');

figure()
skyplot(wrapTo360(Station1.sat1.azimuth(i_visibility_station1)*cspice_dpr()),Station1.sat1.elevation(i_visibility_station1)*cspice_dpr(),MaskElevation=Station1.min_el); 
title(['@',Station1.name])

figure()
skyplot(wrapTo360(Station2.sat1.azimuth(i_visibility_station2)*cspice_dpr()),Station2.sat1.elevation(i_visibility_station2)*cspice_dpr(),MaskElevation=Station2.min_el); 
title(['@',Station2.name])

%With noise
figure()
skyplot(wrapTo360(Station1.sat1.azimuth_noise(i_visibility_noise_station1)),Station1.sat1.elevation_noise(i_visibility_noise_station1),MaskElevation=Station1.min_el); 
title(['@',Station1.name])

figure()
skyplot(wrapTo360(Station2.sat1.azimuth_noise(i_visibility_noise_station2)),Station2.sat1.elevation_noise(i_visibility_noise_station2),MaskElevation=Station2.min_el); 
title(['@',Station2.name])

%% Functions
%
%
function [xf, tt, xx] = keplerian_propagator(t0, x0, t1 , attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 


    GM = cspice_bodvrd(attractor, 'GM', 1);


options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 t1-t0], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end
%
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
%
%
%

function [rho, azimuth, elevation] = pointing(stationName,rr_ECI,vv_ECI,et)
%
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

%
%
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
%
%
%