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
 P0 = [   [5.6e-7  3.5e-7  -7.1e-8 
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


% First state 
[Sat1.xx(1,:), ~, ~] = keplerian_propagator(et_ref, Sat1.x0_mean, et0 , 'Earth');

for j = 2:length(et)

%State 
[ Sat1.xx(j,:), ~,~] = keplerian_propagator(et(j-1), Sat1.xx(j-1,:), et(j) , 'Earth');

end


% Compute antenna angles, satellite range wrt Station 2
[Station2.sat1.rho, Station2.sat1.azimuth, Station2.sat1.elevation] = pointing(Station2.name,Sat1.xx(:,1:3)',Sat1.xx(:,4:6)',et);


%TLE for Mango
longstr1 = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
longstr2 = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';

sat1rec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec] = invjday(sat1rec.jdsatepoch, sat1rec.jdsatepochf);
sat1_tle_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat1_epoch_et = cspice_str2et(sat1_tle_epoch_str);

% Set nutation corrections parameters 
ddpsi = -0.073296 * arcsec2rad; %  [rad]
ddeps = -0.009373 * arcsec2rad; %  [rad]

% Loop over epochs
Sat1.reci = zeros(3,length(et));
Sat1.veci = zeros(3,length(et));

for i = 1:length(et)

    % SGP4 propagation
    tsince = (et(i) - sat1_epoch_et)/60.0;   % time from TLE epoch in minutes
    [~,Sat1.rteme,Sat1.vteme] = sgp4(sat1rec,  tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(et(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [Sat1.reci(:,i), Sat1.veci(:,i)] = teme2eci(Sat1.rteme ,Sat1.vteme, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);

end

% Compute antenna angles, satellite range and range-rate wrt Station 2
[Station2.sat1.rho, Station2.sat1.azimuth, Station2.sat1.elevation] = pointing(Station2.name,Sat1.reci,Sat1.veci,et);

% Station 2 with nosie
i_visibility_noise_station2 = Station2.sat1.elevation*cspice_dpr()  > Station2.min_el ;
Station2.sat1.visibility_noise = et(i_visibility_noise_station2);
visibility(Station2.sat1.visibility_noise,Station2.name,Sat1.name,5)

% Model for noise
% Station 2
Station2.sat1.rho_noise = mvnrnd(Station2.sat1.rho(i_visibility_noise_station2),ones(1,length(Station2.sat1.visibility_noise))*Station2.sigma_rho^2);
Station2.sat1.azimuth_noise = mvnrnd(Station2.sat1.azimuth(i_visibility_noise_station2)*cspice_dpr(),ones(1,length(Station2.sat1.visibility_noise))*Station2.sigma_az^2);
Station2.sat1.elevation_noise = mvnrnd(Station2.sat1.elevation(i_visibility_noise_station2)*cspice_dpr(),ones(1,length(Station2.sat1.visibility_noise))*Station2.sigma_el^2);

y_meas = [Station2.sat1.rho_noise;Station2.sat1.azimuth_noise;Station2.sat1.elevation_noise];

% Plot passes (azimuth and elevation) w noise
figure()
scatter(Station2.sat1.azimuth_noise, Station2.sat1.elevation_noise,'b','filled','DisplayName', Sat1.name)
hold on;
plot([-180 180],[Station2.min_el Station2.min_el],'--')
axis([-180,180,0, 50])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('','Minimum Elevation');


%% 
% Unscented Transform
alpha = 0.1;
n = 6;
lambda = alpha^2 * n;
beta = 2;
R = diag([Station2.sigma_rho^2,Station2.sigma_az^2 , Station2.sigma_el^2]);
time = [et0 ,Station2.sat1.visibility_noise];

% Initialize Sigma points vector
Sat1.x_sigma = zeros(6,13);
  
% Initialize vector containing last propagation of sigma points 
Sat1.xx_sigma = zeros(6,13);

% Initialization vector of mean values 
Sat1.xk_m = zeros(n,length(time)-1); 
Sat1.yk_m = zeros(3,length(time)-1);
Sat1.yyk = zeros(3,2*n+1);
Sat1.xk_up = zeros(6,length(time));


% Initialization vector of covariance matrices for every step
Sat1.Pk_m = zeros(n,n,length(time)-1);
Sat1.Pyy = zeros(3,3,length(time)-1);
Sat1.Pxy = zeros(6,3,length(time)-1);


% Weights for Mean
W_m = zeros(1,2*n + 1);
W_m(1) =  1 - n/(alpha^2*n);
W_m(2:end) = 1/(2 * alpha^2 * n);

% Weights for Covariance
W_c = zeros(1,2*n + 1);
W_c(1) = (2 - alpha^2 + beta) - n/(alpha^2 * n);
W_c(2:end) = 1/(2 * alpha^2 * n);


% Initialize first iteration values
 Sat1.Pk_up(:,:,1) =  1e4 * P0;
 Sat1.xk_up(:,1) = [Sat1.reci(:,1);Sat1.veci(:,1)];

% First iteration
Sat1.P_sq1 = sqrtm(lambda * Sat1.Pk_up(:,:,1));


for j = 2 : length(time)

    for i = 1 : 2 * n + 1

        % Compute Sigma Points
         if(i==1)

            Sat1.x_sigma(:,i) = Sat1.xk_up(:,j-1) ;

         elseif(i>1 && i<8)

            Sat1.x_sigma(:,i) = Sat1.xk_up(:,j-1) + Sat1.P_sq1(:,i-1);
            
         else
    
            Sat1.x_sigma(:,i) = Sat1.xk_up(:,j-1) - Sat1.P_sq1(:,i-n-1);

         end

        % Propagate Sigma points
        [Sat1.xx_sigma(:,i), ~, ~ ] = keplerian_propagator_J2(time(j-1),Sat1.x_sigma(:,i), time(j), 'Earth');
        
        % Sigma points in Measurement space
        [Station2.sat1.rho_meas(i), Station2.sat1.azimuth_meas(i), Station2.sat1.elevation_meas(i)] = pointing(Station2.name,Sat1.xx_sigma(1:3,i),Sat1.xx_sigma(4:6,i),time(j));
       
        % Sample Mean
        Sat1.xk_m(:,j-1) = Sat1.xk_m(:,j-1) + W_m(i) * Sat1.xx_sigma(:,i); 

        Sat1.yyk(:,i) = [Station2.sat1.rho_meas(i);  Station2.sat1.azimuth_meas(i); Station2.sat1.elevation_meas(i)];   

       % Sample Mean of measurements
        Sat1.yk_m(:,j-1) = Sat1.yk_m(:,j-1) + W_m(i) * Sat1.yyk(:,i); 
        
    end

   
    for i = 1: 2*n+1

    Sat1.Pk_m(:,:,j-1) = Sat1.Pk_m(:,:,j-1) + W_c(i) * ((Sat1.xx_sigma(:,i) - Sat1.xk_m(:,j-1)) * (Sat1.xx_sigma(:,i) - Sat1.xk_m(:,j-1))');
    Sat1.Pyy(:,:,j-1) = Sat1.Pyy(:,:,j-1) + W_c(i) * ((Sat1.yyk(:,i) - Sat1.yk_m(:,j-1)) * (Sat1.yyk(:,i) - Sat1.yk_m(:,j-1))');
    Sat1.Pxy(:,:,j-1) = Sat1.Pxy(:,:,j-1) + W_c(i) * ((Sat1.xx_sigma(:,i) - Sat1.xk_m(:,j-1)) * (Sat1.yyk(:,i) - Sat1.yk_m(:,j-1))');

    end
    Sat1.Pyy(:,:,j-1) = Sat1.Pyy(:,:,j-1) + R;

    % Kalman Gain
    Kk = Sat1.Pxy(:,:,j-1) * inv(Sat1.Pyy(:,:,j-1));
    % Update State
    Sat1.xk_up(:,j) = Sat1.xk_m(:,j-1) + Kk * (y_meas(:,j-1) - Sat1.yk_m(:,j-1));
    Sat1.Pk_up(:,:,j) = Sat1.Pk_m(:,:,j-1) - Kk * Sat1.Pyy(:,:,j-1) * Kk';
    Sat1.Pk_up(:,:,j) = (Sat1.Pk_up(:,:,j) + Sat1.Pk_up(:,:,j)')/2;
    Sat1.P_sq1 = sqrtm(lambda * Sat1.Pk_up(:,:,j));
    Sat1.sigma_pos(j-1) = 3 * sqrt(max(eig( Sat1.Pk_up(1:3,1:3,j))));
    Sat1.sigma_vel(j-1) = 3 * sqrt(max(eig( Sat1.Pk_up(4:6,4:6,j)))); 

    
  
 
end

State = [Sat1.reci(:,i_visibility_noise_station2); Sat1.veci(:,i_visibility_noise_station2)];

Sat1.err_pos = Sat1.xk_up(1:3,2:end) - State(1:3,:);
Sat1.norm_err_pos = vecnorm(Sat1.err_pos);

Sat1.err_vel = Sat1.xk_up(4:6,2:end) - State(4:6,:);
Sat1.norm_err_vel = vecnorm(Sat1.err_vel);


figure()
subplot(2,1,1)
plot(Station2.sat1.visibility_noise,Sat1.norm_err_pos,'b','LineWidth',2);
hold on;
plot(Station2.sat1.visibility_noise,Sat1.sigma_pos,'--r','LineWidth',2);
legend('Position error','3 $\sigma$','Interpreter','latex');



subplot(2,1,2)
plot(Station2.sat1.visibility_noise,Sat1.norm_err_vel,'b','LineWidth',2);
hold on;
plot(Station2.sat1.visibility_noise,Sat1.sigma_vel,'--r','LineWidth',2);
legend('Velocity error','3 $\sigma$','Interpreter','latex');




%% Functions
function visibility(et,stationname,satname,tstep)
%
% This function evaluate the boundaries of the visibility windows given as
% input the visibility times with a fixed step
%
%   INPUT:      -et: Time of visibility in ET   [s], [1,n]
%               -tstep: Given step of et, (1 minute if not specified)
%
%
%   AUTHOT: Matteo Bono
%
%   
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
    end
k = k+1;
end
end
%
%
function [xf, tt, xx] = keplerian_propagator(t0, x0, t1 , attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 


    GM = cspice_bodvrd(attractor, 'GM', 1);


options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-9; ones(3,1)*1e-12]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 t1-t0], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end
%
%
function [xf, tt, xx] = keplerian_propagator_J2(t0, x0, t1 , attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 


    GM = cspice_bodvrd(attractor, 'GM', 1);


options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-9; ones(3,1)*1e-12]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs_J2(t,x,GM), [0 t1-t0], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end
%
%
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
%
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
%
