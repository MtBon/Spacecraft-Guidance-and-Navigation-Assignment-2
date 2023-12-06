% Spacecraft Guidance and Navigation (2023/2024)
% Assigment #2
% Exercise #1
% Author: Matteo Bono
%
%% Point 1

clc; clearvars; cspice_kclear; close all

% Replace the path below with your own installation path if not present as
% a subfolder
addpath('functions');
addpath('kernels');
matlab_graphics

% Load spice kernels
cspice_furnsh('assignment02.tm');


% Chaser Satellite
Sat1_name = 'Mango';
Sat1_ID = 36599;

% Target Satellite
Sat2_name = 'Tango';
Sat2_ID = 36827;

% Separation Epoch
et0 = cspice_str2et('2010-08-12T05:27:29.114');   %ET

% Number of orbital periods
N = 10;

% Mean Initial State for Satellite 1 (Mango)
r0_mean_sat1 = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v0_mean_sat1 = [0.812221125483763; -0.721512914578826; 7.42665302729053];
x0_mean_sat1 = [r0_mean_sat1; v0_mean_sat1];

% Mean Initial State for Satellite 2 (Tango)
r0_mean_sat2 = [4621.69343340281; 5399.26386352847; -3.09039248714313];
v0_mean_sat2 = [0.813960847513811; -0.719449862738607; 7.42706066911294];
x0_mean_sat2 = [r0_mean_sat2; v0_mean_sat2];

% Covariance Matrix
 P0 = [   [5.6e-7  3.5e-7  -7.1e-8 
           3.5e-7  9.7e-7   7.6e-8      
          -7.1e-8  7.6e-8   8.1e-8]            zeros(3)

           zeros(3)        diag([2.8e-11 , 2.7e-11 , 9.6e-12])];


% Orbital parameters for Satellite 1
rr0_sat1 = norm(r0_mean_sat1);
vv0_sat1 = norm(v0_mean_sat1);
mu = cspice_bodvrd('Earth','GM',1);

% Semi major axis
a = 1/((2/rr0_sat1) - (vv0_sat1^2/mu));

% eccentricity
e = norm(1/mu * ((vv0_sat1^2 - (mu/rr0_sat1)) * r0_mean_sat1 - dot(r0_mean_sat1,v0_mean_sat1) * v0_mean_sat1));

[a1, e1, i1, Om1, w1, theta1] = rv2kp(r0_mean_sat1, v0_mean_sat1, mu);
[a2, e2, i2, Om2, w2, theta2] = rv2kp(r0_mean_sat2, v0_mean_sat2, mu);
% Orbital Period
T1 = 2 * pi * sqrt(a^3/mu);

% Final epoch after N orbital periods
etf = et0 +  N * T1;

% Time span
et = et0 : T1 : etf;

%Initialization
xx_sat1 = zeros(length(et),6);
xx_sat2 = zeros(length(et),6);
PHI_sat1 = zeros(6,6,length(et));
PHI_sat2 = zeros(6,6,length(et));


% First component of Vectors of Propagated states for N revolutions
xx_sat1(1,:) = x0_mean_sat1';
xx_sat2(1,:) = x0_mean_sat2';

% First STM Mastrices
PHI_sat1(:,:,1) = eye(6);
PHI_sat2(:,:,1) = eye(6);

% First state
x0_sat1 = x0_mean_sat1;
x0_sat2 = x0_mean_sat2;

% Propagation with Linearized Approach (LinCov)
for i = 1 : length(et)-1

% Satellite 1
[xx_sat1(i+1,:), PHI_sat1(:,:,i+1), ~, ~ ] = keplerian_propagator_STM(et(i),x0_sat1, et(i+1) , 'Earth');

% Satellite 2
[xx_sat2(i+1,:), PHI_sat2(:,:,i+1), tt, ~ ] = keplerian_propagator_STM(et(i),x0_sat2, et(i+1) , 'Earth');

% Update State for New Iteration
x0_sat1 = xx_sat1(i+1,:)';
x0_sat2 = xx_sat2(i+1,:)';

end
%%
% Retrieve Propagated States
rr_sat1 = xx_sat1(:,1:3);
vv_sat1 = xx_sat1(:,4:6);

rr_sat2 = xx_sat2(:,1:3);
vv_sat2 = xx_sat2(:,4:6);


% Covariance Matrices
P_sat1 = zeros(6,6,length(et));
P_sat2 = zeros(6,6,length(et));
P0_sat1 = P0;
P0_sat2 = P0;

for i = 1:length(et)

    % Compute new Covariance matrices for every Revolution

    % Sat1
    P_sat1(:,:,i) = PHI_sat1(:,:,i) * P0_sat1 * PHI_sat1(:,:,i)';
    P0_sat1 = P_sat1(:,:,i);

    % Sat 2
    P_sat2(:,:,i) = PHI_sat2(:,:,i) * P0_sat2 * PHI_sat2(:,:,i)';
    P0_sat2 = P_sat2(:,:,i);
end



% Unscented Transform
alpha = 0.1;
n = 6;
lambda = alpha^2 * n;
beta = 2;


% Initialize Sigma points vector
  x_sigma_sat1 = zeros(6,13);
  x_sigma_sat2 = zeros(6,13);

% Initialize vector containing last propagation of sigma points 
xx_sigma_sat1 = zeros(6,13);
xx_sigma_sat2 = zeros(6,13);

% Initialization vector of mean values for every N revolution
sigma_mean_sat1 = zeros(n,length(et)); 
sigma_mean_sat2 = zeros(n,length(et)); 

% Initialization vector of covariance matrices for every step
sigma_cov_sat1 = zeros(n,n,length(et));
sigma_cov_sat2 = zeros(n,n,length(et));


% Weights for Mean
W_m = zeros(1,2*n + 1);
W_m(1) =  1 - n/(alpha^2*n);
W_m(2:end) = 1/(2 * alpha^2 * n);

% Weights for Covariance
W_c = zeros(1,2*n + 1);
W_c(1) = (2 - alpha^2 + beta) - n/(alpha^2 * n);
W_c(2:end) = 1/(2 * alpha^2 * n);


% Initialize first iteration values
 sigma_cov_sat1(:,:,1) = P0;
 sigma_cov_sat2(:,:,1) = P0;
 sigma_mean_sat1(:,1) = x0_mean_sat1;
 sigma_mean_sat2(:,1) = x0_mean_sat2;

% First iteration
P_sq1 = sqrtm(lambda * sigma_cov_sat1(:,:,1)) ;
P_sq2 = sqrtm(lambda * sigma_cov_sat2(:,:,1)) ;

for j = 2 : length(et)

    for i = 1 : 2 * n + 1

        % Compute Sigma Points
         if(i==1)

            x_sigma_sat1(:,i) = sigma_mean_sat1(:,j-1) ;
            x_sigma_sat2(:,i) = sigma_mean_sat2(:,j-1) ;

         elseif(i>1 && i<8)

            x_sigma_sat1(:,i) = sigma_mean_sat1(:,j-1) + P_sq1(:,i-1);
            x_sigma_sat2(:,i) = sigma_mean_sat2(:,j-1) + P_sq2(:,i-1);

         else
    
            x_sigma_sat2(:,i) = sigma_mean_sat2(:,j-1) - P_sq2(:,i-n-1);
            x_sigma_sat1(:,i) = sigma_mean_sat1(:,j-1) - P_sq1(:,i-n-1);
         end

        % Propagate Sigma points
        [xx_sigma_sat1(:,i), ~, ~ ] = keplerian_propagator(et(j-1),x_sigma_sat1(:,i), et(j), 'Earth');
    
        [xx_sigma_sat2(:,i), ~, ~ ] = keplerian_propagator(et(j-1),x_sigma_sat2(:,i), et(j), 'Earth');

        % Sample Mean
        sigma_mean_sat1(:,j) = sigma_mean_sat1(:,j) + W_m(i) * xx_sigma_sat1(:,i); 
        sigma_mean_sat2(:,j) = sigma_mean_sat2(:,j) + W_m(i) * xx_sigma_sat2(:,i); 

    end

   
    for i = 1: 2*n+1

    sigma_cov_sat1(:,:,j) = sigma_cov_sat1(:,:,j) + W_c(i) * ((xx_sigma_sat1(:,i) - sigma_mean_sat1(:,j)) * (xx_sigma_sat1(:,i) - sigma_mean_sat1(:,j))');
    sigma_cov_sat2(:,:,j) = sigma_cov_sat2(:,:,j) + W_c(i) * ((xx_sigma_sat2(:,i) - sigma_mean_sat2(:,j)) * (xx_sigma_sat2(:,i) - sigma_mean_sat2(:,j))');

    end

    % Update Variance matrices
    P_sq1 = sqrtm(lambda * sigma_cov_sat1(:,:,j)) ;
    P_sq2 = sqrtm(lambda * sigma_cov_sat2(:,:,j)) ;

end
            


%% Point 2

% Compute the norm of the relative position

% Initializations
Dr_Lin = zeros(length(et),1);
Dr_UT = zeros(length(et),1);
P_pos_UT_sat1 = zeros(3,3,length(et));
P_pos_UT_sat2 = zeros(3,3,length(et));
P_vel_UT_sat1 = zeros(3,3,length(et));
P_vel_UT_sat2 = zeros(3,3,length(et));
P_pos_sum_UT = zeros(3,3,length(et));
cov_sum_pos_UT = zeros(1,length(et));
cov_pos_sat1_UT = zeros(1,length(et));
cov_vel_sat1_UT = zeros(1,length(et));
cov_pos_sat2_UT = zeros(1,length(et));
cov_vel_sat2_UT = zeros(1,length(et));

P_pos_Lin_sat1 = zeros(3,3,length(et));
P_pos_Lin_sat2 = zeros(3,3,length(et));
P_vel_Lin_sat1 = zeros(3,3,length(et));
P_vel_Lin_sat2 = zeros(3,3,length(et));
P_pos_sum_Lin = zeros(3,3,length(et));
cov_sum_pos_Lin = zeros(1,length(et));
cov_pos_sat1_Lin = zeros(1,length(et));
cov_pos_sat2_Lin = zeros(1,length(et));
cov_vel_sat1_Lin = zeros(1,length(et));
cov_vel_sat2_Lin = zeros(1,length(et));

for i = 1 : length(et)

    % Extract Covariance matrix on position for UT method
    % Sub-matrices of Covariance for Position
    P_pos_UT_sat1(:,:,i) = sigma_cov_sat1(1:3,1:3,i);
    P_pos_UT_sat2(:,:,i) = sigma_cov_sat2(1:3,1:3,i);

    % Sub-matrices of Covariance for Velocity
    P_vel_UT_sat1(:,:,i) = sigma_cov_sat1(4:6,4:6,i);
    P_vel_UT_sat2(:,:,i) = sigma_cov_sat2(4:6,4:6,i);

    % Sum
    P_pos_sum_UT(:,:,i) = P_pos_UT_sat1(:,:,i) + P_pos_UT_sat2(:,:,i);
  

    % Extract P_sum_Lin Covariance matrix on position for LinCov method
    % Sub-matrices of Covariance for Position
    P_pos_Lin_sat1(:,:,i) = P_sat1(1:3,1:3,i);
    P_pos_Lin_sat2(:,:,i) = P_sat2(1:3,1:3,i);

    % Sub-matrices of Covariance for Velocity
    P_vel_Lin_sat1(:,:,i) = P_sat1(4:6,4:6,i);
    P_vel_Lin_sat2(:,:,i) = P_sat2(4:6,4:6,i);

    %Sum
    P_pos_sum_Lin(:,:,i) = P_pos_Lin_sat1(:,:,i) + P_pos_Lin_sat2(:,:,i);


    % Critical condition
    cov_sum_pos_Lin(i) = 3 * sqrt(max(eig(P_pos_sum_Lin(:,:,i))));
    cov_sum_pos_UT(i) = 3 * sqrt(max(eig(P_pos_sum_UT(:,:,i))));

    %Values to plot for Point 3
    %LinCov
    cov_pos_sat1_Lin(i) = 3 * sqrt(max(eig(P_pos_Lin_sat1(:,:,i))));
    cov_vel_sat1_Lin(i) = 3 * sqrt(max(eig(P_vel_Lin_sat1(:,:,i))));
    cov_pos_sat2_Lin(i) = 3 * sqrt(max(eig(P_pos_Lin_sat2(:,:,i))));
    cov_vel_sat2_Lin(i) = 3 * sqrt(max(eig(P_vel_Lin_sat2(:,:,i))));
   
    % UT
    cov_pos_sat1_UT(i) = 3 * sqrt(max(eig(P_pos_UT_sat1(:,:,i))));
    cov_vel_sat1_UT(i) = 3 * sqrt(max(eig(P_vel_UT_sat1(:,:,i))));
    cov_pos_sat2_UT(i) = 3 * sqrt(max(eig(P_pos_UT_sat2(:,:,i))));
    cov_vel_sat2_UT(i) = 3 * sqrt(max(eig(P_vel_UT_sat2(:,:,i))));
    
    
    Dr_Lin(i) = norm(rr_sat1(i,:) - rr_sat2(i,:));
    Dr_UT(i) = norm(sigma_mean_sat1(1:3,i) - sigma_mean_sat2(1:3,i));

    if(Dr_UT(i) < cov_sum_pos_UT(i))
        rev = i-1;
        fprintf('UT : Critical Condition occurs at Revolution n.%d\n',rev)
    end

    if(Dr_Lin(i) < cov_sum_pos_Lin(i))
        rev = i-1;
        fprintf('LinCov : Critical Condition occurs at Revolution n.%d\n',rev)
        fprintf('----------------------------------------------------\n');
    end
end

%% Point 3

num_sim  = 200;


% Initializations
mc_sat1 = zeros(num_sim,n,length(et));
mc_sat2 = zeros(num_sim,n,length(et));
xx_mc_sat1 = zeros(n,num_sim);
xx_mc_sat2 = zeros(n,num_sim);
mc_mean_sat1 = zeros(length(et),n);
mc_mean_sat2 = zeros(length(et),n);
cov_pos_sat1_mc = zeros(1,length(et));
cov_vel_sat1_mc = zeros(1,length(et));
cov_pos_sat2_mc = zeros(1,length(et));
cov_vel_sat2_mc = zeros(1,length(et));

% Sample Covariance
mc_cov_sat1 = zeros(6,6,length(et));
mc_cov_sat2 = zeros(6,6,length(et));


mc_mean_sat1(1,:) = x0_mean_sat1';
mc_mean_sat2(1,:) = x0_mean_sat2';
mc_cov_sat1(:,:,1) = P0;
mc_cov_sat2(:,:,1) = P0;

mc_sat1(:,:,1) = mvnrnd(mc_mean_sat1(1,:),mc_cov_sat1(:,:,1),num_sim);
mc_sat2(:,:,1) = mvnrnd(mc_mean_sat2(1,:),mc_cov_sat2(:,:,1),num_sim);

fprintf('Monte-Carlo simulation Started...\n');
for j = 2 : length(et)

    % Generate points
    

    for i = 1 : num_sim

     % Perform Propagation
    [xx_mc_sat1(:,i),~, ~] = keplerian_propagator(et(j-1), mc_sat1(i,:,j-1), et(j), 'Earth');

    [xx_mc_sat2(:,i), tt, ~] = keplerian_propagator(et(j-1),mc_sat2(i,:,j-1), et(j), 'Earth');

    mc_mean_sat1(j,:) = mc_mean_sat1(j,:) + xx_mc_sat1(:,i)';
    mc_mean_sat2(j,:) = mc_mean_sat2(j,:) + xx_mc_sat2(:,i)';

    end
% 
%      mc_mean_sat1(j,:) = 1/num_sim * mc_mean_sat1(j,:) ;
%      mc_mean_sat2(j,:) = 1/num_sim * mc_mean_sat2(j,:) ;

        % Sample mean
        mc_mean_sat1(j,:) = mean(xx_mc_sat1,2);
        mc_mean_sat2(j,:) = mean(xx_mc_sat2,2);

%         for i = 1 : num_sim
%         % Sample Covariance
%         mc_cov_sat1(:,:,j) = mc_cov_sat1(:,:,j) + ((xx_mc_sat1(:,i) - mc_mean_sat1(j,:)') * (xx_mc_sat1(:,i)' - mc_mean_sat1(j,:)));
%         mc_cov_sat2(:,:,j) = mc_cov_sat2(:,:,j) + ((xx_mc_sat2(:,i) - mc_mean_sat2(j,:)') * (xx_mc_sat2(:,i)' - mc_mean_sat2(j,:))); 
% 
%         end
%         mc_cov_sat1(:,:,j) = sqrtm(1/(num_sim-1) * mc_cov_sat1(:,:,j));
%         mc_cov_sat2(:,:,j) = sqrtm(1/(num_sim-1) * mc_cov_sat2(:,:,j));
          mc_cov_sat1(:,:,j) = cov(xx_mc_sat1');
          mc_cov_sat2(:,:,j) = cov(xx_mc_sat2');

          % Generate points distribution at this time step
          mc_sat1(:,:,j) = mvnrnd(mc_mean_sat1(j,:),mc_cov_sat1(:,:,j),num_sim);
          mc_sat2(:,:,j) = mvnrnd(mc_mean_sat2(j,:),mc_cov_sat2(:,:,j),num_sim);


end

fprintf('Finished\n');



for i = 1 : length(et)
    % Values to plot
    cov_pos_sat1_mc(i) = 3 * sqrt(max(eig(mc_cov_sat1(1:3,1:3,i))));
    cov_vel_sat1_mc(i) = 3 * sqrt(max(eig(mc_cov_sat1(4:6,4:6,i))));

    cov_pos_sat2_mc(i) = 3 * sqrt(max(eig(mc_cov_sat2(1:3,1:3,i))));
    cov_vel_sat2_mc(i) = 3 * sqrt(max(eig(mc_cov_sat2(4:6,4:6,i))));

end





%% Plots

figure(1)
plot3(rr_sat1(:,1), rr_sat1(:,2), rr_sat1(:,3),'.','MarkerSize',30);
grid on;
xlabel('X[Km]');
ylabel('Y[Km]');    
zlabel('Z[Km]');
title('Mean State Propagated with LinCov');

% plot_gaussian_ellipsoid(x0_mean_sat1(1:3), P0(1:3,1:3))
% hold on;
% plot_gaussian_ellipsoid(xf_sat1(1:3), P_f(1:3,1:3))


figure(2)
plot3(rr_sat1(end,1), rr_sat1(end,2), rr_sat1(end,3),'.','MarkerSize',40);
hold on;
plot3(sigma_mean_sat1(1,end),sigma_mean_sat1(2,end),sigma_mean_sat1(3,end),'.k','MarkerSize',40);
%plot3(xf_sigma(end,1),xf_sigma(end,2),xf_sigma(end,3),'r.','MarkerSize',15);
grid on;
xlabel('X[Km]');
ylabel('Y[Km]');    
zlabel('Z[Km]');
legend('LinCov','UT');
title('Mean State Propagated');

tp = linspace(1,11,11);

figure(3)
plot(tp,Dr_UT,tp,Dr_Lin,'LineWidth',2);
xlabel('# Orbital Periods');
ylabel('[Km]')
legend('$\Delta r_{Lin}$','$\Delta r_{UT}$','interpreter','latex')

% Plots for Point 3
figure(4)
subplot(1,2,1)
plot(tp,cov_pos_sat1_Lin,tp,cov_pos_sat1_UT,tp,cov_pos_sat1_mc,'LineWidth',2)
xlabel('#Orbital Period');
title('Position Covariance for Satellite Mango')
legend('LinCov','UT','MC');

subplot(1,2,2)
plot(tp,cov_pos_sat2_Lin,tp,cov_pos_sat2_UT,tp,cov_pos_sat2_mc,'LineWidth',2)
xlabel('#Orbital Period');
title('Position Covariance for Satellite Tango')
legend('LinCov','UT','MC');

figure(5)
subplot(1,2,1)
plot(tp,cov_vel_sat1_Lin,tp,cov_vel_sat1_UT,tp,cov_vel_sat1_mc,'LineWidth',2)
xlabel('#Orbital Period');
title('Velocity Covariance for Satellite Mango')
legend('LinCov','UT','MC');

subplot(1,2,2)
plot(tp,cov_vel_sat2_Lin,tp,cov_vel_sat2_UT,tp,cov_vel_sat2_mc,'LineWidth',2)
xlabel('#Orbital Period');
title('Velocity Covariance for Satellite Tango')
legend('LinCov','UT','MC');


figure(6)
% DisegnaOrbita(a1,e1,w1,i1,Om1);
hold on;
    scatter3(mc_sat1(:,1,1),mc_sat1(:,2,1),mc_sat1(:,3,1),14);
    scatter3(mc_sat1(:,1,5),mc_sat1(:,2,5),mc_sat1(:,3,5),14);


%% Functions

% Keplerian Propagator
function [xf, PHI_f, tt, xx ] = keplerian_propagator_STM( et0,x0, et1 , attractor)

%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 

GM = cspice_bodvrd(attractor, 'GM', 1);

options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12*ones(42,1));


x0Phi0 = [ x0 ; reshape(eye(6),36,1) ];

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_STM_rhs(t,x,GM),[0 et1-et0], x0Phi0, options_STM);

% Extract state vector 
xf = xx(end,1:6);
PHI_f = reshape(xx(end,7:end),6,6);

end
%
%
%
function [dxdt] = keplerian_STM_rhs(~, x, GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian)
%               propagator with STM
%   Evaluates the right-hand-side of a newtonian 2-body propagator with STM.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 11/10/2023
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2023 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2023/2024.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [42, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [42,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(42,1);

% Extract positions
rr = x(1:3);
Phi = reshape(x(7:end),6,6);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Compute the gravitational acceleration using Newton's law
aa_grav =  - GM * rr /(dist*dist2);

% Compute the derivative of the flow
dfdv = 3 * GM / dist^5 * (rr * rr') - GM / dist^3 * eye(3);

% Assemble the matrix A(t)=dfdx
dfdx = [zeros(3),  eye(3);
        dfdv    ,  zeros(3)];
% Compute the derivative of the state transition matrix
Phidot = dfdx * Phi;


dxdt(1:3) = x(4:6);   % Position detivative is object's velocity
dxdt(4:6) = aa_grav;  % Sum up acceleration to right-hand-side
dxdt(7:end) = Phidot(:);

end
%
%
%
function h = plot_gaussian_ellipsoid(m, C, sdwidth, npts, axh)
% PLOT_GAUSSIAN_ELLIPSOIDS plots 2-d and 3-d Gaussian distributions
%
% H = PLOT_GAUSSIAN_ELLIPSOIDS(M, C) plots the distribution specified by 
%  mean M and covariance C. The distribution is plotted as an ellipse (in 
%  2-d) or an ellipsoid (in 3-d).  By default, the distributions are 
%  plotted in the current axes. H is the graphics handle to the plotted 
%  ellipse or ellipsoid.
%
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD) uses SD as the standard deviation 
%  along the major and minor axes (larger SD => larger ellipse). By 
%  default, SD = 1. Note: 
%  * For 2-d distributions, SD=1.0 and SD=2.0 cover ~ 39% and 86% 
%     of the total probability mass, respectively. 
%  * For 3-d distributions, SD=1.0 and SD=2.0 cover ~ 19% and 73%
%     of the total probability mass, respectively.
%  
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD, NPTS) plots the ellipse or 
%  ellipsoid with a resolution of NPTS (ellipsoids are generated 
%  on an NPTS x NPTS mesh; see SPHERE for more details). By
%  default, NPTS = 50 for ellipses, and 20 for ellipsoids.
%
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD, NPTS, AX) adds the plot to the
%  axes specified by the axis handle AX.
%
% Examples: 
% -------------------------------------------
%  % Plot three 2-d Gaussians
%  figure; 
%  h1 = plot_gaussian_ellipsoid([1 1], [1 0.5; 0.5 1]);
%  h2 = plot_gaussian_ellipsoid([2 1.5], [1 -0.7; -0.7 1]);
%  h3 = plot_gaussian_ellipsoid([0 0], [1 0; 0 1]);
%  set(h2,'color','r'); 
%  set(h3,'color','g');
% 
%  % "Contour map" of a 2-d Gaussian
%  figure;
%  for sd = [0.3:0.4:4],
%    h = plot_gaussian_ellipsoid([0 0], [1 0.8; 0.8 1], sd);
%  end
%
%  % Plot three 3-d Gaussians
%  figure;
%  h1 = plot_gaussian_ellipsoid([1 1  0], [1 0.5 0.2; 0.5 1 0.4; 0.2 0.4 1]);
%  h2 = plot_gaussian_ellipsoid([1.5 1 .5], [1 -0.7 0.6; -0.7 1 0; 0.6 0 1]);
%  h3 = plot_gaussian_ellipsoid([1 2 2], [0.5 0 0; 0 0.5 0; 0 0 0.5]);
%  set(h2,'facealpha',0.6);
%  view(129,36); set(gca,'proj','perspective'); grid on; 
%  grid on; axis equal; axis tight;
% -------------------------------------------
% 
%  Gautam Vallabha, Sep-23-2007, Gautam.Vallabha@mathworks.com
%  Revision 1.0, Sep-23-2007
%    - File created
%  Revision 1.1, 26-Sep-2007
%    - NARGOUT==0 check added.
%    - Help added on NPTS for ellipsoids
if ~exist('sdwidth', 'var'), sdwidth = 1; end
if ~exist('npts', 'var'), npts = []; end
if ~exist('axh', 'var'), axh = gca; end
if numel(m) ~= length(m) 
    error('M must be a vector'); 
end
if ~( all(numel(m) == size(C)) )
    error('Dimensionality of M and C must match');
end
if ~(isscalar(axh) && ishandle(axh) && strcmp(get(axh,'type'), 'axes'))
    error('Invalid axes handle');
end
set(axh, 'nextplot', 'add');
switch numel(m)
   case 2, h=show2d(m(:),C,sdwidth,npts,axh);
   case 3, h=show3d(m(:),C,sdwidth,npts,axh);
   otherwise
      error('Unsupported dimensionality');
end
if nargout==0
    clear h;
end
end

%-----------------------------
function h = show2d(means, C, sdwidth, npts, axh)
if isempty(npts), npts=50; end
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
h = plot(bp(1,:), bp(2,:), '-', 'parent', axh);
end

%-----------------------------
function h = show3d(means, C, sdwidth, npts, axh)
if isempty(npts), npts=20; end
[x,y,z] = sphere(npts);
ap = [x(:) y(:) z(:)]';
[v,d]=eig(C); 
if any(d(:) < 0)
   fprintf('warning: negative eigenvalues\n');
   d = max(d,0);
end

d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
xp = reshape(bp(1,:), size(x));
yp = reshape(bp(2,:), size(y));
zp = reshape(bp(3,:), size(z));
h = surf(axh, xp,yp,zp);
end
%
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
%%
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
function [a, e, i, OM, om, theta] = rv2kp(r, v, mu)

%  Conversion from Keplerian elements to Cartesian coordinates. Angles in
%  radians
%
% OUTPUT:
%   a   [1x1] Semi-major axis           [km]
%   e   [1x1] Eccentricity              [-]
%   i   [1x1] Inclination               [rad]
%   OM  [1x1] RAAN                      [rad]
%   om  [1x1] Pericentre anomaly        [rad]
%   th  [1x1] True anomaly              [rad]
%   mu  [1x1] Gravitational parameter   [km^3/s^2]
%
% INPUT:
%   r   [3x1] Position vector           [km]
%   v   [3x1] Velocity vector           [km/s]
%
% 
%

%Reference system
I = [1 0 0]';
J = [0 1 0]';
K = [0 0 1]';


r_norm = norm(r); % magnitude of position vector[km]
 
v_norm = norm(v);                                                           % magnitude of veocity vector [km/s]

a = 1 / (2 / r_norm - v_norm^2 / mu);                                       % semi-major axis [km]

h_vect = cross(r,v);                                                        % specific angular momentum vector [km^2/s]
h = norm(h_vect);

e_vect = 1 / mu * ( cross(v,h_vect) - mu * r / r_norm );                    % eccentricity vector [-]
e = norm(e_vect);
 
i = acos ( dot(h_vect,K) / h );    % orbit's inclination [rad]


%Check on the inclination
if i==0
  N_vect = I;
  n_vect = I; 
else
N_vect = cross(K,h_vect);
n_vect = cross(K,h_vect)/norm(cross(K,h_vect)); % Node line
end 


%check on the eccentricity
if e < 1e-10
    e=0; 
    e_vect = N_vect; 
    om = 0; 
else
    om = acos( dot(n_vect,e_vect) / e );                                    % argument of pericentre [rad]

    if ( dot(e_vect,K) < 0 )

    om = 2*pi - om;

    end

end

%RAAN 

OM = acos ( dot(n_vect,I) ); 


if ( dot(n_vect,J) < 0 )
    OM = 2*pi - OM;
end


check = abs(dot(r,e_vect) - ( r_norm * e )) ; 
    
if check < 1e-8
    theta = 0; 
else 

    if e ==0

    theta = acos( dot(r,e_vect) / (r_norm) ); 

    else

    theta = acos( dot(r,e_vect) / ( r_norm * e ) ); 

    end 

end 

if (dot(v,r) < -1e-8)
    theta = 2*pi - theta;
end
end
%
%
%
function DisegnaOrbita(a,e,w,i,omega)

p=a*(1-e^2);
theta_dis=linspace(0,2*pi,1000);
R= p./(1+e.*cos(theta_dis));
x=R.*cos(theta_dis);
y=R.*sin(theta_dis);
z=0.*theta_dis;
Rw=[cos(w) sin(w) 0;-sin(w) cos(w) 0; 0 0 1];
Ri=[1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];
Romega=[cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0;0 0 1];
Rge2pf=Rw*Ri*Romega;
Rpf2ge=inv(Rge2pf);
for j=1:size(theta_dis,2)
    rpf=[x(j);y(j);z(j)];
    r_dis=Rpf2ge*rpf;
    x_dis(j)=r_dis(1);
    y_dis(j)=r_dis(2);
    z_dis(j)=r_dis(3);
end
plot3(x_dis,y_dis,z_dis,'linewidth',2);
end
%
%
%
% function h=error_ellipse(varargin)
% % ERROR_ELLIPSE - plot an error ellipse, or ellipsoid, defining confidence region
% %    ERROR_ELLIPSE(C22) - Given a 2x2 covariance matrix, plot the
% %    associated error ellipse, at the origin. It returns a graphics handle
% %    of the ellipse that was drawn.
% %
% %    ERROR_ELLIPSE(C33) - Given a 3x3 covariance matrix, plot the
% %    associated error ellipsoid, at the origin, as well as its projections
% %    onto the three axes. Returns a vector of 4 graphics handles, for the
% %    three ellipses (in the X-Y, Y-Z, and Z-X planes, respectively) and for
% %    the ellipsoid.
% %
% %    ERROR_ELLIPSE(C,MU) - Plot the ellipse, or ellipsoid, centered at MU,
% %    a vector whose length should match that of C (which is 2x2 or 3x3).
% %
% %    ERROR_ELLIPSE(...,'Property1',Value1,'Name2',Value2,...) sets the
% %    values of specified properties, including:
% %      'C' - Alternate method of specifying the covariance matrix
% %      'mu' - Alternate method of specifying the ellipse (-oid) center
% %      'conf' - A value betwen 0 and 1 specifying the confidence interval.
% %        the default is 0.5 which is the 50% error ellipse.
% %      'scale' - Allow the plot the be scaled to difference units.
% %      'style' - A plotting style used to format ellipses.
% %      'clip' - specifies a clipping radius. Portions of the ellipse, -oid,
% %        outside the radius will not be shown.
% %
% %    NOTES: C must be positive definite for this function to work properly.
% default_properties = struct(...
%   'C', [], ... % The covaraince matrix (required)
%   'mu', [], ... % Center of ellipse (optional)
%   'conf', 0.5, ... % Percent confidence/100
%   'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
%   'style', '', ...  % Plot style
%   'clip', inf); % Clipping radius
% if length(varargin) >= 1 & isnumeric(varargin{1})
%   default_properties.C = varargin{1};
%   varargin(1) = [];
% end
% if length(varargin) >= 1 & isnumeric(varargin{1})
%   default_properties.mu = varargin{1};
%   varargin(1) = [];
% end
% if length(varargin) >= 1 & isnumeric(varargin{1})
%   default_properties.conf = varargin{1};
%   varargin(1) = [];
% end
% if length(varargin) >= 1 & isnumeric(varargin{1})
%   default_properties.scale = varargin{1};
%   varargin(1) = [];
% end
% if length(varargin) >= 1 & ~ischar(varargin{1})
%   error('Invalid parameter/value pair arguments.') 
% end
% prop = getopt(default_properties, varargin{:});
% C = prop.C;
% if isempty(prop.mu)
%   mu = zeros(length(C),1);
% else
%   mu = prop.mu;
% end
% conf = prop.conf;
% scale = prop.scale;
% style = prop.style;
% if conf <= 0 | conf >= 1
%   error('conf parameter must be in range 0 to 1, exclusive')
% end
% [r,c] = size(C);
% if r ~= c | (r ~= 2 & r ~= 3)
%   error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
% end
% x0=mu(1);
% y0=mu(2);
% % Compute quantile for the desired percentile
% k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)
% hold_state = get(gca,'nextplot');
% if r==3 & c==3
%   z0=mu(3);
%   
%   % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
%   if any(eig(C) <=0)
%     error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
%   end
%   % C is 3x3; extract the 2x2 matricies, and plot the associated error
%   % ellipses. They are drawn in space, around the ellipsoid; it may be
%   % preferable to draw them on the axes.
%   Cxy = C(1:2,1:2);
%   Cyz = C(2:3,2:3);
%   Czx = C([3 1],[3 1]);
%   [x,y,z] = getpoints(Cxy,prop.clip);
%   h1=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
%   [y,z,x] = getpoints(Cyz,prop.clip);
%   h2=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
%   [z,x,y] = getpoints(Czx,prop.clip);
%   h3=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
%   
%   [eigvec,eigval] = eig(C);
%   [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
%   XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
%   
%   X(:) = scale*(k*XYZ(:,1)+x0);
%   Y(:) = scale*(k*XYZ(:,2)+y0);
%   Z(:) = scale*(k*XYZ(:,3)+z0);
%   h4=surf(X,Y,Z);
%   colormap gray
%   alpha(0.3)
%   camlight
%   if nargout
%     h=[h1 h2 h3 h4];
%   end
% elseif r==2 & c==2
%   % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
%   if any(eig(C) <=0)
%     error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
%   end
%   [x,y,z] = getpoints(C,prop.clip);
%   h1=plot(scale*(x0+k*x),scale*(y0+k*y),prop.style);
%   set(h1,'zdata',z+1)
%   if nargout
%     h=h1;
%   end
% else
%   error('C (covaraince matrix) must be specified as a 2x2 or 3x3 matrix)')
% end
% %axis equal
% set(gca,'nextplot',hold_state);
% end
% %---------------------------------------------------------------
% % getpoints - Generate x and y points that define an ellipse, given a 2x2
% %   covariance matrix, C. z, if requested, is all zeros with same shape as
% %   x and y.
% function [x,y,z] = getpoints(C,clipping_radius)
% n=100; % Number of points around ellipse
% p=0:pi/n:2*pi; % angles around a circle
% [eigvec,eigval] = eig(C); % Compute eigen-stuff
% xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
% x = xy(:,1);
% y = xy(:,2);
% z = zeros(size(x));
% % Clip data to a bounding radius
% if nargin >= 2
%   r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
%   x(r > clipping_radius) = nan;
%   y(r > clipping_radius) = nan;
%   z(r > clipping_radius) = nan;
% end
% end
% %---------------------------------------------------------------
% function x=qchisq(P,n)
% % QCHISQ(P,N) - quantile of the chi-square distribution.
% if nargin<2
%   n=1;
% end
% s0 = P==0;
% s1 = P==1;
% s = P>0 & P<1;
% x = 0.5*ones(size(P));
% x(s0) = -inf;
% x(s1) = inf;
% x(~(s0|s1|s))=nan;
% for ii=1:14
%   dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
%   x(s) = x(s)+dx;
%   if all(abs(dx) < 1e-6)
%     break;
%   end
% end
% end
% %---------------------------------------------------------------
% function F=pchisq(x,n)
% % PCHISQ(X,N) - Probability function of the chi-square distribution.
% if nargin<2
%   n=1;
% end
% F=zeros(size(x));
% if rem(n,2) == 0
%   s = x>0;
%   k = 0;
%   for jj = 0:n/2-1;
%     k = k + (x(s)/2).^jj/factorial(jj);
%   end
%   F(s) = 1-exp(-x(s)/2).*k;
% else
%   for ii=1:numel(x)
%     if x(ii) > 0
%       F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
%     else
%       F(ii) = 0;
%     end
%   end
% end
% end
% %---------------------------------------------------------------
% function f=dchisq(x,n)
% % DCHISQ(X,N) - Density function of the chi-square distribution.
% if nargin<2
%   n=1;
% end
% f=zeros(size(x));
% s = x>=0;
% f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
% end
% %---------------------------------------------------------------
% function properties = getopt(properties,varargin)
% %GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
% %
% %   getopt(properties,varargin) returns a modified properties structure,
% %   given an initial properties structure, and a list of paired arguments.
% %   Each argumnet pair should be of the form property_name,val where
% %   property_name is the name of one of the field in properties, and val is
% %   the value to be assigned to that structure field.
% %
% %   No validation of the values is performed.
% %
% % EXAMPLE:
% %   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
% %   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% % would return:
% %   properties = 
% %         zoom: 1
% %       aspect: 0.7600
% %        gamma: 1
% %         file: 'mydata.dat'
% %           bg: []
% %
% % Typical usage in a function:
% %   properties = getopt(properties,varargin{:})
% % Process the properties (optional input arguments)
% prop_names = fieldnames(properties);
% TargetField = [];
% for ii=1:length(varargin)
%   arg = varargin{ii};
%   if isempty(TargetField)
%     if ~ischar(arg)
%       error('Propery names must be character strings');
%     end
%     f = find(strcmp(prop_names, arg));
%     if length(f) == 0
%       error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
%     end
%     TargetField = arg;
%   else
%     % properties.(TargetField) = arg; % Ver 6.5 and later only
%     properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
%     TargetField = '';
%   end
% end
% if ~isempty(TargetField)
%   error('Property names and values must be specified in pairs.');
% end
% end
