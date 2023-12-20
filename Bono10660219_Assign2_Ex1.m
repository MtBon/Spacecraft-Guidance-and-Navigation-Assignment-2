% Spacecraft Guidance and Navigation (2023/2024)
% Assigment #2
% Exercise #1
% Author: Matteo Bono
%
%% Point 1

clc; clearvars; close all

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

% Retrieve Propagated States
rr_sat1 = xx_sat1(:,1:3);
vv_sat1 = xx_sat1(:,4:6);

rr_sat2 = xx_sat2(:,1:3);
vv_sat2 = xx_sat2(:,4:6);

trace_pos_sat1 = zeros(1,length(et));
trace_vel_sat1 = zeros(1,length(et));

trace_pos_sat2 = zeros(1,length(et));
trace_vel_sat2 = zeros(1,length(et));

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
    trace_pos_sat1(i) = sqrt(trace(P_sat1(1:3,1:3,i)));
    trace_vel_sat1(i) = sqrt(trace(P_sat1(4:6,4:6,i)));

    % Sat 2
    P_sat2(:,:,i) = PHI_sat2(:,:,i) * P0_sat2 * PHI_sat2(:,:,i)';
    P0_sat2 = P_sat2(:,:,i);
    trace_pos_sat2(i) = sqrt(trace(P_sat2(1:3,1:3,i)));
    trace_vel_sat2(i) = sqrt(trace(P_sat2(4:6,4:6,i)));

end

%%

% Unscented Transform
alpha = 0.1;
n = 6;
lambda = alpha^2 * n;
beta = 2;


% Initialize Sigma points vector
  x_sigma_sat1 = zeros(6,13);
  x_sigma_sat2 = zeros(6,13);
  sigma_p1 = zeros(6,13,length(et)-1);
  sigma_p2 = zeros(6,13,length(et)-1);

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
    sigma_p1(:,:,j-1) = x_sigma_sat1(:,:);
    sigma_p2(:,:,j-1) = x_sigma_sat2(:,:);

   
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

num_sim  = 400;


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


    end

        % Sample mean
        mc_mean_sat1(j,:) = mean(xx_mc_sat1,2);
        mc_mean_sat2(j,:) = mean(xx_mc_sat2,2);

          mc_cov_sat1(:,:,j) = cov(xx_mc_sat1');
          mc_cov_sat2(:,:,j) = cov(xx_mc_sat2');

          % Generate points distribution at this time step
          mc_sat1(:,:,j) = xx_mc_sat1';
          mc_sat2(:,:,j) = xx_mc_sat2';


end

fprintf('Finished\n');



for i = 1 : length(et)
    % Values to plot
    cov_pos_sat1_mc(i) = 3 * sqrt(max(eig(mc_cov_sat1(1:3,1:3,i))));
    cov_vel_sat1_mc(i) = 3 * sqrt(max(eig(mc_cov_sat1(4:6,4:6,i))));

    cov_pos_sat2_mc(i) = 3 * sqrt(max(eig(mc_cov_sat2(1:3,1:3,i))));
    cov_vel_sat2_mc(i) = 3 * sqrt(max(eig(mc_cov_sat2(4:6,4:6,i))));

end
%%
% Compute Row vectors
r = x0_mean_sat1(1:3);
v = x0_mean_sat1(4:6);

i_ax = r/norm(r);
k_ax = cross(r,v)/norm(cross(r,v));
j_ax = cross(k_ax,i_ax);

% Rotation Matrix
rotm = [i_ax'; j_ax'; k_ax'];

% Initialize
r_mc1 = zeros(num_sim,3,length(et));
r_mc2 = zeros(num_sim,3,length(et));
mc_cov_sat1_lvlh = zeros(3,3,length(et));
mc_cov_sat2_lvlh = zeros(3,3,length(et));
mc_mean_lvlh_sat1 = zeros(3,length(et));
mc_mean_lvlh_sat2 = zeros(3,length(et));
r_lin1 = zeros(3,length(et));
r_lin2 = zeros(3,length(et));
P_lvlh_sat1 = zeros(3,3,length(et));
P_lvlh_sat2 = zeros(3,3,length(et));
sigma_p1_lvlh = zeros(3,13,length(et)-1);
sigma_p2_lvlh = zeros(3,13,length(et)-1);
sigma_mean1_lvlh = zeros(3,length(et)-1);
sigma_mean2_lvlh = zeros(3,length(et)-1);
sigma_cov1_lvlh = zeros(3,3,length(et)-1);
sigma_cov2_lvlh = zeros(3,3,length(et)-1);

% Projection for MC points
for j = 1 : length(et)

    for i = 1 : num_sim
       
        r_mc1(i,:,j) = rotm * mc_sat1(i,1:3,j)';
        r_mc2(i,:,j) = rotm * mc_sat2(i,1:3,j)';
        
    end
    mc_cov_sat1_lvlh(:,:,j) = rotm * mc_cov_sat1(1:3,1:3,j) * rotm';
    mc_cov_sat2_lvlh(:,:,j) = rotm * mc_cov_sat2(1:3,1:3,j) * rotm';

    mc_mean_lvlh_sat1(:,j) = rotm * mc_mean_sat1(j,1:3)';
    mc_mean_lvlh_sat2(:,j) = rotm * mc_mean_sat2(j,1:3)';

    % Projection for LinCov
    r_lin1(:,j) = rotm * rr_sat1(j,:)';
    r_lin2(:,j) = rotm * rr_sat2(j,:)';
    P_lvlh_sat1(:,:,j) = rotm * P_sat1(1:3,1:3,j) * rotm';
    P_lvlh_sat2(:,:,j) = rotm * P_sat2(1:3,1:3,j) * rotm';


end

% Projection for UT
for j = 1 : length(et)-1

    for i = 1:13

        sigma_p1_lvlh(:,i,j) = rotm * sigma_p1(1:3,i,j);
        sigma_p2_lvlh(:,i,j) = rotm * sigma_p2(1:3,i,j);
    end

    sigma_mean1_lvlh(:,j) = rotm * sigma_mean_sat1(1:3,j);
    sigma_mean2_lvlh(:,j) = rotm * sigma_mean_sat2(1:3,j);

    sigma_cov1_lvlh(:,:,j) = rotm * sigma_cov_sat1(1:3,1:3,j) * rotm';
    sigma_cov2_lvlh(:,:,j) = rotm * sigma_cov_sat2(1:3,1:3,j) * rotm';
end



%% Plots

tp = linspace(1,11,11);

% Trace pos and Vel for Point 1
figure()
subplot(1,2,1)
plot(et,trace_pos_sat1,et,trace_pos_sat2,'LineWidth',2);
legend('Mango','Tango');
title('Position')

subplot(1,2,2)
plot(et,trace_vel_sat1,et,trace_vel_sat2,'LineWidth',2);
legend('Mango','Tango');
title('Velocity')

figure()
plot(tp,Dr_UT,tp,Dr_Lin,'LineWidth',2);
xlabel('Orbital Periods');
ylabel('[Km]')
legend('$\Delta r_{Lin}$','$\Delta r_{UT}$','interpreter','latex')

% Plots for Point 3
figure()
subplot(1,2,1)
plot(tp,cov_pos_sat1_Lin,tp,cov_pos_sat1_UT,tp,cov_pos_sat1_mc,'LineWidth',2)
xlabel('Orbital Period');
ylabel('[-]');
title('3$\sigma$ Pos for Mango','Interpreter','Latex');
legend('LinCov','UT','MC');

subplot(1,2,2)
plot(tp,cov_pos_sat2_Lin,tp,cov_pos_sat2_UT,tp,cov_pos_sat2_mc,'LineWidth',2)
xlabel('Orbital Period');
ylabel('[-]');
title('3$\sigma$ Pos for Tango','Interpreter','Latex');
legend('LinCov','UT','MC');

figure()
subplot(1,2,1)
plot(tp,cov_vel_sat1_Lin,tp,cov_vel_sat1_UT,tp,cov_vel_sat1_mc,'LineWidth',2)
xlabel('Orbital Period');
ylabel('[-]');
title('3$\sigma$ Vel for Mango','Interpreter','Latex');
legend('LinCov','UT','MC');


subplot(1,2,2)
plot(tp,cov_vel_sat2_Lin,tp,cov_vel_sat2_UT,tp,cov_vel_sat2_mc,'LineWidth',2)
xlabel('Orbital Period');
ylabel('[-]');
title('3$\sigma$ Vel for Tango','Interpreter','Latex');
legend('LinCov','UT','MC');


% MC Points in Orbital Plane
figure()
subplot(1,2,1)
hold on;
scatter(r_mc1(:,1,1),r_mc1(:,2,1),14);
scatter(r_mc1(:,1,4),r_mc1(:,2,4),14);
scatter(r_mc1(:,1,end),r_mc1(:,2,end),14);

drawEllipse(mc_mean_lvlh_sat1(1:2,1), mc_cov_sat1_lvlh(1:2,1:2,1),3)
drawEllipse(mc_mean_lvlh_sat1(1:2,4), mc_cov_sat1_lvlh(1:2,1:2,4),3)
drawEllipse(mc_mean_lvlh_sat1(1:2,end), mc_cov_sat1_lvlh(1:2,1:2,end),3)
legend('','','','Mango:Covariance @Rev:1','Mango:Covariance @Rev:4','Mango:Covariance at tf');
xlabel('X[Km]');
ylabel('Y[Km]');
title('MC: Mango');

subplot(1,2,2)
hold on;
scatter(r_mc2(:,1,1),r_mc2(:,2,1),14);
scatter(r_mc2(:,1,4),r_mc2(:,2,4),14);
scatter(r_mc2(:,1,end),r_mc2(:,2,end),14);

drawEllipse(mc_mean_lvlh_sat2(1:2,1), mc_cov_sat2_lvlh(1:2,1:2,1),3)
drawEllipse(mc_mean_lvlh_sat2(1:2,4), mc_cov_sat2_lvlh(1:2,1:2,4),3)
drawEllipse(mc_mean_lvlh_sat2(1:2,end), mc_cov_sat2_lvlh(1:2,1:2,end),3)
legend('','','','Tango:Covariance @Rev:1','Tango:Covariance @Rev:4','Tango:Covariance at tf');
title('MC: Tango');
xlabel('X[Km]');
ylabel('Y[Km]');

% UT Points
figure()
subplot(1,2,1)
hold on;
scatter(sigma_p1_lvlh(1,:,1),sigma_p1_lvlh(2,:,1),20)
scatter(sigma_p1_lvlh(1,:,4),sigma_p1_lvlh(2,:,4),20)
scatter(sigma_p1_lvlh(1,:,end),sigma_p1_lvlh(2,:,end),20)

drawEllipse(sigma_mean1_lvlh(1:2,1), sigma_cov1_lvlh(1:2,1:2,1),3)
drawEllipse(sigma_mean1_lvlh(1:2,4), sigma_cov1_lvlh(1:2,1:2,4),3)
drawEllipse(sigma_mean1_lvlh(1:2,end), sigma_cov1_lvlh(1:2,1:2,end),3)
legend('','','','Covariance @Rev:1','Covariance @Rev:4','Covariance at tf');
title('UT: Mango');
xlabel('X[Km]');
ylabel('Y[Km]');

subplot(1,2,2)
hold on;
scatter(sigma_p2_lvlh(1,:,1),sigma_p2_lvlh(2,:,1),20)
scatter(sigma_p2_lvlh(1,:,4),sigma_p2_lvlh(2,:,4),20)
scatter(sigma_p2_lvlh(1,:,end),sigma_p2_lvlh(2,:,end),20)

drawEllipse(sigma_mean2_lvlh(1:2,1), sigma_cov2_lvlh(1:2,1:2,1),3)
drawEllipse(sigma_mean2_lvlh(1:2,4), sigma_cov2_lvlh(1:2,1:2,4),3)
drawEllipse(sigma_mean2_lvlh(1:2,end), sigma_cov2_lvlh(1:2,1:2,end),3)
legend('','','','Covariance @Rev:1','Covariance @Rev:4','Covariance at tf');
title('UT: Tango');
xlabel('X[Km]');
ylabel('Y[Km]');

%LinCov
figure()
subplot(1,2,1)
hold on;
scatter(r_lin1(1,1),r_lin1(2,1),30);
scatter(r_lin1(1,4),r_lin1(2,4),30)
scatter(r_lin1(1,end),r_lin1(2,end),30)

drawEllipse(r_lin1(1:2,1), P_lvlh_sat1(1:2,1:2,1),3)
drawEllipse(r_lin1(1:2,4), P_lvlh_sat1(1:2,1:2,4),3)
drawEllipse(r_lin1(1:2,end), P_lvlh_sat1(1:2,1:2,end),3)
legend('','','','Covariance @Rev:1','Covariance @Rev:4','Covariance at tf');
title('LinCov: Mango');
xlabel('X[Km]');
ylabel('Y[Km]');

subplot(1,2,2)
hold on;
scatter(r_lin2(1,1),r_lin2(2,1),30);
scatter(r_lin2(1,4),r_lin2(2,4),30)
scatter(r_lin2(1,end),r_lin2(2,end),30)

drawEllipse(r_lin2(1:2,1), P_lvlh_sat2(1:2,1:2,1),3)
drawEllipse(r_lin2(1:2,4), P_lvlh_sat2(1:2,1:2,4),3)
drawEllipse(r_lin2(1:2,end), P_lvlh_sat2(1:2,1:2,end),3)

legend('','','','Covariance @Rev:1','Covariance @Rev:4','Covariance at tf');
title('LinCov: Tango');
xlabel('X[Km]');
ylabel('Y[Km]');
%% Functions

function [xf, PHI_f, tt, xx ] = keplerian_propagator_STM( et0,x0, et1 , attractor)

%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 

GM = cspice_bodvrd(attractor, 'GM', 1);

options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12*ones(42,1));


x0Phi0 = [ x0 ; reshape(eye(6),36,1) ];

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_STM_rhs(t,x,GM),[et0 et1], x0Phi0, options_STM);

% Extract state vector 
xf = xx(end,1:6);
PHI_f = reshape(xx(end,7:end),6,6);

end

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

function [xf, tt, xx] = keplerian_propagator(t0, x0, t1 , attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 


    GM = cspice_bodvrd(attractor, 'GM', 1);


options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [t0 t1], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

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

function [varargout] = drawEllipse(mean,P,n_sigma)
% drawEllipse - Draws an ellipse based on mean and covariance matrix.
%
%   Inputs:
%       - mean: Mean vector of the distribution.
%               Column vector.
%       - P: Covariance matrix.
%            Square matrix.
%       - n_sigma: Number of standard deviations for ellipse size.
%                  Positive scalar.
%
%   Outputs:
%       - If called with no output arguments, the function plots the ellipse.
%       - If called with two output arguments, the function returns the x and y
%         coordinates of the ellipse without plotting.
%

% Build circle
n_points = 200;
alpha = 2 * pi / n_points * (0:n_points);
circle = [cos(alpha); sin(alpha)];

% Singular Value Decomposition (SVD) to find ellipse parameters
[R, D] = svd(P);
d = sqrt(D);

% Transform the circle to an aligned ellipse, then rotate it, and finally
% scale it based on the specified number of standard deviations
ellipse = n_sigma * R * d * circle;

% Shift the ellipse to the mean position
x = mean(1) + ellipse(1, :);
y = mean(2) + ellipse(2, :);

if nargout == 2
    varargout{1} = x;
    varargout{2} = y;
else
    plot(x, y);
end

end
