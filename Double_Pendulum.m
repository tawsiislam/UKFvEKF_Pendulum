clear;
close all;
clf;

% Parameters that can be tweaked to test the simulation

R = 0.0001*eye(4); %Process noise %Low 0.0001 High 0.01 Ratio1 R:High Q: Low, Ratio2 R:Low Q:High
Q = 0.0001*eye(2); %Observation noise (we only observe angle) Low 0.0001 High 0.01
theta1_0 = 3*pi/8;    %Low pi/4 High pi/2
theta2_0 = 3*pi/8;    %Low pi/4 High pi/2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterations = 100;
dt = 0.1;
g = 9.82;
L1 = 5;
L2 = 5;
L = [L1 L2];
m1 = 1;
m2 = 1;
m = [m1 m2];

d = 0;

sus_params = {g,L,R,Q,d,m};

omega1_0 = 0;
omega2_0 = 0;


theta1 = zeros(iterations,1);
theta2 = zeros(iterations,1);
omega1 = zeros(iterations,1);
omega2 = zeros(iterations,1);
theta1(1) = theta1_0;
theta2(1) = theta2_0;
omega1(1) = omega1_0;
omega2(1) = omega2_0;

% State variables initialized
state_est = zeros(4,iterations);
sigmas = zeros(4,4,iterations);
%state_est_hope = zeros(4,iterations);
sigmas_hope = zeros(4,4,iterations);

sigmas(:,:,1) = 4*eye(4);
sigma_prev = reshape(sigmas(:,:,1),4,4);
sigmas_hope(:,:,1) = 4*eye(4);
sigma_prev_hope = reshape(sigmas_hope(:,:,1),4,4);
sigmas_ekf = zeros(4,4,iterations);
sigmas_ekf(:,:,1) = 4*eye(4);
Sigma_prev_ekf = sigma_prev;

state_est(1,1) = theta1_0; % + 0*normrnd(0,0.1);
state_est(2,1) = theta2_0; %+ 0*normrnd(0,0.1);
state_est(3,1) = omega1_0; %+ 0*normrnd(0,0.1);
state_est(4,1) = omega2_0; %+ 0*normrnd(0,0.1);
state_est_hope = state_est;

energy_vec = zeros(iterations,1);
PotEng = zeros(iterations,1);
KinEng = zeros(iterations,1);
PotEng(1) = -(m1+m2)*g*L1*cos(theta1(1))-m2*g*L2*cos(theta2(1));
KinEng(1) = 0.5*m1*omega1(1)^2*L1^2+0.5*m2*(omega1(1)^2*L1^2+omega2(1)^2*L2^2 ...
    +2*omega1(1)*omega2(1)*L1*L2*cos(theta1(1)-theta2(1)));
energy_vec(1) = PotEng(1) + KinEng(1);


theta1_ekf = zeros(iterations,1);
theta2_ekf = zeros(iterations,1);
omega1_ekf = zeros(iterations,1);
omega2_ekf = zeros(iterations,1);
theta1_ekf(1) = theta1(1);
theta2_ekf(1) = theta2(1);
omega1_ekf(1) = omega1(1);
omega2_ekf(1) = omega2(1);
error1_ukf = zeros(iterations,1);
error2_ukf = zeros(iterations,1);
error1_ukfhope = zeros(iterations,1);
error2_ukfhope = zeros(iterations,1);
error1_ekf = zeros(iterations,1);
error2_ekf = zeros(iterations,1);
%error_ekf(1) = abs(theta(1)-theta_ekf(1));

%Function to calculate angles and linear momentum for double pendulum
accel1 = @(theta1,theta2,omega1,omega2) (-g*(2*m1+m2)*sin(theta1) - m2*g*sin(theta1-2*theta2)...
        -2*sin(theta1-theta2)*m2*(omega2^2*L2 + omega1^2*L1*cos(theta1 - theta2)))...
        / (L1*(2*m1+m2- m2*cos(2*theta1-2*theta2)));
accel2 = @(theta1,theta2,omega1,omega2) (2*sin(theta1 - theta2)*(omega1^2*L1*(m1+m2)+g*(m1+m2)*cos(theta1)...
        +omega2^2*L2*m2*cos(theta1 - theta2)))/ (L2*(2*m1+m2-m2*cos(2*theta1-2*theta2)));

for n=1:iterations %Runge-Kutta 4 integrator
    curr_theta1 = theta1(n);
    curr_omega1 = omega1(n);
    curr_theta2 = theta2(n);
    curr_omega2 = omega2(n);

    dtheta1_1 = curr_omega1*dt;    %Change in angle using angular velocity
    dtheta2_1 = curr_omega2*dt;
    domega1_1 = accel1(curr_theta1,curr_theta2,curr_omega1,curr_omega2)*dt;    %Change in ang vel using ang acc
    domega2_1 = accel2(curr_theta1,curr_theta2,curr_omega1,curr_omega2)*dt;

    dtheta1_2 = (curr_omega1+domega1_1*0.5)*dt;
    dtheta2_2 = (curr_omega2+domega2_1*0.5)*dt;
    domega1_2 = accel1(curr_theta1+dtheta1_1*0.5,curr_theta2+dtheta2_1*0.5,curr_omega1+domega1_1*0.5,curr_omega2+domega2_1*0.5)*dt;
    domega2_2 = accel2(curr_theta1+dtheta1_1*0.5,curr_theta2+dtheta2_1*0.5,curr_omega1+domega1_1*0.5,curr_omega2+domega2_1*0.5)*dt;

    dtheta1_3 = (curr_omega1+domega1_2*0.5)*dt;
    dtheta2_3 = (curr_omega2+domega2_2*0.5)*dt;
    domega1_3 = accel1(curr_theta1+dtheta1_2*0.5,curr_theta2+dtheta2_2*0.5,curr_omega1+domega1_2*0.5,curr_omega2+domega2_2*0.5)*dt;
    domega2_3 = accel2(curr_theta1+dtheta1_2*0.5,curr_theta2+dtheta2_2*0.5,curr_omega1+domega1_2*0.5,curr_omega2+domega2_2*0.5)*dt;

    dtheta1_4 = (curr_omega1+domega1_3)*dt;
    dtheta2_4 = (curr_omega2+domega2_3)*dt;
    domega1_4 = accel1(curr_theta1+dtheta1_3,curr_theta2+dtheta2_3,curr_omega1+domega1_3,curr_omega2+domega2_3)*dt;
    domega2_4 = accel2(curr_theta1+dtheta1_3,curr_theta2+dtheta2_3,curr_omega1+domega1_3,curr_omega2+domega2_3)*dt;

    theta1(n+1) = theta1(n) + (dtheta1_1+2*dtheta1_2+2*dtheta1_3+dtheta1_4)/6;
    theta2(n+1) = theta2(n) + (dtheta2_1+2*dtheta2_2+2*dtheta2_3+dtheta2_4)/6;
    omega1(n+1) = omega1(n) + (domega1_1+2*domega1_2+2*domega1_3+domega1_4)/6;
    omega2(n+1) = omega2(n) + (domega2_1+2*domega2_2+2*domega2_3+domega2_4)/6;

    PotEng(n+1) = -(m1+m2)*g*L1*cos(theta1(n+1))-m2*g*L2*cos(theta2(n+1));
    KinEng(n+1) = 0.5*m1*omega1(n+1)^2*L1^2+0.5*m2*(omega1(n+1)^2*L1^2+omega2(n+1)^2*L2^2 ...
        +2*omega1(n+1)*omega2(n+1)*L1*L2*cos(theta1(n+1)-theta2(n+1)));
    energy_vec(n+1) = PotEng(n+1) + KinEng(n+1);

    %UKF estimation of the system
    %z = [theta1(n+1); theta2(n+1); omega1(n+1); omega2(n+1)];
    z = [theta1(n+1); theta2(n+1)];
    mu_prev = state_est(:,n);
    [mu, sigma] = UKF(mu_prev, sigma_prev, z, R, Q, m, L, dt);
    state_est(:,n+1) = mu;
    sigmas(:,:,n+1) = sigma;
    sigma_prev = sigma;

    %EKF estimation of system
    mu_prev_ekf = [theta1_ekf(n); theta2_ekf(n); omega1_ekf(n);omega2_ekf(n)];
    [mu_ekf, Sigma_ekf] = EKF(mu_prev_ekf,Sigma_prev_ekf,z,R,Q,dt,m,L);
    theta1_ekf(n+1) = mu_ekf(1);
    theta2_ekf(n+1) = mu_ekf(2);
    omega1_ekf(n+1) = mu_ekf(3);
    omega2_ekf(n+1) = mu_ekf(4);
    sigmas_ekf(:,:,n+1) = Sigma_ekf;
    Sigma_prev_ekf = Sigma_ekf;
    
    ukf_diff1 = abs(theta1(1:n+1) - state_est(1,1:n+1)');
    ukf_diff2 = abs(theta2(1:n+1) - state_est(2,1:n+1)');
    ukfhope_diff1 = abs(theta1(1:n+1) - state_est(1,1:n+1)');
    ukfhope_diff2 = abs(theta2(1:n+1) - state_est(2,1:n+1)');
    ekf_diff1 = abs(theta1(1:n+1) - theta1_ekf(1:n+1));
    ekf_diff2 = abs(theta2(1:n+1) - theta2_ekf(1:n+1));
    error1_ukf(n) = mean(ukf_diff1);
    error2_ukf(n) = mean(ukf_diff2);
    error1_ukfhope(n) = mean(ukfhope_diff1);
    error2_ukfhope(n) = mean(ukfhope_diff2);
    error1_ekf(n) = mean(ekf_diff1);
    error2_ekf(n) = mean(ekf_diff2);


end


t = linspace(0,dt*iterations,iterations);
colorSim = [0.4667    0.6745    0.1882];
colorMass = [0.9373    0.6627    0.3373];
colorLine = [0.3176    0.4784    0.5922];
colorUKF = [0    0.4471    0.7412];
colorEKF = [0.8510    0.3255    0.0980];

% f = figure(1);
%f.Position = [200 100 1200 500];
%h1 = figure(1);
h2 = figure(1);
h3 = figure(2);
h4 = figure(3);


x1_vec = L1*sin(theta1);
x2_vec = L1*sin(theta1)+L2*sin(theta2);
y1_vec = -L1*cos(theta1);
y2_vec = -L1*cos(theta1)-L2*cos(theta2);
x2_vestim = L1*sin(state_est(1,:))+L2*sin(state_est(2,:));
y2_vestim = -L1*cos(state_est(1,:))-L2*cos(state_est(2,:));
x2_vekf = L1*sin(theta1_ekf)+L2*sin(theta2_ekf);
y2_vekf = -L1*cos(theta1_ekf)-L2*cos(theta2_ekf);
error = sqrt((x2_vec-x2_vestim).^2+(y2_vec-y2_vestim).^2);

movieVector = struct('cdata', cell(1, length(t)-1), 'colormap', cell(1, length(t)-1));

min_t = min(t);
max_t = max(t);
max_x = (L1+L2)+1;
min_y = -max_x;
max_y = 0;
max_err1 = max([max(abs(error1_ukf)) max(abs(error1_ekf))])*1.1;
max_err2 = max([max(abs(error2_ukf)) max(abs(error2_ekf))])*1.1;
max_th1 = max(abs(theta1))*1.1;
max_th2 = max(abs(theta2))*1.1;
max_cov1 = max([max(sqrt(reshape(sigmas(1,1,2:end), 1,iterations))) max(sqrt(reshape(sigmas_ekf(1,1,2:end), 1,iterations)))])*1.1;
max_cov2 = max([max(sqrt(reshape(sigmas(2,2,2:end), 1,iterations))) max(sqrt(reshape(sigmas_ekf(2,2,2:end), 1,iterations)))])*1.1;

for k=length(t):length(t)
    figure(1)
    clf;
    t_k = t(k);
    x1_pos = x1_vec(k);
    y1_pos = y1_vec(k);
    x2_pos = x2_vec(k);
    y2_pos = y2_vec(k);
    x_estim = x2_vestim(k);
    y_estim = y2_vestim(k);
    x_ekf = x2_vekf(k);
    y_ekf = y2_vekf(k);

    % figure(h1)
    %subplot(4,2,[1 3 5 7])
%     set(get(gca,'XAxis'),'FontSize', 10)
%     set(get(gca,'YAxis'),'FontSize', 10)
%     plot([0 x1_pos],[0 y1_pos],'color',colorLine,"Linewidth",3)
%     hold on
%     plot([x1_pos x2_pos],[y1_pos y2_pos],'color',colorLine,"Linewidth",3)
%     hold on
%     scatter(x1_pos,y1_pos,70,colorMass,'filled'); %Mass
%     hold on
%     scatter(x2_pos,y2_pos,70,colorMass,'filled'); %Mass
%     hold on
%     scatter(x2_vestim(k),y2_vestim(k),60,colorMarker,'x','LineWidth',2); %Estimated position (UKF)
%     hold on
%     scatter(x2_vekf(k),y2_vekf(k),60,colorEKF,'x','LineWidth',2); %Estimated position (EKF)
%     xlabel('x')
%     ylabel('y')
%     xlim([-max_x max_x])
%     ylim([min_y -min_y])
%     xline(0,'--');
%     yline(0,'--');
%     title(sprintf('Double Pendulum at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',12)
    
    figure(h2)
    subplot(2,1,1)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),theta1(1:k),'--','color',colorSim,"Linewidth",3)
    hold on 
    plot(t(1:k),state_est(1,1:k),'color',colorUKF,"Linewidth",2)
    hold on 
    plot(t(1:k),theta1_ekf(1:k),'color',colorEKF,"Linewidth",2)
    xlabel('time t','FontSize',12)
    ylabel(sprintf('${\\theta_1}$'),'Interpreter','latex','FontSize',18)
    xlim([min_t max_t])
    ylim([-max_th1 max_th1])
    title(sprintf('${\\theta_1}$ at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('Simulation','UKF estim','EKF estim')

    subplot(2,1,2)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),theta2(1:k),'--','color',colorSim,"Linewidth",3)
    hold on 
    plot(t(1:k),state_est(2,1:k),'color',colorUKF,"Linewidth",2)
    hold on 
    plot(t(1:k),theta2_ekf(1:k),'color',colorEKF,"Linewidth",2)
    xlabel('time t','FontSize',12)
    ylabel(sprintf('${\\theta_2}$'),'Interpreter','latex','FontSize',18)
    xlim([min_t max_t])
    ylim([-max_th2 max_th2])
    title(sprintf('${\\theta_2}$ at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('Simulation','UKF estim','EKF estim')
    
    %Figure for error
    figure(h3)
    subplot(2,1,1)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),abs(error1_ukf(1:k)),'color',colorUKF,"Linewidth",3)
    hold on 
    plot(t(1:k),abs(error1_ekf(1:k)),'color',colorEKF,"Linewidth",3)
    xlabel('time t','FontSize',12)
    ylabel('Error','FontSize',12)
    xlim([min_t max_t])
    ylim([0 max_err1])
    %yline(0,'--');
    title(sprintf('Error for ${\\theta_1}$ after t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('Error UKF \theta_1', 'Error EKF \theta_1','Location','southeast');

    subplot(2,1,2)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),abs(error2_ukf(1:k)),'color',colorUKF,"Linewidth",3)
    hold on
    plot(t(1:k),abs(error2_ekf(1:k)),'color',colorEKF,"Linewidth",3)
    xlabel('time t','FontSize',12)
    ylabel('Error','FontSize',12)
    xlim([min_t max_t])
    ylim([0 max_err2])
    %yline(0,'--');
    title(sprintf('Error for ${\\theta_2}$ after t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('Error UKF \theta_2', 'Error EKF \theta_2','Location','southeast');

    figure(h4)
    subplot(2,1,1)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),sqrt(reshape(sigmas(1,1,1:k), 1,k)),'color',colorUKF,"Linewidth",3)
    hold on
    plot(t(1:k),sqrt(reshape(sigmas_ekf(1,1,1:k),1,k)),'color',colorEKF,"Linewidth",3)
    xlabel('time t','FontSize',12)
    ylabel(sprintf('${\\sigma_1}$'),'Interpreter','latex','FontSize',18)
    xlim([min_t max_t])
    ylim([0 max_cov1])
    %yline(0,'--');
    title(sprintf('Standard deviation ${\\sigma_1}$ after t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('St. dev. UKF \theta_1', 'St. dev EKF \theta_1','Location','northeast')

    subplot(2,1,2)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),sqrt(reshape(sigmas(2,2,1:k), 1,k)),'color',colorUKF,"Linewidth",3)
    hold on
    plot(t(1:k),sqrt(reshape(sigmas_ekf(2,2,1:k),1,k)),'color',colorEKF,"Linewidth",3)
    xlabel('time t','FontSize',12)
    ylabel(sprintf('${\\sigma_2}$'),'Interpreter','latex','FontSize',18)
    xlim([min_t max_t])
    ylim([0 max_cov2])
    %yline(0,'--');
    title(sprintf('Standard deviation ${\\sigma_2}$ after t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('St. dev. UKF \theta_2', 'St. dev EKF \theta_2','Location','northeast')
    

    
    %movieVector(k) = getframe(gcf);
    drawnow % Used to run simulation. Comment to only see the end result.
end
%% Comment out this block to make a video
% myWriter = VideoWriter('Double_Pendulum','MPEG-4');
% myWriter.FrameRate= 20;
% open(myWriter);
% writeVideo(myWriter, movieVector);
% close(myWriter);