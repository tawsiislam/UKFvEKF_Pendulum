% Parameters that can be tweaked to test the simulation

R = 0.0001*eye(2); %Process noise %Low 0.0001 High 0.01 Ratio1 R:High Q: Low, Ratio2 R:Low Q:High
Q = 0.0001; %Observation noise (we only observe angle) Low 0.0001 High 0.01
theta0 = 3*pi/8;  %Small pi/4 High pi/2 - Initial angle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterations = 100;

omega0 = 0; % Free fall
alpha0 = 0;
dt = 0.1;
g = 9.82;
L = 10; % Length of the pendulum
m = 1; % Mass at the end of the pendulum

% Initializing state variables
theta = zeros(iterations,1);
omega = zeros(iterations,1);
alpha = zeros(iterations,1); %Vector for angular acceleration
theta(1) = theta0;
omega(1) = omega0;
alpha(1) = alpha0;

PotEng = zeros(iterations,1);
KinEng = zeros(iterations,1);
energy_vec = zeros(iterations,1);
PotEng(1) = -m*g*L*cos(theta(1));
KinEng(1) = 0.5*m*omega(1)^2*L^2;
energy_vec(1) = PotEng(1) + KinEng(1);


state_est = zeros(2,iterations);
sigmas = zeros(2,2,iterations);

sigmas(:,:,1) = 4*eye(2);
sigma_prev = reshape(sigmas(:,:,1),2,2);
sigmas_ekf = zeros(2,2,iterations);
sigmas_ekf(:,:,1) = 4*eye(2);
Sigma_prev_ekf = sigma_prev;

state_est(1,1) = theta0 + 0*normrnd(0,0.1);
state_est(2,1) = omega0 + 0*normrnd(0,0.1);

theta_ekf = zeros(iterations,1);
omega_ekf = zeros(iterations,1);
theta_ekf(1) = theta(1);
omega_ekf(1) = omega(1);
error1_ukf = zeros(iterations,1);
error1_ekf = zeros(iterations,1);
%error_ekf(1) = abs(theta(1)-theta_ekf(1));


% Function to calculate angular acceleration with drag
accel = @(theta,omega) -g*sin(theta)/L-0*omega;

for n=1:iterations %Runge-Kutta integrator
    curr_theta = theta(n);
    curr_omega = omega(n);
    a1 = accel(curr_theta,curr_omega)*dt;
    b1 = curr_omega*dt;
    a2 = accel(curr_theta+b1/2,curr_omega+a1/2)*dt;
    b2 = (curr_omega+a1/2)*dt;
    a3 = accel(curr_theta+b2/2,curr_omega+a2/2)*dt;
    b3 = (curr_omega+a2/2)*dt;
    a4 = accel(curr_theta+b3,curr_omega+a3)*dt;
    b4 = (curr_omega+a3)*dt;

    theta(n+1) = theta(n) + (b1+2*b2+2*b3+b4)/6;
    omega(n+1) = omega(n) + (a1+2*a2+2*a3+a4)/6;

    PotEng(n+1) = -m*g*L*cos(theta(n+1));
    KinEng(n+1) = 0.5*m*omega(n+1)^2*L^2;
    energy_vec(n+1) = PotEng(n+1) + KinEng(n+1);

    %UKF estimation of the system
    z = theta(n+1);
    mu_prev = [theta(n);omega(n)];
    [mu, sigma] = UKF(mu_prev, sigma_prev, z, R, Q,m,L,dt);
    state_est(:,n+1) = mu;
    sigmas(:,:,n+1) = sigma;
    sigma_prev = reshape(sigmas(:,:,n),2,2);

    %EKF estimation of system
    mu_prev_ekf = [theta_ekf(n);omega_ekf(n)];
    [mu_ekf, Sigma_ekf] = EKF(mu_prev_ekf,Sigma_prev_ekf,z,R,Q,dt,m,L);
    theta_ekf(n+1) = mu_ekf(1);
    omega_ekf(n+1) = mu_ekf(2);
    sigmas_ekf(:,:,n+1) = Sigma_ekf;
    Sigma_prev_ekf = Sigma_ekf;

    ukf_diff1 = abs(theta(1:n+1) - state_est(1,1:n+1)');
    ekf_diff1 = abs(theta(1:n+1) - theta_ekf(1:n+1));
    error1_ukf(n) = mean(ukf_diff1);
    error1_ekf(n) = mean(ekf_diff1);

    
end

t = linspace(0,dt*iterations,iterations);
colorSim = [0.4667    0.6745    0.1882];
colorMass = [0.9373    0.6627    0.3373];
colorLine = [0.3176    0.4784    0.5922];
colorUKF = [0    0.4471    0.7412];
colorEKF = [0.8510    0.3255    0.0980];

% f = figure(1);
% f.Position = [200 100 1200 500];    %Sets figure's start position and size
h2 = figure(1);
% h3 = figure(2);
% h4 = figure(3);
x_vec = L*sin(theta);
y_vec = -L*cos(theta);

x_vestim = L*sin(state_est(1,:));
y_vestim = -L*cos(state_est(1,:));
%x_vestim = x_vec+normrnd(0,0.1,[length(x_vec) 1]);  %Simulate an estimated position. This vector will be replaced by iEKF and UKF
%y_vestim = y_vec+normrnd(0,0.1,[length(y_vec) 1]);
error = ((theta-state_est(1,:)'))/length(state_est(1,:));
error_ekf = (theta-theta_ekf)/length(theta_ekf);


%movieVector = struct('cdata', cell(1, length(x_vec)-1), 'colormap',
%cell(1, length(x_vec)-1)); % Needed to create a movie

% Sets the plots sizes
min_t = min(t);
max_t = max(t);
max_x = L;
min_y = -2*L;
max_y = 0;
max_th = max(abs(theta_ekf))*1.1;
max_err = max([max(error1_ukf) max(error1_ekf)])*1.1;
max_cov = max([max(sqrt(reshape(sigmas(1,1,2:end), 1,iterations))) max(sqrt(reshape(sigmas_ekf(1,1,2:end), 1,iterations)))])*1.1;

for k=length(t):length(t)
    clf;
    figure(1)
    t_k = t(k);
    x_pos = x_vec(k);
    y_pos = y_vec(k);
    x_estim = x_vestim(k);
    y_estim = y_vestim(k);

    
%     subplot(4,2,[1 3 5 7])
%     plot([0 x_pos],[0 y_pos],'color',colorLine,"Linewidth",3)
%     hold on
%     scatter(x_pos,y_pos,70,colorMass,'filled'); %Mass
%     hold on
%     scatter(x_estim,y_estim,60,colorMarker,'x','LineWidth',2); %Estimated position
%     xlabel('x')
%     ylabel('y')
%     xlim([-max_x max_x])
%     ylim([min_y 0])
%     title(sprintf('Mass at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',12)

    subplot(3,1,1)
    figure(h2)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),theta(1:k),'--','color',colorSim,"Linewidth",3)
    hold on 
    plot(t(1:k),state_est(1,1:k),'color',colorUKF,"Linewidth",2)
    hold on 
    plot(t(1:k),theta_ekf(1:k),'color',colorEKF,"Linewidth",2)
    xlabel('time t', 'FontSize', 12)
    ylabel(sprintf('${\\theta}$'),'Interpreter','latex','FontSize',18)
    xlim([min_t max_t])
    ylim([-max_th max_th])
    title(sprintf('${\\theta}$ at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('Sim','UKF estim','EKF estim')

%     subplot(5,1,2)
%     plot(t(1:k),omega(1:k),'--','color','#de425b',"Linewidth",3)
%     hold on 
%     plot(t(1:k),state_est(2,1:k),'color',colorUKF,"Linewidth",2)
%     hold on 
%     plot(t(1:k),omega_ekf(1:k),'color',colorEKF,"Linewidth",2)
%     xlabel('time t')
%     ylabel('\omega')
%     xlim([min_t max_t])
%     ylim([-max_th max_th])
%     title(sprintf('${\\omega}$ at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',20)
%     legend('Sim','UKF estim','EKF estim')
    
    %Figure for error
    subplot(3,1,2)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),abs(error1_ukf(1:k)),'color',colorUKF,"Linewidth",3)
    hold on
    plot(t(1:k),abs(error1_ekf(1:k)),'color',colorEKF,"Linewidth",3)
    xlabel('time t','FontSize',12)
    ylabel('Error','FontSize',12)
    xlim([min_t max_t])
    ylim([0 max_err])
    yline(0,'--');
%     title(sprintf('Error at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',20)
%     legend('Error UKF','Error EKF')
    title(sprintf('Error for ${\\theta}$ after t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('Error UKF \theta', 'Error EKF \theta','Location','southeast');

    %Figure for error covariance
    subplot(3,1,3)
    set(get(gca,'XAxis'),'FontSize', 12)
    set(get(gca,'YAxis'),'FontSize', 12)
    plot(t(1:k),sqrt(reshape(sigmas(1,1,1:k), 1,k)),'color',colorUKF,"Linewidth",3)
    hold on
    plot(t(1:k),sqrt(reshape(sigmas_ekf(1,1,1:k),1,k)),'color',colorEKF,"Linewidth",3)
    xlabel('time t','FontSize',12)
    ylabel(sprintf('${\\sigma}$'),'Interpreter','latex','FontSize',16)
    xlim([min_t max_t])
    ylim([0 max_cov])
    %ylim([-max_err max_err])
    yline(0,'--');
%     title(sprintf('Error std. dev. ${\\theta}$ at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',20)
%     legend('Error stand. dev. UKF','Error stand. dev. EKF')
    title(sprintf('Standard deviation ${\\theta}$ after t = %.2f seconds', t_k),'Interpreter','latex','FontSize',16)
    legend('St. dev. UKF \theta', 'St. dev. EKF \theta','Location','northeast')

    %Figure for error covariance
%     subplot(2,1,2)
%     figure(h4)
%     plot(t(1:k),sqrt(reshape(sigmas(2,2,1:k),1,k)),'color',colorUKF,"Linewidth",3)
%     hold on
%     plot(t(1:k),sqrt(reshape(sigmas_ekf(2,2,1:k),1,k)),'color',colorEKF,"Linewidth",3)
%     xlabel('time t')
%     ylabel('Error standard deviation')
%     xlim([min_t max_t])
%     %ylim([-max_err max_err])
%     yline(0,'--');
%     title(sprintf('Error std. dev. ${\\omega}$ at t = %.2f seconds', t_k),'Interpreter','latex','FontSize',20)
%     legend('Error stand. dev. UKF','Error stand. dev. EKF')

    
    %movieVector(k) = getframe(gcf);
    %drawnow
end

%% Uncomment block to create movie
% myWriter = VideoWriter('Pendulum','MPEG-4');
% myWriter.FrameRate= 20;
% 
% open(myWriter);
% writeVideo(myWriter, movieVector);
% close(myWriter);