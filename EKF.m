function [mu,Sigma] = EKF(mu_prev,Sigma_prev,z,R,Q,dt,m,L)
% EKF This function uses Extended Kalman Filter algorithm to estimate
% states given previous state, uncertainty, process and measurement noise
%   Input arguments:
%       mu_prev is ndim x 1 vector consisting of previous time step's state
%       Sigma_prev is ndim x ndim uncertainty matrix for previous time step
%       z is ndim x 1 vector with measurements
%       R is ndim x ndim process noise matrix
%       Q is ndim x ndim measurement noise matrix
%       ndim is number of state that is estimated, 2 for simple and 4 for
%       double
%       dt is timestep
%       m is mass. Will be a vector for double pendulum with m1 and m2
%       L is length of the bar(s). Will be a vector for double pendulum
%       with L1 and L2
%   Output arguments:
%       mu is ndim x 1 vector consisting of estimated states
%       Sigma is ndim x ndim uncertainty matrix

g = 9.82;

if length(L) == 1   % Runs for single pendulum

    %Jacobian of motion model, SIMPLE pendulum:
    jacob = @(theta,omega) -9.82*dt/L;
    G = [1 dt;
        jacob(mu_prev(1),mu_prev(2)) 1];

    H = [1 0];

elseif length(L) == 2   % Runs for double pendulum

    L_1 = L(1);
    L_2 = L(2);
    m_1 = m(1);
    m_2 = m(2);
    %Jacobian of motion model, DOUBLE pendulum:
    
    jacob_11 = mu_prev(3);
    jacob_12 = 0;
    jacob_13 = dt;
    jacob_14 = 0;
    
    jacob_21 = 0;
    jacob_22 = mu_prev(4);
    jacob_23 = 0;
    jacob_24 = dt;
    
    
    jacob_31 = @(theta_1, theta_2,omega_1, omega_2)((2*L_1*m_2*omega_1^2*sin(theta_1-theta_2)^2-2*m_2*cos(theta_1-theta_2)*(L_1*omega_1^2*cos(theta_1-theta_2)+L_2*omega_2^2)...
        -g*m_2*cos(theta_1-2*theta_2)-g*(m_2+2*m_1)*cos(theta_1))*(-m_2*cos(2*theta_1-2*theta_2)+m_2+2*m_1)-2*m_2*(-2*m_2*(L_1*omega_1^2*cos(theta_1-theta_2)+L_2*omega_2^2)*sin(theta_1-theta_2)...
        -g*m_2*sin(theta_1-2*theta_2)-g*(m_2+2*m_1)*sin(theta_1))*sin(2*theta_1-2*theta_2))...
        /(L_1*(-m_2*cos(2*theta_1-2*theta_2)+m_2+2*m_1)^2);
    
    jacob_32 = @(theta_1, theta_2,omega_1, omega_2)((m_2*(cos(2*(theta_2-theta_1))-1)-2*m_1)*(-2*g*m_2*cos(2*theta_2-theta_1)+2*L_1*m_2*omega_1^2*sin(theta_2-theta_1)^2-2*m_2*cos(theta_2-theta_1)...
        *(L_1*omega_1^2*cos(theta_2-theta_1)+L_2*omega_2^2))-2*m_2*sin(2*(theta_2-theta_1))*(g*m_2*sin(2*theta_2-theta_1)+2*m_2*(L_1*omega_1^2*cos(theta_2-theta_1)+L_2*omega_2^2)...
        *sin(theta_2-theta_1)-g*(m_2+2*m_1)*sin(theta_1)))...
        /(L_1*(m_2*(cos(2*(theta_2-theta_1))-1)-2*m_1)^2);
    jacob_33 = @(theta_1, theta_2,omega_1, omega_2) -(2*m_2*sin(2*(theta_2-theta_1))*omega_1)/(m_2*(cos(2*(theta_2-theta_1))-1)-2*m_1);
    jacob_34 = @(theta_1, theta_2,omega_1, omega_2) -(4*L_2*m_2*sin(theta_2-theta_1)*omega_2)/(L_1*(m_2*(cos(2*(theta_2-theta_1))-1)-2*m_1));
    
    jacob_41 = @(theta_1, theta_2,omega_1, omega_2) ((-2*sin(theta_1-theta_2)*(-L_2*m_2*omega_2^2*sin(theta_1-theta_2)-g*(m_2+m_1)*sin(theta_1))-2*cos(theta_1-theta_2)*(L_2*m_2*omega_2^2*cos(theta_1-theta_2)...
        +g*(m_2+m_1)*cos(theta_1)+L_1*(m_2+m_1)*omega_1^2))*(m_2*(cos(2*(theta_1-theta_2))-1)-2*m_1)-4*m_2*(L_2*m_2*omega_2^2*cos(theta_1-theta_2)+g*(m_2+m_1)*cos(theta_1)+L_1*(m_2+m_1)*omega_1^2)*sin(theta_1-theta_2)*sin(2*(theta_1-theta_2)))...
        /(L_2*(m_2*(cos(2*(theta_1-theta_2))-1)-2*m_1)^2);
    
    jacob_42 = @(theta_1, theta_2,omega_1, omega_2) (4*m_2*(L_2*m_2*omega_2^2*cos(theta_2-theta_1)+g*(m_2+m_1)*cos(theta_1)+L_1*(m_2+m_1)*omega_1^2)*sin(theta_2-theta_1)*sin(2*(theta_2-theta_1))...
        +(2*cos(theta_2-theta_1)*(L_2*m_2*omega_2^2*cos(theta_2-theta_1)+g*(m_2+m_1)*cos(theta_1)+L_1*(m_2+m_1)*omega_1^2)-2*L_2*m_2*omega_2^2*sin(theta_2-theta_1)^2)*(m_2*(cos(2*(theta_2-theta_1))-1)-2*m_1))...
        /(L_2*(m_2*(cos(2*(theta_2-theta_1))-1)-2*m_1)^2);
    
    jacob_43 = @(theta_1, theta_2,omega_1, omega_2) -(4*L_1*(m_2+m_1)*sin(theta_2-theta_1)*omega_1)/(L_2*(-m_2*cos(2*theta_2-2*theta_1)+m_2+2*m_1));
    jacob_44 = @(theta_1, theta_2,omega_1, omega_2) -(4*m_2*cos(theta_2-theta_1)*sin(theta_2-theta_1)*omega_2)/(-m_2*cos(2*theta_2-2*theta_1)+m_2+2*m_1);
    
    G = [jacob_11, jacob_12, jacob_13, jacob_14 ;
         jacob_21, jacob_22, jacob_23, jacob_24 ;
         jacob_31(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt, jacob_32(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt, jacob_33(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt + 1, jacob_34(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt ;
         jacob_41(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt, jacob_42(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt, jacob_43(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt, jacob_44(mu_prev(1),mu_prev(2),mu_prev(3),mu_prev(4))*dt + 1 ];

%     H = [1 0 0 0;
%          0 1 0 0;
%          0 0 1 0;
%          0 0 0 1];
    H = [1 0 0 0;
         0 1 0 0];
end

%prediction step
mu_bar = state_trans(mu_prev,m,L,dt) + 1*normrnd(0,0.1,size(mu_prev,1),1); %+ 0*normrnd(0,0.2,2,1); %TODO make it adapt to pendulum
Sigma_bar = G * Sigma_prev * G' + R;

%update step (measurement model is linear, so H = C)
z_bar = H*mu_bar; %measurment with predicted state
K = Sigma_bar*H'/(H*Sigma_bar*H' + Q); %Kalman gain

mu = mu_bar + K*(z - z_bar);
Sigma = (eye(size(Sigma_prev)) - K*H)*Sigma_bar;


end