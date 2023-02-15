function [X] = state_trans(X_prev,m,L,dt)
%STATE_TRANS This function should be the state transistion function (model)
% where the mass in the pendulum goes from X_t-1 to X_t.
%   Input arguments:
%         X_prev is 2x1 vector consisting of angle and angular velocity for
%         the mass in the simple pendulum for time t-1
%   Output arguments:
%         X_t is 2x1 vector consisting of angle and angular velocity for
%         the mass in the simple pendulum for time t

g = 9.82;

if length(L) == 1   % Runs for single pendulum

    % Model assume no drag, TODO: input model to the function so params can be changed in one file
    accel = @(theta,omega) -g*sin(theta)/L; 

    % theta = theta_prev + dt*omega
    % omega = omega_prega + dt*accel Simple step in time with Euler forward
    X = zeros(size(X_prev));
    A = [1 dt;
         0 1];
    B = [0;
         dt];
    for i=1:size(X_prev,2)
        u = accel(X_prev(1,i),X_prev(2,i));
        X(:,i) = A*X_prev(:,i)+B*u;
    end

elseif length(L) == 2   % Runs for double pendulum
    X = zeros(size(X_prev));
    L1 = L(1);
    L2 = L(2);
    m1 = m(1);
    m2 = m(2);

    % Acceleration on mass 1 respectively mass 2 for double pendulum
    accel1 = @(theta1,theta2,omega1,omega2) (-g*(2*m1+m2)*sin(theta1) - m2*g*sin(theta1-2*theta2)...
            -2*sin(theta1-theta2)*m2*(omega2^2*L2 + omega1^2*L1*cos(theta1 - theta2)))...
            / (L1*(2*m1+m2 - m2*cos(2*theta1-2*theta2)));
    accel2 = @(theta1,theta2,omega1,omega2) (2*sin(theta1 - theta2)*(omega1^2*L1*(m1+m2)+g*(m1+m2)*cos(theta1)...
            +omega2^2*L2*m2*cos(theta1 - theta2)))/ (L2*(2*m1+m2-m2*cos(2*theta1-2*theta2)));


    A = [1 0 dt 0;
         0 1 0 dt;
         0 0 1 0;
         0 0 0 1];
    B = [0 0;
         0 0;
         dt 0;
         0 dt];
    for i=1:size(X_prev,2)
        u = [accel1(X_prev(1,i),X_prev(2,i),X_prev(3,i),X_prev(4,i));
             accel2(X_prev(1,i),X_prev(2,i),X_prev(3,i),X_prev(4,i))];
        X(:,i) = A*X_prev(:,i)+B*u;
    end
end


end

