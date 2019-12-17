close all;
clear;
clc;

support_y1 = [ -3 -5 0 ]; 
support_y2 = [ -2.68 -1 0];

y0 = [ 0 -3 -5 0 5 -2.68 -1 0 0];
v = [ support_y1;
    support_y2];
fv = @(v) -1*max_distance(v(1,:),v(2,:));
options = optimset('Display','iter');
x0 = [-1,2];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-5.01 -5.01 -5.01 -5.01 -5.01 -5.01];
ub = [10 10 10 10 10 10];
nonlcon = [];
topt = fmincon(fv,v,A,b,Aeq,beq,lb,ub,nonlcon,options);

% figure('OuterPosition', [1 1 1200 800]);

% *********************** DT1 ***************************
% Set the car and motor parameters
% v0 = 0;  % m/s
% m = 100;  % kg
% Fm = 800;  % Newtons
% time_span_m = [0.8, 1.8];  % 1 second, at the start
% v0 = 0;  % initial speed = 0
% 
% % Load the% Plot the track
% % clf; hold on;
% pp = spline(support_x, support_y);
% xx = linspace(0, max(support_x), 200);
% yy = ppval(pp, xx);
% % plot(xx, yy, 'k', 'linewidth', 3);
% % plot(x, y, '.b', 'markersize', 20);
% 
% % Gravity
% g = -9.8;
% 
% 
% % initial position and velocities
% x0 = support_x(1);
% vx0 = v0/(sqrt(1 + splineDeriv1(pp, 0)^2));
% y0 = support_y(1);
% vy0 = v0*splineDeriv1(pp, 0)/(sqrt(1 + splineDeriv1(pp, 0)^2));
% 
% % DE function
% f = @(t, w) rollerCoasterDE(t, w, pp, m, Fm, time_span_m);
% time_span = 0:0.05:5;
% w0 = [x0, vx0, y0, vy0];
% 
% threshold = 10;
% ev = @(t, y) events_func(t, y, threshold);
% options = odeset('AbsTol', 1e-10, 'RelTol', 1e-8, 'events', ev);
% [t, w] = ode45(f, time_span, w0, options); support locations
% A = csvread('DesignTask.csv');
% support_x = A(:,1);  % first column in the file
% support_y = A(:,2);  % second column in the file
% 
% 
% 
% %% Position and Velocities
% 
% xpos = w(:,1);
% xvel = w(:,2);
% ypos = w(:,3);
% yvel = w(:,4);
% 
% xend = xpos(end);
% xvel_end = xvel(end);
% yend = ypos(end);
% yvel_end = yvel(end);
% 
% mag_vel = norm([xvel_end yvel_end]);
% theta = atand(yvel_end/xvel_end);
% 
% %% Trajectory
% t2 = (0:0.05:3);
% 
% x2 = xvel_end.*t2 + xend;
% y2 = yvel_end.*t2 + yend + 0.5.*g.*(t2.^2);
% 
% xvel_t = zeros(1, length(x2));
% xvel_t(1,:) = xvel_end;
% 
% yvel_t = zeros(1, length(x2));
% yvel_t(1,:) = yvel_end + g.*t2;
% 
% landing = landing_x(yend, mag_vel, theta) + xend;
% 
% for (i = 1:length(x2))
%     if (y2(i)>=0)
%         x_out(i,1) = x2(i);
%         y_out(i,1) = y2(i);
%     end
% end
% 
% % figure('OuterPosition', [0 1 1200 800 ])
% % for (i = 1:length(xpos))
% %     % Plot the specific point on the animation
% %     plot(xx, yy, 'b');
% %     hold on;
% %     plot(w(:,1), w(:,3), 'r');
% %     plot(x_out, y_out, 'g');
% %     plot(x_out(end), y_out(end), 'xr', 'markersize', 20);
% %     ylim([-10 17]);
% %     xlim([-2 30]);
% %     
% %     if (t(i)>time_span_m(1) && t(i)<time_span_m(2))
% %         plot(xpos(i), ypos(i), 's', 'markersize', 10, 'MarkerFaceColor', 'blue');
% %         hold on;
% %         quiver(xpos(i), ypos(i), xvel(i), yvel(i));
% %         text(1, 1, sprintf('Time: %.2f s', t(i)) );
% %         text(1, 0, sprintf('X Velocity: %.2f m/s', xvel(i)));
% %         text(1, -1, sprintf('Y Velocity: %.2f m/s', yvel(i)) );
% %         title(sprintf('Motor: ON'));
% %         getframe;
% %     elseif (i < length(xpos))
% %         plot(xpos(i), ypos(i), 's', 'markersize', 10, 'MarkerFaceColor', 'blue');
% %         text(1, 1, sprintf('Time: %.2f s', t(i)) );
% %         text(1, 0, sprintf('X Velocity: %.2f m/s', xvel(i)));
% %         text(1, -1, sprintf('Y Velocity: %.2f m/s', yvel(i)) );
% %         title(sprintf('Motor: OFF'));
% %         getframe;
% %     else
% %         for j = (1:length(x_out))
% %             plot(xx, yy, 'b');
% %             hold on;
% %             title(sprintf('Take off: %0.2f m/s at %0.1f degrees', mag_vel, theta));
% %             plot(w(:,1), w(:,3), 'r');
% %             plot(x_out, y_out, 'g');
% %             plot(x_out(end), y_out(end), 'xr', 'markersize', 20);
% %             plot(x_out(j), y_out(j), 's', 'markersize', 10, 'MarkerFaceColor', 'blue');
% %             ylim([-10 17]);
% %             xlim([-2 30]);
% %             text(1, 1, sprintf('Time: %.2f s', t2(j) + t(end) ));
% %             text(1, 0, sprintf('X Velocity: %.2f m/s', xvel_t(j)));
% %             text(1, -1, sprintf('Y Velocity: %.2f m/s', yvel_t(j)) );
% %             getframe;
% %             hold off;
% %         end
% %         title(sprintf('Impact Location: %.2f m', landing));
% %     end
% %     hold off;
% % end