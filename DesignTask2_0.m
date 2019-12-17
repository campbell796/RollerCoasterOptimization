close all;
clear;
clc;


% % Define a sample track with 21 supports
support_x = 0:1:20;
support_y1 = [-2	-3.25	-4.25	-4.75	-5	-4.5	-3.5	-2	1];
support_y2 = [2.605991303	1.514476327	-0.023803668	-1.655477203	-3.005718566	-3.764352242	-3.757105988	-2.985644404	-1.62718654];

% Define a sample track with 21 supports
% A = xlsread('DesignTask2.xlsx');
% support_x = A(:,1);  % first column in the file
% support_y = A(:,2);  % second column in the file
% 
y0 = [0	-2	-3.25	-4.25	-4.75	-5	-4.5	-3.5	-2	1	3	2.605991303	1.514476327	-0.023803668	-1.655477203	-3.005718566	-3.764352242	-3.757105988	-2.985644404	-1.62718654	0];
v = [support_y1;
    support_y2]
fv = @(v) min_time(v(1,:), v(2,:));
options = optimset('Display', 'iter');
% lb = [-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5];
% up = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10];
% topt = fmincon(fv, v, [], [], [], [], lb, up, [], options);
topt = fminsearch(fv,v,options);
% 
% % Set the car and motor parameters
% speed0 = 0;  % m/s
% m = 100;  % kg
% time_span_m = [0.67, 1.67];  % 1 second, at the star
% Fm = 800;  % Newtons
% v0 = 0;
% 
% % Plot the track, colouring the moveable and immovable supports
% % hold on;
% % plot(support_x, support_y, '.b', 'markersize', 15);
% % plot(support_x([1,11,21]), support_y([1,11,21]), '.r', 'markersize', 20)
% pp = spline(support_x, support_y);
% xx = linspace(0, max(support_x), 200);
% yy = ppval(pp, xx);
% % plot(xx, yy, 'k');
% 
% % Gravity
% g = -9.8;
% 
% % initial position and velocities
% x0 = support_x(1);
% vx0 = v0/(sqrt(1 + splineDeriv1(pp, 0)^2));
% y0 = support_y(1);
% vy0 = v0*splineDeriv1(pp, 0)/(sqrt(1 + splineDeriv1(pp, 0)^2));
% 
% % DE function
% f = @(t, w) rollerCoasterDE(t, w, pp, m, Fm, time_span_m);
% time_span = 0:0.05:10;
% w0 = [x0, vx0, y0, vy0];
% 
% threshold = 40;
% ev = @(t, y) events_func(t, y, threshold);
% options = odeset('AbsTol', 1e-10, 'RelTol', 1e-8, 'events', ev);
% [t, w] = ode45(f, time_span, w0, options);
% 
% % Position and Velocities
% 
% xpos = w(:,1);
% xvel = w(:,2);
% ypos = w(:,3);
% yvel = w(:,4);
% 
% 
% plot(support_x, support_y, '.r');
% hold on;
% plot(xx, yy);
% 
% figure('OuterPosition', [0 1 1200 800 ])
% for (i = 1:length(xpos))
% %     Plot the specific point on the animation
%     plot(xx, yy, 'b');
%     hold on;
%     plot(w(:,1), w(:,3), 'r');
%     ylim([-10 10]);
%     xlim([-2 45]);
%     
%     if (t(i)>time_span_m(1) && t(i)<time_span_m(2))
%         plot(xpos(i), ypos(i), 's', 'markersize', 10, 'MarkerFaceColor', 'blue');
%         hold on;
%         quiver(xpos(i), ypos(i), xvel(i), yvel(i));
%         text(1, 1, sprintf('Time: %.2f s', t(i)) );
%         text(1, 0, sprintf('X Velocity: %.2f m/s', xvel(i)));
%         text(1, -1, sprintf('Y Velocity: %.2f m/s', yvel(i)) );
%         title(sprintf('Motor: ON'));
%         getframe;
%     else (i < length(xpos))
%         plot(xpos(i), ypos(i), 's', 'markersize', 10, 'MarkerFaceColor', 'blue');
%         text(1, 1, sprintf('Time: %.2f s', t(i)) );
%         text(1, 0, sprintf('X Velocity: %.2f m/s', xvel(i)));
%         text(1, -1, sprintf('Y Velocity: %.2f m/s', yvel(i)) );
%         title(sprintf('Motor: OFF'));
%         getframe;
%     end
%     hold off;
% end
% 
% time = t(end)