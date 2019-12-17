function [time, support_y] = min_time(support_y1, support_y2)
for i = 1:length(support_y1)
    if support_y1(i) < -5
        support_y1(i) = -5;
    elseif support_y1(i) > 10
        support_y1(i) = 10;
    else 
        support_y1(i) = support_y1(i);
    end 
end 

for i = 1:length(support_y2)
    if support_y2(i) < -5
        support_y2(i) = -5;
    elseif support_y2(i) > 10
        support_y2(i) = 10;
    else 
        support_y2(i) = support_y2(i);
    end 
end 
% Define the locations of the cities.
support_x = 0:2:40;
support_y = [0, support_y1, 3, support_y2, 0];
% Set the car and motor parameters
m = 100;  % kg
time_span_m = [0.67, 1.67];  % 1 second, at the star
Fm = 800;  % Newtons
v0 = 0;

% Plot the track, colouring the moveable and immovable supports
% hold on;
% plot(support_x, support_y, '.b', 'markersize', 15);
% plot(support_x([1,11,21]), support_y([1,11,21]), '.r', 'markersize', 20)
pp = spline(support_x, support_y);
xx = linspace(0, max(support_x), 200);
yy = ppval(pp, xx);
% plot(xx, yy, 'k');

% Gravity
g = -9.8;

% initial position and velocities
x0 = support_x(1);
vx0 = v0/(sqrt(1 + splineDeriv1(pp, 0)^2));
y0 = support_y(1);
vy0 = v0*splineDeriv1(pp, 0)/(sqrt(1 + splineDeriv1(pp, 0)^2));

% DE function
f = @(t, w) rollerCoasterDE(t, w, pp, m, Fm, time_span_m);
time_span = 0:0.05:10;
w0 = [x0, vx0, y0, vy0];

threshold = 40;
ev = @(t, y) events_func(t, y, threshold);
options = odeset('AbsTol', 1e-10, 'RelTol', 1e-8, 'events', ev);
[t, w] = ode45(f, time_span, w0, options);

support_y

time = t(end);