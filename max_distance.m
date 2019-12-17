function [dist ] = max_distance(support_y1, support_y2)
support_x = [ 0 1.25 2.5 3.75 5 6.25 7.5 8.75 10];
% Set the car and motor parameters
support_y = [0,support_y1,5,support_y2,0];

v0 = 0;  % m/s
m = 100;  % kg
Fm = 800;  % Newtons
time_span_m = [0.8, 1.8];  % 1 second, at the start
v0 = 0;  % initial speed = 0

% Load the% Plot the track
% clf; hold on;
pp = spline(support_x, support_y);
xx = linspace(0, max(support_x), 200);
yy = ppval(pp, xx);
% plot(xx, yy, 'k', 'linewidth', 3);
% plot(x, y, '.b', 'markersize', 20);

% Gravity
g = -9.8;


% initial position and velocities
x0 = support_x(1);
vx0 = v0/(sqrt(1 + splineDeriv1(pp, 0)^2));
y0 = support_y(1);
vy0 = v0*splineDeriv1(pp, 0)/(sqrt(1 + splineDeriv1(pp, 0)^2));

% DE function
f = @(t, w) rollerCoasterDE(t, w, pp, m, Fm, time_span_m);
time_span = 0:0.05:5;
w0 = [x0, vx0, y0, vy0];

threshold = 10;
ev = @(t, y) events_func(t, y, threshold);
options = odeset('AbsTol', 1e-10, 'RelTol', 1e-8, 'events', ev);
[t, w] = ode45(f, time_span, w0, options);




%% Position and Velocities

xpos = w(:,1);
xvel = w(:,2);
ypos = w(:,3);
yvel = w(:,4);

xend = xpos(end);
xvel_end = xvel(end);
yend = ypos(end);
yvel_end = yvel(end);

mag_vel = norm([xvel_end yvel_end]);
theta = atand(yvel_end/xvel_end);

%% Trajectory
t2 = (0:0.05:3);

x2 = xvel_end.*t2 + xend;
y2 = yvel_end.*t2 + yend + 0.5.*g.*(t2.^2);

xvel_t = zeros(1, length(x2));
xvel_t(1,:) = xvel_end;

yvel_t = zeros(1, length(x2));
yvel_t(1,:) = yvel_end + g.*t2;

landing = landing_x(yend, mag_vel, theta) + xend;

for (i = 1:length(x2))
    if (y2(i)>=0)
        x_out(i,1) = x2(i);
        y_out(i,1) = y2(i);
    end
end
support_y
dist = landing;
