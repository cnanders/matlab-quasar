close all
clc
clear

% add src to path
[path, name, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(path, '..', 'src')))

out = quasar2(...
    'numPoles', uint8(11), ...
    'theta', 10 ... 
);

x = out.x;
y = out.y;
r = out.r;
theta = out.theta;
t = out.t;

figure('Color', 'white')
subplot(121)
plot(r, '.-b')
title('r')

subplot(122)
plot(theta, '.-b')
title('theta')


figure('Color', 'white')
plot(x, y, '.-b')
xlabel('x')
ylabel('y')
axis image
xlim([-1 1])
ylim([-1 1])


figure('Color', 'white')
hold on
plot(t, x, '.-r');
plot(t, y, '.-b')
legend({'x', 'y'})
xlabel('time');
ylabel('amp [arb]')

figure('Color', 'white')
plot3(x, y, t, '.-b')
xlabel('x')
ylabel('y')
zlabel('time (s)')
xlim([-1 1])
ylim([-1 1])