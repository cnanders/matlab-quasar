close all
clc
clear

out = quasar(...
    'numPoles', uint8(2), ...
    'theta', 120, ...
    'numR', uint8(7) ...
);

x = out.x;
y = out.y;
r = out.r;
theta = out.theta;
t = out.t;

figure
subplot(121)
plot(r, '.-b')
title('r')

subplot(122)
plot(theta, '.-b')
title('theta')


figure
plot(x, y)
axis image
xlim([-1 1])
ylim([-1 1])


figure
hold on
plot(t, x, '.-r');
plot(t, y, '.-b')
legend({'x', 'y'})
xlabel('time');
ylabel('amp [arb]')

figure
plot3(x, y, t, '.-b')
xlabel('x')
ylabel('y')
zlabel('time (s)')
xlim([-1 1])
ylim([-1 1])