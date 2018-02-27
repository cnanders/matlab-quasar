# MATLAB Quasar Illumination Waveform Generator

Generates  x(t) y(t) continuous waveforms suitable for a tip/tilt scanner to “paint” an arbitrary quasar illumination pattern in time.  

# Examples

## Dipole

x(t) vs y(t)

![Quasar Dipole Illumination](img/dipole.jpg?raw=true)

x(t) and y(t)

![Quasar Dipole Illumination](img/dipole-waveforms.jpg?raw=true)

## Dipole With Different Paramaters

x(t) vs y(t)

![Quasar Dipole Illumination](img/dipole-2.jpg?raw=true)

x(t) and y(t)

![Quasar Dipole Illumination](img/dipole-2-waveforms.jpg?raw=true)

## Quadrupole

x(t) vs y(t)

![Quasar Quadrupole Illumination](img/quad-pole.jpg?raw=true)

x(t) and y(t)

![Quasar Quadrupole Illumination](img/quad-pole-waveforms.jpg?raw=true)

## 11 Poles 

This would probably never be used in practice but shows that `quasar` is generalized and works with any number of poles. 

x(t) vs y(t)

![Quasar 11-Pole Illumination](img/11-pole.jpg?raw=true)

x(t) and y(t)


## Out and Back Version (quasar2)

The `quasar2` function is an altered version that performs an “out and back” on every pole and does all pole transitions on the central connecting arcs.

quasar2 x(t) vs y(t)

![Quasar Quadrupole Illumination](img/quasar-2-4-pole.jpg?raw=true)

quasar2 x(t) vs y(t) vs t shows the trace in time

![Quasar Quadrupole Illumination](img/quasar2-4-pole.gif?raw=true)


# Documentation



The `quasar` function accepts six optional parameters in any order

```matlab
% @param {double 1x1} [radiusPoleInner = 0.5] - inner radius of poles [0 : 1]
% @param {double 1x1} [radiusPoleOuter = 0.7] - outer radius of poles [0 : 1]
% @param {uint8 1x1} [numArcs = 9] - number of arcs per pole (> 3 and odd)
% @param {double 1x1} [theta = 30] - angle subtended by each pole (deg) (> 0)
% @param {uint8 1x1} [numPoles = 4] - number of poles (> 0)
% @param {double 1x1} [dt = 10e-6] - time separation of samples (sec)
% @param {double 1x1} [period = 100e-3] - period of one full cycle (sec)
```

The `quasar` function returns a structure with the following properties:

```matlab
% @property {double 1xm} x - x values [-1 : 1]
% @property {double 1xm} y - y values [-1 : 1]
% @property {double 1xm} r - r values [0 : 1]
% @property {double 1xm} theta - theta values [0, 2pi]
% @property {double 1xm} t - time values
```
Parameters are supplied to `quasar` using the `varargin` syntax.  E.g.,

```matlab
out = quasar(...
    'numPoles', uint8(2), ...
    'numArcs', uint8(17) ...
    'theta', 120, ...
);
```

# Rotations, Offsets, and Scaling

Global rotations, offsets, and scaling can be applied to the data returned by `quasar`.  For example, a rotation of 45 degrees can be applied to a dipole quasar pattern.

```matlab
out = quasar(...
    'numPoles', uint8(2), ...
    'theta', 40, ...
    'numArcs', uint8(13) ...
);

x = out.x;
y = out.y;
r = out.r;
theta = out.theta;
t = out.t;

[theta, rho] = cart2pol(x, y)
theta = theta + pi/4;
[x, y] = pol2cart(theta, rho);
```
Dipole x(t) vs y(t) after 45-degree rotation

![Quasar Rotated Dipole Illumination](img/dipole-rotated-45.jpg?raw=true)

If the tip/tilt mirror is not at normal incidence, the amplitude of the resulting waveforms will need to be scaled by sin(AOI), for example.


# Connecting Paths

The connecting straight paths between the arcs of each a pole and the connecing curved paths between each pole are part of the output waveforms.  This enables smooth, continuous mechanical motion of the tip/tilt scanner  
