% @typedef {struct 1x1} quasar_data
% @property {double 1xm} x - x values [-1 : 1]
% @property {double 1xm} y - y values [-1 : 1]
% @property {double 1xm} r - r values [0 : 1]
% @property {double 1xm} theta - theta values [0, 2pi]
% @property {double 1xm} t - time values

% @param {double 1x1} [radiusPoleInner = 0.5] - inner radius of poles [0 : 1]
% @param {double 1x1} [radiusPoleOuter = 0.7] - outer radius of poles [0 : 1]
% @param {uint8 1x1} [numArcs = 9] - number of arcs per pole (> 3 and odd)
% @param {double 1x1} [theta = 30] - angle subtended by each pole (deg) (>
% 0)
% @param {uint8 1x1} [numPoles = 4] - number of poles (> 0)
% @param {double 1x1} [dt = 10e-6] - separation of time samples (sec)
% @param {double 1x1} [period = 100e-3] - period of 1 full cycle (sec)
% @return {quasar_data 1x1}

% This differs from quasar (original) because each pole goes out and back
% and all transitions occur in the ctner

function [out] = quasar2(varargin)
    
    %% Input Validation and Parsing
    
    p = inputParser;
    
    iseven = @(x) mod(x, 2) == 0;
    isodd = @(x) mod(x, 2) ~= 0;

    addParameter(p, 'radiusPoleInner', 0.5, @(x) isscalar(x) && isnumeric(x) && (x > 0) && (x <= 1))
    addParameter(p, 'radiusPoleOuter', 0.7, @(x) isscalar(x) && isnumeric(x) && (x > 0) && (x <= 1))
    addParameter(p, 'numArcs', 9, @(x) isscalar(x) && isinteger(x) && (x >= 3) && isodd(x))
    addParameter(p, 'numPoles', 4, @(x) isscalar(x) && isinteger(x) && (x > 0))
    addParameter(p, 'theta', 30, @(x) isscalar(x) && isnumeric(x) && (x > 0))
    addParameter(p, 'dt', 10e-6, @(x) isscalar(x) && isnumeric(x) && (x > 0))
    addParameter(p, 'period', 100e-3, @(x) isscalar(x) && isnumeric(x) && (x > 0))

    parse(p, varargin{:});

    radiusPoleInner = p.Results.radiusPoleInner;
    radiusPoleOuter = p.Results.radiusPoleOuter;
    numArcs = double(p.Results.numArcs);
    theta = p.Results.theta;
    numPoles = double(p.Results.numPoles);
    dt = p.Results.dt;
    period = p.Results.period;
    
    
    %% Begin
    
    samples = period / dt;
    theta_rad = theta * pi / 180;  
        

    % the general strategy is to compute the length of the drawn path assuming
    % that the circle has a radius of 1.  The path is comprised of:
    % 1)    arcs within each pole
    % 2)    radial connectors between arcs within each pole
    % 3)    arcs that connect poles
    % once the total length of the drawn path is known, the number of 
    % samples in each subpath can be determined.  There are probably other
    % more elegant ways to build this code involveing reflections and such
    % but I didn't care about elegance; I wanted to get it working quickly.

    % each pole has numArcs - 1 radial connectors
    % their combined length is radiusPoleOuter - radiusPoleInner
    length_radial_connectors = numPoles * (radiusPoleOuter - radiusPoleInner) * 2;
     % there are (numArcs - 1) * 2 radial connectors per pole (out and back)
    length_radial_connector = length_radial_connectors / ((numArcs - 1) * 2 * numPoles);
    
    
    % by subtracting the angle subtended by all of the poles from 360, the 
    % angle for all arc connectors remains.  
    angle_for_arc_connectors = 360 - numPoles * theta; % degrees

    % there are an even number of inner and outer arc connectors
    % This assumes numPoles is even.  Eventually, I
    % generalized this to work with numPoles odd as well.  When numPoles is odd
    % the last arc connector has to go from radiusPoleOuter to radiusPoleInner so its length is only
    % approximated by the math here.

    length_arc_connectors = angle_for_arc_connectors * pi / 180 * radiusPoleInner;
    length_arc_connector = length_arc_connectors / numPoles;

    % list of radius values of the arcs within each pole
    r = linspace(radiusPoleInner, radiusPoleOuter, numArcs);
    r = [r fliplr(r(1 : end - 1))];
    

    % length of the {numArcs} archs of {numPoles} poles
    length_arcs = sum(theta_rad .* r) * numPoles;

    % length of the entire path
    length_period = length_arcs + length_arc_connectors + length_radial_connectors;


    % number of samples in each arc of a pole (orded by increasing radius)
    samples_of_arc =  theta_rad .* r / length_period * samples;
    samples_of_arc_connector = length_arc_connector / length_period * samples;
    % there are (numArcs - 1) * 2 radial connectors per pole (out and back)
    
    samples_of_radial_connector = length_radial_connector / length_period * samples;

    % the above samples_* variables are not integers.  In the code below I
    % round to integer numbers.  Where different lines join,
    % there will be some small blips but this will be smoothed after the
    % butterworth filter is applied.

    % Build everything in polar coordinates since that is easier to reason
    % about then transform to cartesian at the end.
    
    % Storage for (r, theta) samples along each line segment of the total
    % path.  
    r_out = [];
    theta_out = [];

    thetaPoleSep = 360 / numPoles;
    
    for n = 1 : numPoles % poles

        theta_center = (n - 1) * thetaPoleSep;
        theta_start = theta_center - theta / 2;
        theta_end = theta_center + theta / 2;
        
        for m = 1 : length(r)

            % special case, if monopole, don't draw last arc
            
            if numPoles == 1 && m == length(r)
                continue
            end
            
            theta_line = linspace(theta_start, theta_end, round(samples_of_arc(m)));
            
            if mod(m, 2) == 0 % even
                theta_line = fliplr(theta_line);
            end
            r_line = r(m) * ones(size(theta_line));

            r_out = [r_out, r_line];
            theta_out = [theta_out, theta_line];

            % Add a connector from in to out at the last theta value

            if (m < length(r))
                r_line = linspace(r(m), r(m + 1), round(samples_of_radial_connector));
                theta_line = theta_line(end) * ones(size(r_line));
                r_out = [r_out, r_line];
                theta_out = [theta_out, theta_line];
            end

        end
        
        % add arc connector to next pole
        % RECALL thetaPoleSep is angle between the center of adjacent
        % poles (360 / numPoles)

        % Special case for monopole, no connector
        if numPoles == 1
            continue
        end
        
        theta_center = thetaPoleSep/2 + (n - 1) * thetaPoleSep;
        theta_start = theta_center - (thetaPoleSep - theta) / 2;
        theta_end = theta_center + (thetaPoleSep - theta) / 2;
        theta_line = linspace(theta_start, theta_end, round(samples_of_arc_connector));
        r_line = r(m) * ones(size(theta_line));

        r_out = [r_out, r_line];
        theta_out = [theta_out, theta_line];
    end

    % theta_out = theta_out + 45;
    
    t = 0 : dt : (length(r_out) - 1) * dt;
    x = r_out .* cos(theta_out * pi / 180);
    y = r_out .* sin(theta_out * pi / 180);

    out = struct();
    out.x = x;
    out.y = y;
    out.r = r_out;
    out.theta = theta_out * pi / 180;
    out.t = t;

end