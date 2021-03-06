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

function [out] = quasar(varargin)
    
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
    length_radial_connectors = numPoles * (radiusPoleOuter - radiusPoleInner);

    % by subtracting the angle subtended by all of the poles from 360, the 
    % angle for all arc connectors remains.  
    angle_for_arc_connectors = 360 - numPoles * theta;

    % there are an even number of inner and outer arc connectors
    % This assumes numPoles is even.  Eventually, I
    % generalized this to work with numPoles odd as well.  When numPoles is odd
    % the last arc connector has to go from radiusPoleOuter to radiusPoleInner so its length is only
    % approximated by the math here.

    length_arc_connectors_inner = (angle_for_arc_connectors / 2) * pi / 180 * radiusPoleInner;
    length_arc_connectors_outer = (angle_for_arc_connectors / 2) * pi / 180 * radiusPoleOuter;

    length_arc_connector_inner = length_arc_connectors_inner / (numPoles / 2);
    length_arc_connector_outer = length_arc_connectors_outer / (numPoles / 2);

    length_arc_connectors = length_arc_connectors_inner + length_arc_connectors_outer;

    % list of radius values of the arcs within each pole
    r = linspace(radiusPoleInner, radiusPoleOuter, numArcs);

    % length of the {numArcs} archs of {numPoles} poles
    length_arcs = sum(theta_rad .* r) * numPoles;

    % length of the entire path
    length_period = length_arcs + length_arc_connectors + length_radial_connectors;


    % number of samples in each arc of a pole (orded by increasing radius)
    samples_of_arc =  theta_rad .* r / length_period * samples;

    samples_of_arc_connector_outer = length_arc_connector_outer / length_period * samples;
    samples_of_arc_connector_inner = length_arc_connector_inner / length_period * samples;

    samples_of_radial_connector = (radiusPoleOuter - radiusPoleInner) / (numArcs - 1) / length_period * samples;

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

        if mod(n, 2) == 0
            % even pole (out to in)

            for m = length(r) : -1 : 1

                theta_line = linspace(theta_start, theta_end, round(samples_of_arc(m))); % counter clockwise
                if mod(m, 2) == 0 % even === clockwise
                    theta_line = fliplr(theta_line);
                end
                r_line = r(m) * ones(size(theta_line));

                r_out = [r_out, r_line];
                theta_out = [theta_out, theta_line];

                % Add a connector from out to in at the last theta value

                if (m > 1)
                    r_line = linspace(r(m), r(m - 1), round(samples_of_radial_connector));
                    theta_line = theta_line(end) * ones(size(r_line));
                    r_out = [r_out, r_line];
                    theta_out = [theta_out, theta_line];
                end


            end

            % add arc connector to next pole
            % RECALL thetaPoleSep is angle between the center of adjacent
            % poles (360 / numPoles)

            theta_center = thetaPoleSep/2 + (n - 1) * thetaPoleSep;
            theta_start = theta_center - (thetaPoleSep - theta) / 2;
            theta_end = theta_center + (thetaPoleSep - theta) / 2;
            theta_line = linspace(theta_start, theta_end, round(samples_of_arc_connector_inner));
            r_line = r(m) * ones(size(theta_line));

            r_out = [r_out, r_line];
            theta_out = [theta_out, theta_line];

        % end even pole (out to in)
        else
            % odd pole in to out
            for m = 1 : length(r)

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

            % SPECIAL CASE when numPoles === 1, copy a fliplr version of r_line
            % and theta_line to r_out and theta_out

            if numPoles == 1

                r_out = [r_out, fliplr(r_out)]
                theta_out = [theta_out, fliplr(theta_out)]
            else

                theta_center = thetaPoleSep/2 + (n - 1) * thetaPoleSep;
                theta_start = theta_center - (thetaPoleSep - theta) / 2;
                theta_end = theta_center + (thetaPoleSep - theta) / 2;
                theta_line = linspace(theta_start, theta_end, round(samples_of_arc_connector_outer));

                if n == numPoles
                    % last arc connector needs to radially move from out to in
                    % since there is an even number of poles and we ended on the 
                    % outside but need to get back to the inside
                    r_line = linspace(radiusPoleOuter, radiusPoleInner, length(theta_line));
                else

                    r_line = r(m) * ones(size(theta_line));
                end

                r_out = [r_out, r_line];
                theta_out = [theta_out, theta_line];

            end

        end % odd pole in to out
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