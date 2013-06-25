classdef ReactionDiffusion < handle
    properties
        xlim
        h
        diffusion
        y0
        boundary
        kinetics_fcn
        Tspan
        solver
        animation_speedup
        animation_framerate
        varnames
    end
    properties (SetAccess = private)
        x
        xa
        D
        odeopts
        odesol
        m
        n
        JPattern
    end
    
    methods
        
        function obj = ReactionDiffusion(varargin)
        %REACTIONDIFFUSION  Set up a reaction diffusion simulation
        %    Usage:
        %        SIM = REACTIONDIFFUSION(param1, val1, param2, val2...)
        %    Parameters required for simulation:
        %
        %        'kinetics_fun'    Function handle to a function of the
        %                          form f(t, x, y). If the kinetics system 
        %                          has m equations, then inputs will have t
        %                          scalar, x an n-element row vector and y
        %                          an m x n matrix, one column per spatial
        %                          location
        %        'diffusion'       A vector of diffusion coefficients, or a
        %                          cell vector of a mixture of scalars and
        %                          function handles. e.g
        %                          diffusion = {1, @(x) x.^2}
        %                          has a diffusion coefficient of 1 for the
        %                          first variable and x for the second. If
        %                          the entry is a function of x, it needs
        %                          to be vectorised (callable with a
        %                          vector).
        %        'y0'              A vector of initial conditions. This
        %                          the same format as 'diffusion', either
        %                          scalars or functions of x
        %    Optional parameters:
        %
        %        'Tspan'           Two-element vector, default [0 10]
        %        'xlim'            Two-element vector, default [0 1]
        %        'h'               Spatial step size, default 0.01
        %        'boundary'        {|'noflux'|, 'periodic}
        %        'animation_speedup' Speedup factor for the animation. E.g.
        %                          if speedup is 2, 2 seconds of simulation
        %                          time will display in 1 second. Default 1
        %        'animation_framerate' Default 20
        %        'varnames'        Cell array of variable names as strings.
        %                          This is just for displaying labels on
        %                          graphs, can be left undefined
        %        'solver'          function handle to a MATLAB ode solver,
        %                          default @ode15s
        %
        %    All parameters can also be set after initialisation by the
        %    syntax (for example)
        %
        %        sim = simulate()
        %        sim.Tspan = [0 100];
            
            params = parse_inputs(varargin{:});
            names = fieldnames(params);
            for i = 1:numel(names)
                obj.(names{i}) = params.(names{i});
            end
        end
        
        % Main simulation method
        function simulate(obj)
            %SIMULATE    Perform a simulation
            %    Syntax:
            %        sim = ReactionDiffusion(...)
            %           ...
            %        sim.simulate()
            %
            %    As long as parameters are all defined, simulation will
            %    proceed with stored parameters.
            check_parameters(obj);
            fprintf('\tComputing discretisation ... ')
            [obj.x, obj.xa, obj.D, obj.m, obj.n] = discretise(obj);
            obj.JPattern = compute_jpattern(obj);
            fprintf('done\n')
            obj.odeopts = odeset('JPattern', obj.JPattern);
            tic
            fprintf('\tSolving ode system ... ')
            obj.odesol = obj.solver(@(t, y) obj.rhsfun(t, y), obj.Tspan, ...
                obj.y0vals, obj.odeopts);
            fprintf('done (elapsed time %.2f seconds)\n', toc)
        end
        
        % Utility routines
        function image(obj, idx)
            %IMAGE   Display an image of one variable of the simulation
            %
            %    sim.image(1) would display an image of the first variable
            %        with time on the x-axis, and space on the y-axis.
            %
            %    The time resolution is chosen to be in proportion with the
            %    spatial resolution. The ability to change this may be
            %    added in future
            %
            
            nT = round(obj.n / 3 * 4);
            T = linspace(obj.Tspan(1), obj.Tspan(end), nT);
            Y = deval(obj.odesol, T, idx:obj.m:(obj.m*obj.n));
            imagesc(T, obj.x, Y)
            colormap hot
            axis tight
            xlabel('t')
            if ~isempty(obj.varnames)
                title(obj.varnames{idx})
            else
                title(sprintf('y_%d', idx))
            end
            ylabel('x')
            set(gca, 'ydir', 'normal')
            
        end
        function animate(obj, idx)
            %ANIMATE    Perform an animation of a solution
            %    sim.animate([1 2])   Animate the solution, showing the
            %        first and second variables (in their own axes)
            %
            %    sim.animate()        Animate all variables
            %    sim.animate(1)       Animate just the first variable
            %
            %    Animation speed and frame rate can be controlled by
            %    setting the parameters 'animation_speedup' and
            %    'animation_framerate'. By default there is no speedup, and
            %    the framerate is 20. If your system can't display this
            %    smoothly, try reducing the framerate. Note that the
            %    framerate corresponds to real time not simulated time, so
            %    if the framerate is 20, it will attempt to display 20
            %    frames per second on the screen, irrespective of what the
            %    speedup is.
            %        e.g
            %    sim.animation_framerate = 10
            %    sim.animate()
            if nargin == 1
                idx = 1:obj.m;
            end
            wt_end = diff(obj.Tspan) / obj.animation_speedup;
            dt = 1 / obj.animation_framerate;
            nFrames = obj.animation_framerate * wt_end;
            T = linspace(obj.Tspan(1), obj.Tspan(end), nFrames);
            Y = deval(obj.odesol, T);
            ni = numel(idx);
            hy = zeros(ni, 1);
            ha = zeros(ni, 1);
            for i = 1:numel(idx)
                subplot(ni, 1, i)
                hy(i) = plot(obj.x, Y(idx(i):obj.m:end, 1));
                ha(i) = gca;
                yy = Y(idx(i):obj.m:end, :);
                set(ha(i), 'ylim', [min(yy(:)), max(yy(:))]);
                if (i == 1)
                    ht = title(sprintf('t = %.2f', T(1)));
                end
                if ~isempty(obj.varnames)
                    ylabel(obj.varnames{i});
                else
                    ylabel(sprintf('y_%d', idx(i)));
                end
            end
            xlabel('x');
            linkaxes(ha, 'x');
            pause(dt);
            for iFrame = 2:nFrames
                tic
                for i = 1:ni
                    set(hy(i), 'YData', Y(idx(i):obj.m:end, iFrame))
                end
                set(ht, 'String', sprintf('t = %.2f', T(iFrame)))
                drawnow
                pause(dt - toc);
            end
        end
        function [t, y] = export(obj, t, idx)
            %EXPORT    Export solution of selected indices at selected time
            %    values
            % 
            %    [T, Y] = sim.export()   Exports the entire solution with
            %                            200 evenly-spaced time points
            %                            spanning sim.Tspan. Y{i}
            %                            corresponds to the solution of the
            %                            i-th variable
            %
            %    [T, Y] = sim.export(Ts)  Exports entire solution at
            %                            selected times. 
            %
            %        If Ts is scalar, it's treated as a time step, and Y 
            %        corresponds to the solution at 
            %            sim.Tspan(1):Ts:sim.Tspan(2)
            %        If Ts has two elements, it's treated to be an initial
            %        and final time (both of which must be within
            %        sim.Tspan), and the solution is interpolated to a
            %        uniform grid of 200 time steps within this range
            %        If Ts is a vector of length greater than two, the
            %        solution is evaluated for each time in Ts.
            %
            %    [T, Y] = sim.export(Ts, idx)
            %
            %        idx is a vector (can be length one) of indices. Only
            %        the the solution elements of idx are exported.
            %    e.g.
            %    [T, Y] = sim.export(Ts, [2 3])
            %        Y{1} will correspond to the second variable
            %        Y{2} will correspond to the third variable
            %
            %        if idx is a single index, then Y is a matrix, not a
            %        cell array
            % Figure out t first
            if nargin == 1
                t = linspace(obj.Tspan(1), obj.Tspan(2), 200);
            elseif isscalar(t)
                t = obj.Tspan(1):t:obj.Tspan(2);
            elseif numel(t) == 2
                t = linspace(t(1), t(2), 200);
            end
            % And now idx
            if nargin < 3
                idx = 1:obj.m;
            end
            if isempty(obj.odesol)
                error('Unable to export: no simulation data')
            end
            
            if numel(idx) == 1
                y = deval(obj.odesol, t, idx:obj.m:(obj.m*obj.n));
            else
                y = cell(numel(idx), 1);
                for i = 1:numel(idx)
                    y{i} = deval(obj.odesol, t, idx(i):obj.m:(obj.m*obj.n));
                end
            end
            
            
        end
        
        % Specialised setters
        function set.y0(obj, val)
            if isnumeric(val)
                val = num2cell(val);
            end
            obj.y0 = val;
        end
        function set.diffusion(obj, val)
            if isnumeric(val)
                val = num2cell(val);
            end
            obj.diffusion = val;
        end
    end
    methods (Access = protected)
        function dy = rhsfun(obj, t, y)
            dy = obj.kinetics_fcn(t, obj.x, reshape(y, obj.m ,obj.n));
            dy = dy(:) + obj.D * y(:);
        end
        function y0 = y0vals(obj)
            y0 = zeros(obj.m, obj.n);
            for i = 1:obj.m
                if isnumeric(obj.y0{i})
                    y0(i,:) = obj.y0{i};
                else
                    y0(i,:) = obj.y0{i}(obj.x);
                end
            end
            y0 = y0(:);            
        end
        function J = compute_jpattern(obj)
            C = repmat({spones(ones(obj.m))}, obj.n, 1);
            J = spones(spones(obj.D) + blkdiag(C{:}));
        end
        function check_parameters(obj)
            %TODO: Add more than simple existence checks
            if isempty(obj.diffusion)
                error('Diffusion coefficients/functions undefined')
            end
            if isempty(obj.y0)
                error('Initial conditions undefined')
            end
            if isempty(obj.kinetics_fcn)
                error('Kinetics function undefined')
            end
            
        end
        function [x, xa, D, m, n] = discretise(obj)
            % First perform the spatial discretisation of the domain
            x = obj.xlim(1):obj.h:obj.xlim(2);
            xa = [obj.xlim(1) - obj.h/2, x + obj.h/2];
            m = numel(obj.diffusion);
            n = numel(x);
            
            % Now construct the diffusion operator block by block
            Dop = cell(m, 1);
            for i = 1:m
                if ~isnumeric(obj.diffusion{i})
                    a = obj.diffusion{i}(xa);
                else
                    a = obj.diffusion{i};
                end
                Dop{i} = diffop(n, obj.h, a, obj.boundary);
            end
            % Compute permutation matrix to get localised ordering
            I = 1:(m*n);
            J = reshape(I, n, m)';
            J = J(:);
            Q = sparse(I, J, 1., m*n, m*n);
            D = Q * blkdiag(Dop{:}) * Q.';
        end
        
    end
    
end

function params = parse_inputs(varargin)
p = inputParser();
p.addParamValue('h', 1e-2);
p.addParamValue('xlim', [0 1]);
p.addParamValue('diffusion', []);
p.addParamValue('y0', []);
p.addParamValue('Tspan', [0 10]);
p.addParamValue('boundary', 'noflux');
p.addParamValue('solver', @ode15s);
p.addParamValue('kinetics_fcn', [])
p.addParamValue('animation_speedup', 1);
p.addParamValue('animation_framerate', 20);
p.addParamValue('varnames', {});
p.parse(varargin{:});

params = p.Results;
end
function D = diffop(N, h, a, varargin)
%DIFFOP Create the square diffusion operator (including the h^2) for a
%single block
%
%    D = DIFFOP(N, h, a)  Create a diffusion operator for a grid of N
%        points with mesh size h. By default zero flux boundaries are used.
%        If a is a scalar, it is treated as a constant diffusion
%        coefficient
%        If a is a vector, it is assumed to be a spatially varying
%        diffusion coefficient evaluated on the mesh -h/2:h:N*h+h/2.
%
%    D = DIFFOP(N, h, a, boundary)
%        Boundary can be either 'noflux' (default) or 'periodic'
%        if the boundary conditions are periodic, a(1) and a(end) must match.

if (nargin == 3)
    bc = 'noflux';
else
    bc = varargin{1};
end

% Create the sub, main, super diagonals
if isscalar(a)
    a = repmat(a, N+1, 1);
else
    a = a(:);
end

A = 1/h^2 * [[a(2:end-1); 0], -(a(1:end-1) + a(2:end)), [0; a(2:end-1)]];
D = spdiags(A, [-1 0 1], N, N);

% Fix up the boundary conditions
switch bc
    case 'noflux'
        D(1, 2)       = 1/h^2 * (a(1) + a(2));
        D(end, end-1) = 1/h^2 * (a(end-1) + a(end));
    case 'periodic'
        if (a(1) ~= a(end))
            warning('reactiondiffusion:notperiodic', ...
                'diffusion coefficient function not periodic, using a(1) for a(end)')
        end
        D(1, end) = 1/h^2 * a(1);
        D(end, 1) = 1/h^2 * a(1);
    otherwise
        error('Unknown boundary type  use ''noflux'' or ''periodic''')
end
end
