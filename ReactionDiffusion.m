classdef ReactionDiffusion < handle
    properties
        % Setable properties. These define everything we need to know
        xlim         %
        n            % Spatial discretisation size
        diffusion    % diffusion coefficients. Some may be zero
        y0
        kinetics_fcn
        Tspan
        solver
        method
        varnames
        boundary
        odeopts
        image_nx
        image_nt
        animation_speedup
        animation_framerate
    end
    properties (SetAccess = private)
        % Workspace properties
        m
        fd
        spectral
        JPattern
        odesol
        Tsol
        Ysol
    end
    
    methods
        function obj = ReactionDiffusion(varargin)
            %ReactionDiffusion Construct reaction diffusion object.
            
            % This function is to set the parameters of the problem
            params = parse_inputs(varargin{:});
            names = fieldnames(params);
            for i = 1:numel(names)
                obj.(names{i}) = params.(names{i});
            end
            
        end
        function init_spectral(obj)
            % compute discretisation etc for spectral method
            N = obj.n - 1;
            [obj.spectral.D1, obj.spectral.x] = cheb(N);
            obj.spectral.x = obj.spectral.x(:).';
            obj.spectral.D2 = obj.spectral.D1^2;
            obj.spectral.BC = -obj.spectral.D1([1 N+1], [1 N+1]) \ ...
                obj.spectral.D1([1 N+1], 2:N);
            obj.m = numel(obj.diffusion);
            
            Dop = cell(obj.m, 1);
            is_diff = false(obj.m, 1);
            for i = 1:obj.m
                if ~isnumeric(obj.diffusion{i})
                    a = obj.diffusion{i}(obj.spectral.x);
                    Dop{i} = obj.spectral.D1 * diag(a) * obj.spectral.D1;
                    is_diff(i) = true;
                else
                    a = obj.diffusion{i};
                    Dop{i} = a * obj.spectral.D2;
                    if any(a)
                        is_diff(i) = true;
                    end
                end
            end
            obj.spectral.i_diff = find(is_diff);
            I = 1:(obj.m * obj.n);
            J = reshape(I, obj.n, obj.m)';
            J = J(:);
            Q = sparse(I, J, 1., obj.m*obj.n, ...
                obj.m*obj.n);
            obj.spectral.D = Q * blkdiag(Dop{:}) * Q.';
            
            obj.spectral.i_map = setdiff(1:obj.m*obj.n, ...
                [obj.spectral.i_diff;
                obj.m*(obj.n-1) + obj.spectral.i_diff]);
            obj.spectral.ntot = obj.m * obj.n;
            obj.spectral.y = zeros(obj.spectral.ntot, 1);
            
            
        end
        function init_fd(obj)
            % compute discretisation etc for fd method
            obj.fd.x  = linspace(obj.xlim(1), obj.xlim(2), obj.n);
            obj.fd.h  = (obj.xlim(2) - obj.xlim(1)) / (obj.n - 1);
            obj.fd.xa = [obj.xlim(1) - obj.fd.h/2, obj.fd.x + obj.fd.h/2];
            obj.m     = numel(obj.diffusion);
            obj.n     = numel(obj.fd.x);
            
            % Now construct the diffusion operator block by block
            Dop = cell(obj.m, 1);
            for i = 1:obj.m
                if ~isnumeric(obj.diffusion{i})
                    a = obj.diffusion{i}(obj.fd.xa); % if it is a function handle
                else
                    a = obj.diffusion{i};
                end
                Dop{i} = diffop(obj.n, obj.fd.h, a, obj.boundary);
            end
            % Compute permutation matrix to get localised ordering
            I = 1:(obj.m*obj.n);
            J = reshape(I, obj.n, obj.m)';
            J = J(:);
            Q = sparse(I, J, 1., obj.m*obj.n, obj.m*obj.n);
            obj.fd.D = Q * blkdiag(Dop{:}) * Q.';
        end
        function simulate(obj)
            %SIMULATE    Perform a simulation
            %    Syntax:
            %        sim = ReactionDiffusion(...)
            %           ...
            %        sim.simulate()
            %
            %    As long as parameters are all defined, simulation will
            %    proceed with stored parameters.
            obj.check_parameters();
            fprintf('\tComputing discretisation ... ')
            switch obj.method
                case 'fd'
                    obj.init_fd();
                    obj.JPattern = compute_jpattern(obj);
                    obj.odeopts = odeset(obj.odeopts, 'JPattern', obj.JPattern);
                    f = @(t, y) obj.rhs_fd(t, y);
                case 'spectral'
                    obj.init_spectral();
                    f = @(t, y) obj.rhs_spectral(t, y);
            end
            y0v = obj.y0vals();
            fprintf('done\n')
            tic
            fprintf('\tSolving ode system ... ')
            if numel(obj.Tspan) == 2 
                obj.odesol = obj.solver(f, obj.Tspan, y0v, obj.odeopts);
                obj.Tsol = obj.odesol.x;
                obj.Ysol = obj.odesol.y.';
            else
                [obj.Tsol, obj.Ysol] = ...
                    obj.solver(f, obj.Tspan, y0v, obj.odeopts);
                obj.odesol = [];
            end
            fprintf('done (elapsed time %.2f seconds)\n', toc)
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
        function y = unpack_spectral(obj, yreduced)
            ncols = size(yreduced, 2);
            y = zeros(obj.spectral.ntot, ncols);
            y(obj.spectral.i_map, :) = yreduced;
            for i = 1:numel(obj.spectral.i_diff)
                ii = obj.spectral.i_diff(i);
                
                y([ii; ii + obj.m*(obj.n-1)], :) = ...
                    obj.spectral.BC * y(obj.m*(1:(obj.n-2)) + ii, :);
            end
        end
        function Y = soln(obj, idx, t, x)
            %SOLN Extract solution at specified time and spatial locations
            if isempty(t)
                t = obj.Tsol;
            end
            if isempty(x)
                switch obj.method
                    case 'fd'
                        x = obj.fd.x;
                    case 'spectral'
                        x = obj.spectral.x;
                end
            end
            % Bounds checking
            if any(t < obj.Tspan(1)) || any(t > obj.Tspan(end))
                error('Rdsolve:t_out_of_range', ...
                    'specified t values outside of simulation range');
            end
            if any(x < obj.xlim(1)) || any(x > obj.xlim(end))
                error('Rdsolve:x_out_of_range', ...
                    'specified x values outside of simulation range');
            end
            
            % Get solutions at node points at each time value
            if ~isempty(obj.odesol)
                Y = deval(obj.odesol, t);
            else
                % Check to see if the times are actually simulation values
                [tf, loc] = ismember(t, obj.Tspan);
                % If not, interpolate
                if all(tf)
                    Y = obj.Ysol(loc, :).';
                else
                    Y = obj.Ysol.';
                    Y = spline(obj.Tsol, Y, t);
                end
            end
            
            % Fill out boundary conditions if using spectral method
            if isequal(obj.method, 'spectral')
                Y = obj.unpack_spectral(Y);
            end
            % Extract the variable out that we need
            Y = Y(idx:obj.m:end, :);
            
            % Now spatial interpolation
            switch obj.method
                case 'fd'
                    xs = obj.fd.x;
                case 'spectral'
                    xs = obj.spectral.x;
            end
            [tf, loc] = ismember(x, xs);
            if all(tf)
                Y = Y(loc, :);
            else
                Y = spline(xs, Y.', x).';
            end       
            
        end
        
        function varargout = image(obj, idx)
            %IMAGE   Display an image of one variable of the simulation
            %
            %    sim.image(1) would display an image of the first variable
            %        with time on the x-axis, and space on the y-axis.
            %
            %    The time resolution is chosen to be in proportion with the
            %    spatial resolution. The ability to change this may be
            %    added in future
            %
            
            % Work out the times
            if numel(obj.Tspan) == 2
                nt = obj.image_nt;
                T = linspace(obj.Tspan(1), obj.Tspan(2), nt);
            else
                nt = numel(obj.Tspan);
                warning('Rdsolve:ignore_nt', ...
                    'Ignoring specified image_nt value, using number of elements in Tspan instead');
                T = obj.Tspan;
            end
            % Fill in missing entries if we're using spectral method
            if isequal(obj.method, 'spectral')
                if isempty(obj.image_nx)
                    nx = round(3/4 * nt);
                else
                    nx = obj.image_nx;
                end
                x = linspace(obj.xlim(1), obj.xlim(end), nx);
            else
                if isempty(obj.image_nx)
                    x = obj.fd.x;
                else
                    x = linspace(obj.xlim(1), obj.xlim(end), obj.image_nx);
                end
            end
            Y = obj.soln(idx, T, x);
                        
            hI = imagesc(T, x, Y);
            colormap(hot(16384));
            axis tight
            xlabel('t')
            if ~isempty(obj.varnames)
                title(obj.varnames{idx})
            else
                title(sprintf('y_%d', idx))
            end
            ylabel('x')
            set(gca, 'ydir', 'normal')
            switch nargout
                case 0
                    varargout = {};
                case 1
                    varargout = {hI};
                case 2
                    varargout = {hI, Y};
            end
                    
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
            wt_end = diff(obj.Tspan([1 end])) / obj.animation_speedup;
            dt = 1 / obj.animation_framerate;
            nFrames = obj.animation_framerate * wt_end;
            if isempty(obj.image_nx)
                switch obj.method
                    case 'fd'
                        x = obj.fd.x;
                    case 'spectral'
                        x = obj.spectral.x;
                end
            else
                x = linspace(obj.xlim(1), obj.xlim(end), obj.image_nx);
            end
            T = linspace(obj.Tspan(1), obj.Tspan(end), nFrames);
            
            ni = numel(idx);
            Y = cell(ni, 1);
            hy = zeros(ni, 1);
            ha = zeros(ni, 1);
            for i = 1:numel(idx)
                Y{i} = obj.soln(idx(i), T, x);
                subplot(ni, 1, i)
                hy(i) = plot(x, Y{i}(:, i));
                ha(i) = gca;
                set(ha(i), 'ylim', [min(Y{i}(:)), max(Y{i}(:))]);
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
                    set(hy(i), 'YData', Y{i}(:, iFrame))
                end
                set(ht, 'String', sprintf('t = %.2f', T(iFrame)))
                drawnow
                pause(dt - toc);
            end
        end
        
    end
    
    
    methods (Access = protected)
        function dy = rhs_spectral(obj, t, yreduced)
            y = obj.unpack_spectral(yreduced);
            dy = obj.kinetics_fcn(t, obj.spectral.x, ...
                reshape(y, obj.m, obj.n));
            dy = dy(:) + obj.spectral.D * y(:);
            dy = dy(obj.spectral.i_map);
        end
        function dy = rhs_fd(obj, t, y)
            dy = obj.kinetics_fcn(t, obj.fd.x, reshape(y, obj.m, obj.n));
            dy = dy(:) + obj.fd.D * y(:);
        end
        function y0 = y0vals(obj)
            y0 = zeros(obj.m, obj.n);
            for i = 1:obj.m
                if isnumeric(obj.y0{i})
                    y0(i,:) = obj.y0{i};
                else
                    switch obj.method
                        case 'fd'
                            y0(i,:) = obj.y0{i}(obj.fd.x);
                        case 'spectral'
                            y0(i,:) = obj.y0{i}(obj.spectral.x);
                    end
                end
            end
            y0 = y0(:);
            if isequal(obj.method, 'spectral')
                y0 = y0(obj.spectral.i_map);
            end
        end
        
        function J = compute_jpattern(obj)
            C = repmat({spones(ones(obj.m))}, obj.n, 1);
            J = spones(spones(obj.fd.D) + blkdiag(C{:}));
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
    end
    
    
    
end

function params = parse_inputs(varargin)
p = inputParser();
p.addParamValue('n', 100);
p.addParamValue('xlim', [-1 1]);
p.addParamValue('diffusion', []);
p.addParamValue('y0', []);
p.addParamValue('Tspan', [0 10]);
p.addParamValue('boundary', 'noflux');
p.addParamValue('method', 'fd'); % fd or spectral
p.addParamValue('solver', @ode15s);
p.addParamValue('odeopts', odeset());
p.addParamValue('kinetics_fcn', [])
p.addParamValue('varnames', {});
p.addParamValue('image_nx', []);
p.addParamValue('image_nt', 2000);
p.addParamValue('animation_framerate', 20);
p.addParamValue('animation_speedup', 5);
p.parse(varargin{:});

params = p.Results;
end
function D = diffop(N, h, a, varargin)
%DIFFOP Create the square finite difference diffusion operator
% (including the h^2) for a single block
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
function [D,x] = cheb(N)
if N == 0
    D = 0;
    x = 1;
    return
end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+eye(N+1));
D = D - diag(sum(D,2));
end