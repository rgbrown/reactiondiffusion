classdef Rdsolve < handle
    properties
        % Setable properties. These define everything we need to know
        xlim         % 
        n            % Number of spatial points
        diffusion    % diffusion coefficients. Some may be zero
        y0
        kinetics_fcn  
        Tspan
        solver
        method
        varnames
        boundary
    end
    properties (SetAccess = private)
        % Workspace properties
        fd
        spectral
        odeopts
        JPattern
        odesol
    end
    
    methods
        function obj = Rdsolve(varargin)
            %RDSOLVE Construct reaction diffusion object.
            
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
            obj.spectral.D2 = obj.spectral.D1^2;
            obj.spectral.BC = -obj.spectral.D1([1 N+1], [1 N+1]) \ ...
                obj.spectral.D1([1 N+1], 2:N);
            obj.spectral.
            
        end
        function init_fd(obj)
            % compute discretisation etc for fd method
            obj.fd.x  = linspace(obj.xlim(1), obj.xlim(2), obj.n);
            obj.fd.h  = (obj.xlim(2) - obj.xlim(1)) / (obj.n - 1);
            obj.fd.xa = [obj.xlim(1) - obj.fd.h/2, obj.fd.x + obj.fd.h/2];
            obj.fd.m     = numel(obj.diffusion);
            obj.fd.n     = numel(obj.fd.x);
            
            % Now construct the diffusion operator block by block
            Dop = cell(obj.fd.m, 1);
            for i = 1:obj.fd.m
                if ~isnumeric(obj.diffusion{i})
                    a = obj.diffusion{i}(obj.fd.xa); % if it is a function handle
                else
                    a = obj.diffusion{i};
                end
                Dop{i} = diffop(obj.fd.n, obj.fd.h, a, obj.boundary);
            end
            % Compute permutation matrix to get localised ordering
            I = 1:(obj.fd.m*obj.fd.n);
            J = reshape(I, obj.fd.n, obj.fd.m)';
            J = J(:);
            Q = sparse(I, J, 1., obj.fd.m*obj.fd.n, obj.fd.m*obj.fd.n);
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
                    obj.odeopts = odeset('JPattern', obj.JPattern);
                    y0v = obj.y0vals_fd();
                    f = @(t, y) obj.rhs_fd(t, y);
                case 'spectral'
                    obj.init_spectral();
                    obj.odeopts = odeset();
                    y0v = obj.y0vals_spectral();
                    f = @(t, y) obj.rhs_spectral(t, y);
            end
            fprintf('done\n')
            tic
            fprintf('\tSolving ode system ... ')
            obj.odesol = obj.solver(f, obj.Tspan, y0v, obj.odeopts);
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
    end
    methods (Access = protected)
        function dy = rhs_spectral(obj, t, y)
            
        end
        function dy = rhs_fd(obj, t, y)
            dy = obj.kinetics_fcn(t, obj.fd.x, reshape(y, obj.fd.m, obj.fd.n));
            dy = dy(:) + obj.fd.D * y(:);            
        end
        function y0 = y0vals_fd(obj)
            y0 = zeros(obj.fd.m, obj.fd.n);
            for i = 1:obj.fd.m
                if isnumeric(obj.y0{i})
                    y0(i,:) = obj.y0{i};
                else
                    y0(i,:) = obj.y0{i}(obj.fd.x);
                end
            end
            y0 = y0(:);            
        end
        function y0 = y0vals_spectral(obj)
        end
        function J = compute_jpattern(obj)
            C = repmat({spones(ones(obj.fd.m))}, obj.fd.n, 1);
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
p.addParamValue('kinetics_fcn', [])
p.addParamValue('varnames', {});
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