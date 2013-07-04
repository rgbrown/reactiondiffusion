function spectralrd(kinetics_fun)
%% Reaction diffusion by spectral collocation -- simple version

% We want to return a right hand side function, and a function for
% returning the data structure

% Construct differentiation matrix and spatial grid
[D, x] = cheb(N);

% For now we'll look at constant diffusion coefficients
D2 = D^2;

% For Neumann boundary conditions, the end points are determined by the
% rest of the values on the grid -- i.e. u1 and un are determined as
%  [u(1); u(N+1)] = B * u(2:N)
B = - D([1 N+1], [1 N+1]) \ D([1 N+1], 2:N);

% Extra things we need
n_vars = 3; % Number of variables in problem
i_diffusive = [1 2];
diffusion_coeff = [1e-3, 2e-4];
idx = reshape(1:(n_vars*(N+1)), n_vars, []); 

% Initialisation phase: work out which variables need to diffuse

    function dy = rhs(t, y)
        % Expand y into its full version
        
        
        % Evaluate kinetics
        dy = kinetics_fun(t, x, reshape(y, n_vars, []));
        dy = dy(:);
            
        % Do the diffusion thing
        for i = i_diffusive
            ii = idx(i, :);
            dy = dy;
            
            
            
        end
    end

    function inflate()
        
    end

    function deflate()
        
    end
end