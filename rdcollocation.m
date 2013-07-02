function f = rdcollocation(kinfun, diffusion, y0_loc, varargin)
%RDCOLLOCATION  Solve reaction diffusion problem by collocation
%
%    Can solve multi-variable reaction diffusion equations by collocation
%    where some or all of the variables diffuse
p = parseinputs(varargin{:});
N = p.N;

% Sort out the Chebyshev differentiation matrix and boundary condition
% stuff
[D1,x] = cheb(N);
D2 = D1^2;
BC = -D1([1 N+1], [1 N+1]) \ D1([1 N+1], 2:N);

i_diff = find(diffusion ~= 0);
n_vars = numel(diffusion);
D_op   = diff_matrix(D2, i_diff, n_vars);
i_map  = create_mapping(N, i_diff, n_vars);
y_full = zeros(n_vars*(N+1), 1);

%y0 = repmat
f = @rhs;

% Construct mapping from reduced to full set of variables

    function dy = rhs(t, y_reduced)
        % Each diffusing quantity has two fewer mesh points which we need
        % to fill in
        y_full(i_map) = y_reduced;
        for i = i_diff
            y_full([1 N+1] + (i-1)*(N+1)) = BC * y_full((2:N) + (i-1)*(N+1));
        end
        dy = kinfun(t, x, reshape(y_full, n_vars, []));
        for i = i_diff
            idx = (1:N+1) + (i-1)*(N+1);
            dy(idx) = dy(idx) + diffusion(i) * D2 * y_full(idx);
        end
        dy = dy(i_map);        
    end

end

function params = parseinputs(varargin)
p = inputParser();
p.addParamValue('xlim', [-1 1]);
p.addParamValue('N', 20); % Number of Chebyshev points
p.parse(varargin{:});
params = p.Results;
end

function D = diff_matrix(D2, i_diff, n_vars)
N = size(D2, 1) - 1;
Dop = cell(n_vars, 1);
Z = sparse(N+1, N+1);
n_diff = numel(i_diff);
Dop(i_diff) = repmat({sparse(D2)}, n_diff, 1);
Dop(setdiff(1:n_vars, i_diff)) = repmat({Z}, n_vars - n_diff, 1); 
D = blkdiag(Dop{:});
end

function i_map = create_mapping(N, I, n)
i_full = 1:((N+1)*n);
for i = I
    i_full([1 N+1] + (i-1)*(N+1)) = nan;    
end
i_map = i_full(~isnan(i_full));

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