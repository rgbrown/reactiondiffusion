function [fk, varnames, xlim, y0] = ermentrout(varargin)
%ERMENTROUT: Set up Ermentrout reaction diffusion kinetics function
p = inputParser();
p.addParamValue('v1fun', @(x) repmat(-22.5e-3, size(x))); 
p.addParamValue('xlim', [0 1]);
p.parse(varargin{:});
fv1 = p.Results.v1fun;
% Parameters
v2   = 25.0e-3;  %V
v4   = 14.5e-3;  %V
v5   = 8.0e-3;   %V
v6   = -15.0e-3; %V
Ca3  = 400e-9;   %M
Ca4  = 150e-9;   %M
phi_n = 2.664;    %s^-1
vL = -70e-3;     %V
vK = -90e-3;     %V
vCa = 80.0e-3;   %V

C = 19.635e-12;  %F
gL = 78.54e-12;  %S
gK = 314.16e-12; %S
gCa = 157e-12;   %S
Kd = 1.0e-6;     %M
BT = 100e-6;     %M
alpha = 7.9976e6; %M/C
% beta = 0.55;
kCa = 135.67537; %s^-1


ca0 = 400e-9;
v0 = -60e-3;
n0 = 0.3;
y0 = [v0; ca0; n0];
fk = @rhs;
varnames = {'v', 'ca', 'n'};
xlim = p.Results.xlim;

    function dy = rhs(t, x, y)
        v  = y(1, :); 
        ca = y(2, :); 
        n  = y(3, :);
        
        dy = zeros(size(y));
        
        v1 = fv1(x);
        
        v3 = -0.5*v5*tanh((ca - Ca3)/Ca4) + v6;
        m_infty = 0.5*(1 + tanh((v - v1)/v2));
        n_infty = 0.5*(1 + tanh((v - v3)/v4));
        rho = (Kd + ca).^2 ./ ((Kd + ca).^2 + Kd*BT);
        lambda_n = phi_n*cosh((v - v3)/(2*v4));
        
        dy(1, :)  = (-gL*(v - vL) - gK*n.*(v - vK) - gCa*m_infty.*(v - vCa))/C;
        dy(2, :) = rho.*(-alpha*gCa*m_infty.*(v - vCa) - kCa*ca);
        dy(3, :)  = lambda_n.*(n_infty - n);
    end

end


