xlim = [-1 1];
%v1fun = @(x) 5e-3 * (1 + tanh((x - 0.5)/0.125)) - 29.6e-3;
v1fun = @(x) -22.5e-3;

%y0fun = @(x) exp(-(x-0.5).^2 / 100);
[fn, varnames, xlim, y0_loc] = ermentrout('xlim', xlim, 'v1fun', v1fun); 

% Put diffusion on the membrane potential equation
diffusion = [1e-4, 0, 0];
[fode] = rdcollocation(fn, [1e-3, 0, 0], y0_loc);

%[T, Y] = ode15s(fode, [0 20], 