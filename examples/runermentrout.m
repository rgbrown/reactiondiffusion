%% Ermentrout Reaction Diffusion simulation
% Initialise the Ermentrout kinetics
xlim = [-1 1];
v1fun = @(x) 5e-3 * (1 + tanh(x/0.25)) - 29.6e-3;
%v1fun = @(x) -22.5e-3;

%y0fun = @(x) exp(-(x-0.5).^2 / 100);
[fn, varnames, xlim, y0] = ermentrout('xlim', xlim, 'v1fun', v1fun); 

% Put diffusion on the membrane potential equation
diffusion = [1e-4, 0, 0];

% Initialise the Reaction Diffusion simulation
sim = Rdsolve('kinetics_fcn', fn, 'xlim', xlim', ...
    'diffusion', diffusion, 'varnames', varnames, ...
    'method', 'spectral', 'n', 100, 'y0', y0);

sim.Tspan = [0 10];
sim.simulate()
sim.image(2)
%%
% Plot the time variation of membrane potential at the middle of the domain
[T, Y] = sim.export();
[~, idx] = min(abs(sim.x - 0.5));
plot(T, Y{1}(idx, :))

%%
% Now try decreasing the diffusion
sim.diffusion{1} = 1e-3;
sim.simulate()
sim.image(1)
%%
% Now plot the membrane potential at the middle of the grid over time
[T, Y] = sim.export();
[~, idx] = min(abs(sim.x - 0.5));
plot(T, Y{1}(idx, :))
%%
% And let's see an animation, first fast
sim.animation_speedup = 5;
sim.animate()
%%
% And now a slow version
sim.animation_speedup = 0.5;
sim.animate()