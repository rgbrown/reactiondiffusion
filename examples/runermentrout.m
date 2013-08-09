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
sim = ReactionDiffusion('kinetics_fcn', fn, ...
    'xlim', xlim', ...
    'diffusion', diffusion, ...
    'varnames', varnames, ...
    'method', 'spectral', ...
    'n', 200, ...
    'y0', y0);

sim.Tspan = linspace(0, 40, 1000);
sim.simulate()
sim.image(1)
%%
% Plot the time variation of membrane potential at the middle of the domain
t = linspace(sim.Tspan(1), sim.Tspan(2), 1000);
y = sim.soln(1, t, 0);
plot(t, y)

%%
% Now try decreasing the diffusion
sim.diffusion{1} = 1e-3;
sim.simulate()
sim.image(1)
%%
% Now plot the membrane potential at the middle of the grid over time
y = sim.soln(1, t, 0);
plot(t, y);
%%
% And let's see an animation, first fast
sim.image_nx = 1000;
sim.animation_speedup = 5;
sim.animate()
%%
% And now a slow version
sim.animation_speedup = 0.5;
sim.animate()