%% Ermentrout Reaction Diffusion simulation
% This script is to demonstrate the use of the ReactionDiffusion solver

%% Example 1: Ermentrout model
% This example uses the additional file |ermentrout.m|

%%
% First, we need to initialise the kinetics. We define a spatially varying
% v1 parameter for the model.
xlim = [-1 1];
v1fun = @(x) 5e-3 * (1 + tanh(x/0.25)) - 29.6e-3;
%%
% The function ermentrout provides us with a function handle to simulate
% the kinetics with, together with a few extra things, like the names of
% the variables in case you forget which is which
[fn, varnames, xlim, y0] = ermentrout('xlim', xlim, 'v1fun', v1fun); 

%%
% We want to run a simulation where we put a constant diffusion coefficient
% on the membrane potential equation, and no diffusion on the other
% variables
diffusion = [1e-4, 0, 0];

%%
% We now initialise the Reaction Diffusion simulation
sim = ReactionDiffusion('kinetics_fcn', fn, ...
    'xlim', xlim', ...
    'diffusion', diffusion, ...
    'varnames', varnames, ...
    'method', 'fd', ...
    'n', 2000, ...
    'y0', y0);

%%
% Let's do our first simulation. We provide a set of time values to
% interpolate the solution onto, simulate, and then plot the membrane
% potential
sim.Tspan = linspace(0, 40, 1000);
sim.simulate()
sim.image(1)
%%
% Plot the time variation of membrane potential at the middle of the domain
t = linspace(0, 40, 500);
y = sim.soln(1, t, 0);
plot(t, y)

%%
% Now try increasing the diffusion coefficient
sim.diffusion{1} = 1e-3;
sim.simulate()
sim.image(1)
%%
% Now plot the membrane potential at the middle of the grid over time
y = sim.soln(1, t, 0);
plot(t, y);
%%
% And let's see an animation
sim.image_nx = 1000;
sim.animation_speedup = 2;
sim.animate()