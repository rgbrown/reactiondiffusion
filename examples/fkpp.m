%% FKPP Equation with variable diffusion
% Simple use of ReactionDiffusion to model a problem with spatially varying
% diffusion

%%
% We use standard FKPP reaction kinetics
f = @(t, x, y) y .* (1 - y);

%% Diffusion coefficient function
% We'll use a comb-like diffusion function
alpha = 15;
f1 = @(x) exp(-alpha*x.^2);
d = @(x) 1e-3 - 9.9e-4*f1(mod(5*x + 1, 2) - 1);
x = linspace(-1, 1, 1000);
plot(x, d(x));

%% Initial conditions
% We'll start from a smooth initial profile
y0 = @(x) 0.5 * (1 - tanh(200*(x +0.8)));
plot(x, y0(x))

%% Simulate and plot
sim = ReactionDiffusion('kinetics_fcn', f, 'diffusion', {d}, ...
    'method', 'fd', 'n', 2000, 'Tspan', [0 60])
sim.y0 = {y0};
sim.animation_speedup = 6;
sim.simulate()

sim.image(1)

%%
sim.animate()