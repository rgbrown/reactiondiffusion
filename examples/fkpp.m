%% FKPP Equation with variable diffusion
f = @(t, x, y) y .* (1 - y);
alpha = 100;
y0 = @(x) 0.5 * (1 - tanh(alpha*(x - 0.2)));
d  = @(x) 0.001*(1 + 5*exp(-500*(x - 0.5).^2));
x = linspace(0, 1);

plot(x, d(x))
%%
sim = ReactionDiffusion('kinetics_fcn', f, 'diffusion', {d}, 'method', 'spectral', 'n', 1000)
sim.y0 = {y0}
sim.animation_speedup = 2;
sim.simulate()
sim.animate()
sim.diffusion{1} = 1e-3;
sim.simulate
sim.animate() 
