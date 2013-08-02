%% Reaction Diffusion Example: The Brusselator
% We want to simulate the Brusselator as an example of how to use the 
% somewhat general purpose reaction-diffusion simulator.

%% Initialise the problem
% First specify the problem domain. Just the boundaries, not the internal
% mesh points.
xlim = [-1 1];

%%
% Now specify the parameters. We'll make b our spatially varying parameter.
% The best way to do this is to make it a function of x using anonymous
% functions
a = 1;
b = @(x) 1.5 * x + 2.5;

%%
% Now we'll define the kinetics. The Brusselator is simple enough that you
% can just use an anonymous function. For a more complex problem you'd
% write an m file whose signature is of the form
%
%      |function f = mykinetics(t, x, y)|
%
% The kinetics can be space and time-varying
f = @(~, x, y) ...
    [a + y(1, :).^2 .* y(2, :) - (b(x) + 1) .* y(1, :);
    b(x) .* y(1, :) - y(1, :).^2 .* y(2, :)];

%%
% A couple of extra things
varnames = {'\phi', '\psi'}; % only need these to make graphs look nice
%h = 5e-2;                    % spatial step
n = 200;
diffusion = {@(x) 1e-2 * exp(-x.^2), 0};  % initially we'll do phi-diffusion

%%
% define some initial conditions as functions of x. They can also be scalars
y0 = {1, @(x) b(x) + 0.2 * sin(4*pi*x)}; 
%% 
% Now we'll create our reaction-diffusion system
sim = Rdsolve('kinetics_fcn', f, 'diffusion', diffusion, 'n', n, ...
    'xlim', xlim, 'y0', y0, 'varnames', varnames, 'Tspan', [0 10], ...
    'method', 'spectral');

%% Time to simulate the system
sim.simulate()

%% Plot results
figure(1)
sim.image(1); % Display the first variable in image form

%%
% Hmmm, looks like we didn't simulate for long enough, and could use better
% resolution. Let's try again
sim.Tspan = [0 40];
sim.h = 1e-3;
sim.simulate();
figure(2)
sim.image(2);


%% Animation time
% Sometimes it's better to look at a movie than at a picture. (you need to
% run the code not look at the static HTML page to see this)
figure(3)
sim.animation_speedup = 2; % speed up simulation time by a factor 2 in real time
sim.animate(); % Animate both variables together. Could also do e.g.
% sim.animate(1); % to animate the first one

%% Funky spikes
% Let's put diffusion on \psi instead to see the funky spikes. Best viewed
% at lower resolution 
figure(4)
sim.h = 1e-2;
sim.diffusion{1} = 0;
sim.diffusion{2} = 1e-4;
sim.simulate();
sim.image(1);
%%
% Yup, there they go


%% Now let's do an exercise to show visually that our method converges with small enough mesh
sim.Tspan = [0 10];
sim.diffusion = {1e-2, 0};
sim.h = 0.1;
sim.simulate();
figure(5)
subplot(2,2,1)
sim.image(1);
title(sprintf('h = %.2g', sim.h));

sim.h = 0.05;
sim.simulate()
subplot(2,2,2)
sim.image(1);
title(sprintf('h = %.2g', sim.h));


sim.h = 0.01;
sim.simulate()
subplot(2,2,3)
sim.image(1);
title(sprintf('h = %.2g', sim.h));


sim.h = 0.005;
sim.simulate()
subplot(2,2,4)
sim.image(1);
title(sprintf('h = %.2g', sim.h));

%% Now we'll try comparing different solvers
sim.h = 5e-3;
sim.Tspan = [0 10];
fprintf('Testing the stiff solvers:\n')
fprintf('ode15s:\n')
sim.solver = @ode15s;
sim.simulate();
fprintf('ode23s:\n')
sim.solver = @ode23s;
sim.simulate();
fprintf('\nTesting the explicit solvers:\n')
fprintf('ode23:\n')
sim.solver = @ode23;
sim.simulate();
fprintf('ode45:\n')
sim.solver = @ode45;
sim.simulate();






