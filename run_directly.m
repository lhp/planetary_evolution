%% Must be set in the following lines:
%   run_timer
%
%   d
%   H
%   normalized_core_C
%   Ra
%   L_l
%   kappa_l
%   solid
%   
%   viz_every
%   num_checkpoints
%   max_clock_hours
%   output_folder
run_timer = tic;

% Physical configuration of the system to be modeled.
% Example configuration: a 600 km thick mantle with all HPE in bottom 1/8.
d = [75e3 525e3];
H = [6e-11*600/75 0];
normalized_core_C = 2;
Ra = 1e5;
L_l = 160;
kappa_l = 1/6;
solid = [1 0];  % Bottom layer is solid.
total_time = 1 * 1e9*365*24*60*60;  % In this example, run for 1 Gy.

% Numerical parameters and configuration of visualization/data output.
velocity_skip = 20;  % Ratio between timesteps between momentum and heat solvers.
settling_tolerance = 1e-4;  % How close to steady state velocities need to be.
viz_every = 10;  % Visualize every 10 seconds
max_clock_hours = Inf;
num_checkpoints = 10;
output_folder = "test";

mkdir(output_folder);
save(fullfile(output_folder, "meta.mat"), ...
    "d", "H", "normalized_core_C", "Ra", "L_l", "kappa_l", "solid", ...
    "velocity_skip", "settling_tolerance", ...
    "num_checkpoints", "viz_every", "max_clock_hours");

[model, scaling, total_ticks]=setup_model(...
    d, H, normalized_core_C, Ra, L_l, kappa_l, solid, total_time, ...
    velocity_skip, settling_tolerance);

starting_tick = 1;
run_guts;