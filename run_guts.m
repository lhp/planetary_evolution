%% Not intended to run on its own. Should be called from run.m or run_from_checkpoint.m
% Required:
%   run_timer           Timer tracking wall time since beginning of job
%   starting_tick       Current iteration
%   total_ticks         Total goal iterations
%   model               Model to run
%   scaling             Relationship between model and reality       
%   viz_every           How often (in wall time seconds) to visualize results
%   num_checkpoints     How many checkpoints to make
%   output_folder       Where to write checkpoints and final state
%   max_clock_hours     How long to run (WT) before checkpointing and exiting

run_total_time = toc(run_timer);
run_timer = tic;
prediction_timer = tic;

% Set up checkpointing and visualization.
if ~isempty(output_folder)
    mkdir(output_folder);
    checkpoint_folder = fullfile(output_folder, "checkpoints");
    mkdir(checkpoint_folder);
elseif num_checkpoints > 0
    fprintf("Checkpointing is requested, but no output folder is provided!\n");
    return
end
checkpoint_interval = max(floor(total_ticks / num_checkpoints), 1);
if viz_every < Inf
    viz = Visualizer(scaling);
end
last_viz_time = -Inf;

for i=starting_tick:total_ticks
    model = model.tick();

    tick_elapsed_time = toc(run_timer);
    run_total_time = run_total_time + tick_elapsed_time;
    if run_total_time - last_viz_time > viz_every
        viz = viz.update_display(model.macro, model.history);
        last_viz_time = run_total_time;
    end
    if mod(i, checkpoint_interval) == 1
        starting_tick = i + 1;
        checkpoint_file = fullfile(checkpoint_folder, ...
            sprintf("checkpoint_%s.mat", int2str(i)));
        save(checkpoint_file, "starting_tick", "total_ticks", "model", "scaling", ...
            "viz_every", "num_checkpoints", "output_folder", "max_clock_hours");
        prediction_elapsed_time = toc(prediction_timer);
        prediction_timer = tic;
        fprintf("Estimated remaining time : %g h.\n", ...
            prediction_elapsed_time/(60*60) * (total_ticks-i)/checkpoint_interval);
    end

    if run_total_time/(60*60) > max_clock_hours
        break;
    end
    run_timer = tic;
end

final_exit_code = 0;
% Checkpoint and visualize final state no matter what.
if viz_every < Inf
    viz = viz.update_display(model.macro, model.history);
    last_viz_time = run_total_time;
end
if ~isempty(output_folder)
    starting_tick = i + 1;
    if i == total_ticks
        checkpoint_file = fullfile(checkpoint_folder, ...
            sprintf("checkpoint_%s_final.mat", int2str(i)));
        final_exit_code = 0;
    else
        checkpoint_file = fullfile(checkpoint_folder, ...
            sprintf("checkpoint_%s.mat", int2str(i)));
        final_exit_code = 3;  % Incomplete.
    end

    save(checkpoint_file, "starting_tick", "total_ticks", "model", "scaling", ...
        "viz_every", "num_checkpoints", "output_folder");

    fprintf("Final state at %s saved to %s\n", int2str(i), checkpoint_file);
end

exit_code = final_exit_code;