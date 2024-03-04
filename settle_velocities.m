function [ns, force_scale]=settle_velocities(...
    ns_in, no_slip, free_slip, force_in, force_scale, tolerance, dt)
    % Expects that ns_in is already scaled by force_scale. The ns returned
    % will similarly be scaled by the force_scale returned.

    ns = ns_in;
    force.x = force_scale * force_in.x;
    force.y = force_scale * force_in.y;
    [~, ux, uy] = ns.get_macro();
    u_mag = sqrt(ux.^2 + uy.^2);

    max_du = double(ux(1:0));  % Maintain the type of the array...
    max_u = double(ux(1:0));
    net_u = double(ux(1:0));
    net_mom = double(ux(1:0));
    window = round(length(mean(uy)) / sqrt(ns_in.lattice.cs2));
    it = 0;
    all_it = 0;
    while 1
        if max(max(abs(u_mag))) > 1e-3
            % We are exiting stokes regime. Scale down the forces more.
            scale_by = 1e-1;
            force_scale = force_scale * scale_by;

            ns = NavierStokes(ns_in.nu, ns_in.ux*force_scale, ...
                ns_in.uy*force_scale, 1 + (ns_in.rho-1)*force_scale);
            force.x = force_scale * force_in.x;
            force.y = force_scale * force_in.y;
            ns.update_macro(no_slip, free_slip, force);
            [rho, ux, uy] = ns.get_macro();
            u_mag = sqrt(ux.^2 + uy.^2);
            it = 0;
            fprintf("It %i, decreasing force scale: now %g\n", all_it, force_scale);
        elseif it > window && max(max_u(all_it-window:all_it)) < 1e-6
        end

        % Stop when the average maximum change over the window is below the
        % specified tolerances:
        %   1. as a fraction of the maximum velocity
        %   2. as an absolute quantity (in case velocities are ~0)
        if it > window
            max_du_over_window = max(max_du(all_it-window:all_it));
%             max_net_u_over_window = max(net_u(all_it-window:all_it));
            max_net_mom_over_window = max(net_mom(all_it-window:all_it));
            max_u_over_window = max(max_u(all_it-window:all_it));
            % Temporarily...
            if (it/window > 10) & (max_u_over_window < 1e-8)
                break;
            end
            if (max_du_over_window < 10*eps || ...
                    max_du_over_window < tolerance*max_u_over_window) && ...
                (max_net_mom_over_window < 10*eps || ...
                    max_net_mom_over_window < tolerance*max_u_over_window)
                break;
            end

            % For now, if we've been trying long enough, forget about the
            % net velocities, they might be steady state...
%             if (max_du_over_window < 10*eps || ...
%                     max_du_over_window < tolerance*max_u_over_window) && ...
%                 (all_it > 3000)
%                 break;
%             end
        end

        it = it + 1;
        all_it = all_it+1;
        ns = ns.enforce_bc_free_slip(free_slip);
        ns = ns.update_macro(no_slip, free_slip, force);
        ns = ns.tick(no_slip, free_slip, force, dt);
        [rho, ux, uy] = ns.get_macro();

        last_u_mag = u_mag;
        u_mag = sqrt(ux.^2 + uy.^2);
        if all_it > length(max_du)
            grow_by = ceil(0.25 * length(max_du));
            max_du = [max_du, nan(1, grow_by)];
            max_u = [max_u, nan(1, grow_by)];
            net_u = [net_u, nan(1, grow_by)];
            net_mom = [net_mom, nan(1, grow_by)];

%             semilogy(max_du, 'k--', 'DisplayName', '\Delta u', 'linewidth', 1);
%             hold on; 
%             semilogy(net_u, 'r-', 'DisplayName', 'mean u', 'linewidth', 2);
%             semilogy(net_mom, 'g-', 'DisplayName', 'mean mom', 'linewidth', 2);
%             semilogy(max_u, 'k-', 'DisplayName', 'u', 'linewidth', 2);
%             yline(eps, 'DisplayName', 'double resolution');
%             legend('location', 'sw')
%             xlabel("Iteration");
%             title("Velocity settling results");
%             set(gca, 'fontsize', 12)
%             hold off; 
%             drawnow;
        end
        max_du(all_it) = max(max(abs(u_mag - last_u_mag)));
        max_u(all_it) = max(max(abs(double(u_mag))));
        net_u(all_it) = max(abs(mean(uy)));
        net_mom(all_it) = max(abs(mean(uy.*rho)));
    end

    fprintf("Settled after %i its or %.4g wavelengths\n", all_it, all_it/window);
end