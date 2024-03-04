classdef NavierStokes
    % Storing all variables flat substantially improves memory access
    % patterns, improving performance of the model.
    properties
        lattice
        
        f_in_1
        f_in_2
        f_in_3
        f_in_4
        f_in_5
        f_in_6
        f_in_7
        f_in_8
        f_in_9
        
        f_out_1
        f_out_2
        f_out_3
        f_out_4
        f_out_5
        f_out_6
        f_out_7
        f_out_8
        f_out_9
        
        f_eq_1
        f_eq_2
        f_eq_3
        f_eq_4
        f_eq_5
        f_eq_6
        f_eq_7
        f_eq_8
        f_eq_9
        
        nu
        vx
        vy

        rho
        ux
        uy
    end

    methods
        function obj=NavierStokes(nu, ux_0, uy_0, rho_0)
            d2q9 = get_lattice_d2q9();
            obj.lattice = d2q9;
            obj.nu = nu;
            
            for i=1:length(ux_0(:, 1))
                for j=1:length(ux_0(1, :))
                    obj.vx(i, j, :) = d2q9.vx;
                    obj.vy(i, j, :) = d2q9.vy;
                end
            end

            obj = obj.compute_feq(ux_0, uy_0, rho_0);
            for vi=1:9
                obj.(sprintf("f_in_%i", vi)) = obj.(sprintf("f_eq_%i", vi));
                obj.(sprintf("f_out_%i", vi)) = 0 * obj.(sprintf("f_eq_%i", vi));
            end

            obj.rho = rho_0;
            obj.ux = ux_0;
            obj.uy = uy_0;
        end

        function obj=enforce_bc_free_slip(obj, free_slip)  % [bottom top left right]
            bottom = free_slip(1);
            top = free_slip(2);
            left = free_slip(3);
            right = free_slip(4);
            if bottom
                obj.f_in_3(:, 1) = obj.f_in_5(:, 1);
                obj.f_in_6(:, 1) = obj.f_in_9(:, 1);
                obj.f_in_7(:, 1) = obj.f_in_8(:, 1);
            end
            if top
                obj.f_in_5(:, end) = obj.f_in_3(:, end);
                obj.f_in_9(:, end) = obj.f_in_6(:, end);
                obj.f_in_8(:, end) = obj.f_in_7(:, end);
            end
            if left
                obj.f_in_2(1, :) = obj.f_in_4(1, :);
                obj.f_in_6(1, :) = obj.f_in_7(1, :);
                obj.f_in_9(1, :) = obj.f_in_8(1, :);
            end
            if right
                obj.f_in_4(end, :) = obj.f_in_2(end, :);
                obj.f_in_7(end, :) = obj.f_in_6(end, :);
                obj.f_in_8(end, :) = obj.f_in_9(end, :);
            end
        end

        function obj=update_macro(obj, no_slip, free_slip, force)
            obj.rho = obj.f_in_1 + obj.f_in_2 + obj.f_in_3 + obj.f_in_4 + ...
                obj.f_in_5 + obj.f_in_6 + obj.f_in_7 + obj.f_in_8 + obj.f_in_9;
            vx_dot_f_in = ...
                obj.lattice.vx(1) * obj.f_in_1 + ...
                obj.lattice.vx(2) * obj.f_in_2 + ...
                obj.lattice.vx(3) * obj.f_in_3 + ...
                obj.lattice.vx(4) * obj.f_in_4 + ...
                obj.lattice.vx(5) * obj.f_in_5 + ...
                obj.lattice.vx(6) * obj.f_in_6 + ...
                obj.lattice.vx(7) * obj.f_in_7 + ...
                obj.lattice.vx(8) * obj.f_in_8 + ...
                obj.lattice.vx(9) * obj.f_in_9;
            vy_dot_f_in = ...
                obj.lattice.vy(1) * obj.f_in_1 + ...
                obj.lattice.vy(2) * obj.f_in_2 + ...
                obj.lattice.vy(3) * obj.f_in_3 + ...
                obj.lattice.vy(4) * obj.f_in_4 + ...
                obj.lattice.vy(5) * obj.f_in_5 + ...
                obj.lattice.vy(6) * obj.f_in_6 + ...
                obj.lattice.vy(7) * obj.f_in_7 + ...
                obj.lattice.vy(8) * obj.f_in_8 + ...
                obj.lattice.vy(9) * obj.f_in_9;

            % Normal forces on free-slip boundaries should not lead to
            % acceleration. Remove them.
            if free_slip(1); force.y(:, 1) = 0; end  % bottom
            if free_slip(2); force.y(:, end) = 0; end  % top
            if free_slip(3); force.x(1, :) = 0; end  % left
            if free_slip(4); force.x(end, :) = 0; end  % right

            obj.ux = (1 - no_slip) .* (vx_dot_f_in ./ obj.rho + force.x/2);
            obj.uy = (1 - no_slip) .* (vy_dot_f_in ./ obj.rho + force.y/2);
        end
        
        function obj=tick(obj, no_slip, free_slip, force, dt)
            obj = obj.compute_feq(obj.ux, obj.uy, obj.rho);
            obj = obj.collide(no_slip, free_slip, force, dt);
            obj = obj.stream();
        end

        function [rho, ux, uy]=get_macro(obj)
            rho = obj.rho;
            ux = obj.ux;
            uy = obj.uy;
        end

    end

    methods (Access = private)

        function obj=compute_feq(obj, ux, uy, rho)
            cs2 = obj.lattice.cs2;
            udotu = ux.*ux + uy.*uy;
            for vi=1:length(obj.lattice.vx)
                vdotu = ux*obj.lattice.vx(vi) + uy*obj.lattice.vy(vi);
                obj.(sprintf("f_eq_%i", vi)) = obj.lattice.w(vi) * rho .* ...
                    (1 + vdotu/cs2 + vdotu.^2/(2*cs2^2) - udotu/(2*cs2));
            end
        end

        function obj=collide(obj, no_slip, free_slip, force, dt)
            tau = obj.nu / obj.lattice.cs2 + .5;
            d2q9 = obj.lattice;
            for vi=1:length(d2q9.vx)
                f_in = obj.(sprintf("f_in_%i", vi));
                f_in_opp = obj.(sprintf("f_in_%i", d2q9.opp(vi)));
                f_eq = obj.(sprintf("f_eq_%i", vi));
                
                % Normal forces on free-slip boundaries should not lead to
                % acceleration. Remove them.
                if free_slip(1); force.y(:, 1) = 0; end  % bottom
                if free_slip(2); force.y(:, end) = 0; end  % top
                if free_slip(3); force.x(1, :) = 0; end  % left
                if free_slip(4); force.x(end, :) = 0; end  % right

                % from He et al., 1998
                forcedotvmu = force.x .* (d2q9.vx(vi)-obj.ux) + ...
                    force.y .* (d2q9.vy(vi)-obj.uy);
                Fi = (dt * (tau-.5))/tau * ...
                    forcedotvmu/d2q9.cs2 .* f_eq;

                obj.(sprintf("f_out_%i", vi)) = ...
                    (1 - no_slip) .* (f_in ...
                        - 1./tau .* (f_in - f_eq) ...
                        + Fi) + ...
                    no_slip .* f_in_opp;

            end
        end

        function obj=stream(obj)
            for vi=1:length(obj.lattice.vx)
                obj.(sprintf("f_in_%i", vi)) = circshift(...
                    obj.(sprintf("f_out_%i", vi)), [obj.lattice.vx(vi) obj.lattice.vy(vi) 0]);
            end
        end
    end
end
