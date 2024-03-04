classdef AdvectionDiffusion
    properties
        lattice
        f_in
        f_out
        f_eq
        vx
        vy
        kappa

        T_top = nan
        T_bot = nan

        % Macro variables
        T
        dTdz
    end

    methods
        function obj=AdvectionDiffusion(kappa, ux_0, uy_0, T_0)
            d2q5 = get_lattice_d2q5();
            obj.lattice = d2q5;
            obj.kappa = kappa;
            obj = obj.compute_feq(ux_0, uy_0, T_0);
            obj.f_in = obj.f_eq;
            obj.f_out = zeros(size(obj.f_in));
            obj.T = T_0;
            obj.dTdz = zeros(size(T_0));
        end

        function obj=set_dirichlet_BCs(obj, T_top, T_bot)
            obj.T_top = T_top;
            obj.T_bot = T_bot;
        end

        function obj=tick_all(obj, ux, uy, heating_dTdt, dt)
            if ~isnan(obj.T_top)
                obj = obj.enforce_bc_top_T(obj.T_top);
            end
            if ~isnan(obj.T_bot)
                obj = obj.enforce_bc_bot_T(obj.T_bot);
            end

            obj = obj.update_macro(dt);

            obj = obj.tick(ux, uy, heating_dTdt, dt);
        end

        function [T, dTdz]=get_macro(obj)
            T = obj.T;
            dTdz = obj.dTdz;
        end

    end

    methods (Access = private)

        function obj=enforce_bc_top_T(obj, T_top)
            obj.f_in(:, end, 5) = ...
                T_top - sum(obj.f_in(:, end, [1 2 3 4]), 3);
        end
        function obj=enforce_bc_bot_T(obj, T_bot)
            obj.f_in(:, 1, 3) = ...
                T_bot - sum(obj.f_in(:, 1, [1 2 4 5]), 3);
        end

        function obj=update_macro(obj, dt)
            obj.T = sum(obj.f_in, 3);

            tau = obj.kappa / obj.lattice.cs2 + .5;
            obj.dTdz = ((obj.f_in(:, :, 3) - obj.f_eq(:, :, 3)) - ...
                (obj.f_in(:, :, 5) - obj.f_eq(:, :, 5))) / ...
                (- tau * dt * obj.lattice.cs2);
        end

        function obj=tick(obj, ux, uy, heating_dTdt, dt)
            obj = obj.compute_feq(ux, uy, obj.T);
            obj = obj.collide(heating_dTdt, dt);
            obj = obj.stream();
        end

        function obj=compute_feq(obj, ux, uy, T)
            cs2 = obj.lattice.cs2;
            for vi=1:length(obj.lattice.vx)
                vdotu = ux*obj.lattice.vx(vi) + uy*obj.lattice.vy(vi);
                obj.f_eq(:, :, vi) = obj.lattice.w(vi) * T .* ...
                    (1 + vdotu/cs2);
            end
        end

        function obj=collide(obj, heating_dTdt, dt)
            tau = obj.kappa / obj.lattice.cs2 + .5;
            for vi=1:length(obj.lattice.vx)
                obj.f_out(:, :, vi) = obj.f_in(:, :, vi) - ...
                    1/tau * (obj.f_in(:, :, vi) - obj.f_eq(:, :, vi)) + ...
                    obj.lattice.w(vi) * heating_dTdt * dt;
            end
        end

        function obj=stream(obj)
            for vi=1:length(obj.lattice.vx)
                obj.f_in(:, :, vi) = circshift(...
                    obj.f_out(:, :, vi), [obj.lattice.vx(vi) obj.lattice.vy(vi) 0]);
            end
        end
    end
end
