classdef PlanetModel
    properties
        time
        history

        adv_diff
        navier_stokes
        core
        macro

        pp_l
        dTdt_l
        solid

        free_slip
        force_scale
        u_settling_tolerance
        velocity_skip
        dt
        ref_T_for_force;
    end
    methods
        function obj=PlanetModel(mantle_T_0, core_T_0, core_mantle_Cp_ratio, ...
                physical_props, dTdt, solid, total_time, velocity_skip, settling_tolerance)
            [nx, ny] = size(mantle_T_0);
            
            obj.macro.T = mantle_T_0;
            obj.macro.ux = zeros(nx, ny);
            obj.macro.uy = zeros(nx, ny);
            obj.macro.rho = physical_props.rho * ones(nx, ny);
            obj.adv_diff = AdvectionDiffusion(physical_props.kappa, obj.macro.ux, obj.macro.uy, obj.macro.T);
            obj.adv_diff = obj.adv_diff.set_dirichlet_BCs(0, core_T_0);
            obj.navier_stokes = NavierStokes(physical_props.nu, obj.macro.ux, obj.macro.uy, obj.macro.rho);

            % V * rho = mass --> M * C_p = total heat capacity
            mantle_cp = nx * ny * physical_props.rho * physical_props.cp;
            core_cp = core_mantle_Cp_ratio * mantle_cp;
            obj.core = Core(core_T_0, core_cp);

            obj.pp_l = physical_props;
            obj.dTdt_l = dTdt;
            obj.solid = solid;
            obj.solid(:, end) = 1;  % Ensure no-slip top boundary.

            obj.free_slip = [1 0 0 0];
            obj.force_scale = 1;
            obj.velocity_skip = velocity_skip;
            obj.u_settling_tolerance = settling_tolerance;
            obj.dt = 1;
            obj.ref_T_for_force = core_T_0;

            obj.time = 0;
            profile_every = ceil(total_time / 10000);  % At most 10,000 profiles
            obj.history = ModelHistory(profile_every, physical_props);
        end
        
        function obj=tick(obj)
            % Update adv-diff boundary condition and heating rate.
            obj.adv_diff = obj.adv_diff.set_dirichlet_BCs(...
                obj.adv_diff.T_top, obj.core.get_T());
            dTdt_i = obj.dTdt_l * .5^(obj.time/obj.pp_l.H_halflife);

            obj.adv_diff = obj.adv_diff.tick_all(obj.macro.ux, obj.macro.uy, dTdt_i, obj.dt);
            [obj.macro.T, obj.macro.dTdz] = obj.adv_diff.get_macro();
            
            Q_cmb = sum(-obj.pp_l.k * obj.macro.dTdz(:, 1));
            obj.core = obj.core.tick(Q_cmb, 1);
        
            if mod(obj.time, obj.velocity_skip) == 0
                force.y = obj.pp_l.alpha_g * obj.pp_l.rho * (obj.macro.T-obj.ref_T_for_force);
                force.x = 0 * force.y;

                [obj.navier_stokes, obj.force_scale] = settle_velocities(...
                    obj.navier_stokes, obj.solid, obj.free_slip, force, ...
                    obj.force_scale, obj.u_settling_tolerance, obj.dt);
                [obj.macro.rho, obj.macro.ux, obj.macro.uy] = obj.navier_stokes.get_macro();
                obj.macro.ux = obj.macro.ux / obj.force_scale;
                obj.macro.uy = obj.macro.uy / obj.force_scale;
            end

            obj.history = obj.history.update_history(obj.time, obj.macro, obj.core);
            obj.time = obj.time + obj.dt;
        end
    end
end
