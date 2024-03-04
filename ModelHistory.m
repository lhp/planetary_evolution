classdef ModelHistory < handle
    % It is necessary to derive from handle so that this object is not
    % copied every time it is modified.
    properties
        pp_l
        time_i
        profile_i

        times

        mean_T
        T_surface
        T_cmb
        core_T

        q_surface
        q_cmb
        max_q_adv
        max_mean_q_adv

        max_uy
        max_u

        mantle_E
        core_E
        total_E

        profile_times
        profile_every

        T_profiles
        q_profiles
        q_adv_profiles
    end
    methods
        function obj=ModelHistory(profile_every, pp_l)
            obj.pp_l = pp_l;
            obj.time_i = 0;
            obj.profile_i = 0;
    
            obj.times = [];
    
            obj.mean_T = [];
            obj.T_surface = [];
            obj.T_cmb = [];
            obj.core_T = [];
    
            obj.q_surface = [];
            obj.q_cmb = [];
            obj.max_q_adv = [];
            obj.max_mean_q_adv = [];
    
            obj.max_uy = [];
            obj.max_u = [];
    
            obj.mantle_E = [];
            obj.core_E = [];
            obj.total_E = [];
    
            obj.profile_times = [];
            obj.profile_every = profile_every;
    
            obj.T_profiles = [];
            obj.q_profiles = [];
            obj.q_adv_profiles = [];
        end

        function obj=update_history(obj, time, macro, core)
            obj.time_i = obj.time_i + 1;
            obj = obj.grow_if_full();

            obj.times(obj.time_i) = time;

            T_prof = mean(macro.T);
            obj.mean_T(obj.time_i) = mean(T_prof);
            obj.T_surface(obj.time_i) = mean(macro.T(:, end));
            obj.T_cmb(obj.time_i) = mean(macro.T(:, 1));
            obj.core_T(obj.time_i) = core.get_T();
    
            q = -obj.pp_l.k * mean(macro.dTdz);
            obj.q_surface(obj.time_i) = mean(q(:, end));
            obj.q_cmb(obj.time_i) = mean(q(:, 1));
            E = macro.T * obj.pp_l.rho * obj.pp_l.cp;  % assumes dx=1
            obj.max_q_adv(obj.time_i) = max(max(abs(E .* macro.uy)));
            obj.max_mean_q_adv(obj.time_i) = max(mean(E .* macro.uy));

            obj.max_uy(obj.time_i) = max(max(abs(macro.uy)));
            obj.max_u(obj.time_i) = max(max(sqrt(macro.ux.^2 + macro.uy.^2)));

            obj.mantle_E(obj.time_i) = sum(E, 'all');
            obj.core_E(obj.time_i) = core.E;
            obj.total_E(obj.time_i) = obj.mantle_E(obj.time_i) + obj.core_E(obj.time_i);


            if mod(obj.time_i, obj.profile_every) == 0
                obj.profile_i = obj.profile_i + 1;
                obj.profile_times(obj.profile_i) = time;

                obj.T_profiles(:, obj.profile_i) = mean(macro.T);
                obj.q_profiles(:, obj.profile_i) = -obj.pp_l.k * mean(macro.dTdz);
                obj.q_adv_profiles(:, obj.profile_i) = mean(E .* macro.uy);
            end
        end

    end

    methods (Access = private)
        function obj=grow_if_full(obj)
            % Necessary for reasonable amortized cost of dynamic growth.
            if ~isempty(obj.times) && ~isnan(obj.times(end))
                grow_by = round(0.25 * length(obj.times));
                obj.times = [obj.times, nan(1, grow_by)];
                obj.mean_T = [obj.mean_T, nan(1, grow_by)];
                obj.T_surface = [obj.T_surface, nan(1, grow_by)];
                obj.T_cmb = [obj.T_cmb, nan(1, grow_by)];
                obj.core_T = [obj.core_T, nan(1, grow_by)];
        
                obj.q_surface = [obj.q_surface, nan(1, grow_by)];
                obj.q_cmb = [obj.q_cmb, nan(1, grow_by)];
                obj.max_q_adv = [obj.max_q_adv, nan(1, grow_by)];
                obj.max_mean_q_adv = [obj.max_mean_q_adv, nan(1, grow_by)];
        
                obj.max_uy = [obj.max_uy, nan(1, grow_by)];
                obj.max_u = [obj.max_u, nan(1, grow_by)];
        
                obj.mantle_E = [obj.mantle_E, nan(1, grow_by)];
                obj.core_E = [obj.core_E, nan(1, grow_by)];
                obj.total_E = [obj.total_E, nan(1, grow_by)];
            end

            if ~isempty(obj.profile_times) && ~isnan(obj.profile_times(end))
                grow_by = round(0.25 * length(obj.profile_times));
                obj.profile_times = [obj.profile_times, nan(1, grow_by)];

                profile_height = length(obj.T_profiles(:, 1));
                obj.T_profiles = [obj.T_profiles, nan(profile_height, grow_by)];
                obj.q_profiles = [obj.q_profiles, nan(profile_height, grow_by)];
                obj.q_adv_profiles = [obj.q_adv_profiles, nan(profile_height, grow_by)];
            end
        end
    end
end
