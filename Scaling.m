classdef Scaling
    properties
        kappa
        nu
        x       % L
        rho     % M / L^3
        cp      % E / (M * T)
        T       % T
        T_offset

        k       % E / (t * L * T)
        q       % E / (t * L^2)
        H       % E / (t * M)

        t
        alpha_g
    end
    methods
        function obj=Scaling(x_scale, kappa_scale, nu_scale, rho_scale, cp_scale, T_scale, T_offset)
            obj.kappa = kappa_scale;
            obj.nu = nu_scale;
            obj.x = x_scale;
            obj.rho = rho_scale;
            obj.cp = cp_scale;
            obj.T = T_scale;
            obj.T_offset = T_offset;

            % Two derived scales.
            obj.t = obj.x ^2 / obj.kappa;
            % equivalent to nu / (T*t*x), and larger alpha_g scale gives
            % smaller lattice alpha_g.
            obj.alpha_g = obj.kappa * obj.nu / (obj.T * obj.x^3);

            % Compound useful scales.
            obj.k = obj.kappa * obj.rho * obj.cp;
            obj.q = obj.k * obj.T / obj.x;
            obj.H = obj.cp * obj.T / obj.t;
        end

        function T_r=to_T_r(obj, T_l)
            T_r = T_l * obj.T + obj.T_offset;
        end
        function T_l=to_T_l(obj, T_r)
            T_l = (T_r - obj.T_offset) / obj.T;
        end
    end
end
