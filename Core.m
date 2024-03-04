classdef Core
    properties
        Cp
        E
    end
    methods
        function obj=Core(T_0, Cp)
            obj.Cp = Cp;
            obj.E = T_0 * obj.Cp;
        end

        function T=get_T(obj)
            T = obj.E / obj.Cp;
        end
        
        function obj=tick(obj, Q, dt)
            % Positive q is outward heat flux.
            dE = -Q * dt;
            obj.E = obj.E + dE;
        end
    end
end
