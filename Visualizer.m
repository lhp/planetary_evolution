classdef Visualizer < handle
    properties
        width;
        height;
        hT; hU;

        scaling
    end
    methods
        function obj=Visualizer(scaling)
            obj.scaling = scaling;
            obj.width = -1;
            obj.height = -1;
        end

        function obj=update_display(obj, macro, history)
            if obj.width < 1
                obj = obj.init_display(macro);
            end

            subplot(obj.height, obj.width, [1:3 5:7 9:11]);
            T_r = obj.scaling.to_T_r(macro.T);
            obj.hT = surf(T_r');
            obj.hU = streamslice(squeeze(macro.ux)', ...
                squeeze(macro.uy)', 'color', 'k');
            for i=1:length(obj.hU)
                set(obj.hU(i), 'zdata', 2*max(max(T_r))*ones(size(obj.hU(i).XData)));
            end
            colorbar('location', 'northoutside');
            shading interp
            view(2);
            [nx, ny] = size(T_r);
            xlim([1 nx]);
            ylim([1 ny]);
            
            subplot(obj.height, obj.width, [4 8 12]);
            yyaxis left
            plot(obj.scaling.to_T_r(mean(macro.T)), (1:ny)*obj.scaling.x/1e3);
            xlabel("T (K)");
            ylabel("y (km)");
                
            gy_in_s = 1e9 * (365*24*60*60);
            times_r = history.times * obj.scaling.t/gy_in_s;

            subplot(obj.height, obj.width, 13:16);
            yyaxis left
            plot(times_r, history.max_uy);
            ylabel("max uy (lattice units)");
            limleft = ylim;
            ticksleft = yticks;
            yyaxis right
            plot(times_r, history.max_uy);
            ylim(limleft);
            yticklabels(ticksleft * obj.scaling.x / obj.scaling.t * 1e2*365*24*60*60);
            ylabel("max uy (cm/year)");
            xlabel("t (Gy)");
            
            subplot(obj.height, obj.width, 17:20);
            plot(times_r, history.total_E/history.total_E(1), 'k-', 'DisplayName', 'total');
            hold on
            plot(times_r, history.mantle_E/history.mantle_E(1), 'r-', 'DisplayName', 'mantle');
            plot(times_r, history.core_E/history.core_E(1), 'b-', 'DisplayName', 'core');
            hold off
            yline(1, 'k-', 'HandleVisibility', 'off');
            yline(.8, 'k--', 'HandleVisibility', 'off');
            ylabel("E (normalized)");
            xlabel("t (Gy)");

            subplot(obj.height, obj.width, 21:24);
            scale_q = obj.scaling.k * obj.scaling.T / obj.scaling.x;
            plot(times_r, scale_q * history.q_cmb);
            yline(0);
            ylabel("q_cmb");
            xlabel("t (Gy)");
            
            drawnow;

        end
    end

    methods (Access=private)

        function obj=init_display(obj, macro)
            obj.width = 4;
            obj.height = 6;
            figure();
            set(gcf, 'position', [0 0 1000 800]);

            subplot(obj.height, obj.width, [1:3 5:7 9:11]);
            obj.hT = surf(obj.scaling.to_T_r(macro.T)');
            obj.hU = streamslice(squeeze(macro.ux)', ...
                squeeze(macro.uy)', 'color', 'k');
        end

    end
end
