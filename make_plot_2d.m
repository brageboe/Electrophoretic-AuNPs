% Written by Olle Hennert and Brage Bøe Svendsen
% makes log log plot by default
% res = resolution of the plot, e.g. res = 50 samples 50*50 points
function make_plot_2d(fun, xinterval, yinterval, plot_settings)
    arguments
        fun
        xinterval
        yinterval
        plot_settings.res = 50
        plot_settings.lim = 0
        plot_settings.x_log = true
        plot_settings.y_log = true
        plot_settings.zlabel = ""
        plot_settings.colorbar = true
        plot_settings.minvalue = 1e-12
        plot_settings.ticks = []
        plot_settings.cscale = false
        plot_settings.isoline1 = 1
        plot_settings.isoline2 = 1e3
    end
    res = plot_settings.res;
    
    if plot_settings.x_log == false
        x = linspace(xinterval(1), xinterval(2), res);
    else
        x = logspace(log10(xinterval(1)), log10(xinterval(2)), res);
    end
    if plot_settings.y_log == false
        y = linspace(yinterval(1), yinterval(2), res);
    else
        y = logspace(log10(yinterval(1)), log10(yinterval(2)), res);
    end
    
    z = zeros(res);
    for i = 1:res
        for j = 1:res
            value = max(fun(x(i), y(j)), plot_settings.minvalue);
            z(j, i) = value; 
        end
    end

    a = min(z(:));
    b = max(z(:));
    [a_ind_x, a_ind_y] = find(z(:,:)==a);
    disp("min = "+a+" @ ["+x(a_ind_x)+", "+y(a_ind_y)+"]")
    [b_ind_x, b_ind_y] = find(z(:,:)==b);
    disp("max = "+b+" @ ["+x(b_ind_x)+", "+y(b_ind_y)+"]")
    
    levels = logspace(log10(a), log10(b), 256);
    %contourf(x, y, z, levels, "LineColor", "none");
    if plot_settings.lim == 0
        lim = [a, b];
    else
        lim = plot_settings.lim;
    end
    
    imagesc(xinterval, yinterval, z, 'Interpolation', 'bilinear', lim);
    
    %färg intervall på plotten
    if plot_settings.colorbar
        %"northoutside"
        c = colorbar();
        c.Label.String = plot_settings.zlabel;
        %c.Position(4) = 1.5*c.Position(4);
        if ~isempty(plot_settings.ticks)
            c.Ticks = plot_settings.ticks;
        end
        %c.Ticks = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7];
        %c.Position = c.Position + 1e-10;
    end
    if plot_settings.cscale
        cmap = cm_inferno();
        clen = length(cmap);
        cportion = 0.7;
        cstart = round((1-cportion)*clen) + 1;
        cscaled = cmap(cstart:clen, :);
        colormap(cscaled);
    else
        colormap(cm_inferno());
        % colormap(cm_viridis());
        % colormap("turbo");
        % colormap("default");
    end
    
    set(gca,'ColorScale','log');

    set(gca,'YDir','normal')
    if plot_settings.x_log ~= false
        set(gca, 'XScale','log');
    end
    if plot_settings.y_log ~= false
        set(gca, 'YScale','log');
    end
    xlim(xinterval);
    ylim(yinterval);
    
    if num2str(plot_settings.isoline1) == "max"
        plot_settings.isoline1 = 0.9 * b; % isoline @ 90% of max
    end
    if num2str(plot_settings.isoline2) == "max"
        plot_settings.isoline2 = 0.9 * b; % isoline @ 90% of max
    end

    C = contourc(x, y, z, [plot_settings.isoline1 plot_settings.isoline1]); % draw isoline at [Fabs Fabs]
    hold on;
    if (length(C) > 1)
        N = C(2,1);
        start = 2;
        while true
            endpos = start + N - 1;
            contour_x = C(1, start:endpos);
            contour_y = C(2, start:endpos);
            plot(contour_x, contour_y, 'k--');
            
            if endpos >= length(C)
                break;
            end
            
            start = endpos + 2;
            N = C(2, endpos + 1);
        end
    end
    hold off;   
    C = contourc(x, y, z, [plot_settings.isoline2 plot_settings.isoline2]); % draw isoline at [Fabs Fabs]
    hold on;
    if (length(C) > 1)
        N = C(2,1);
        start = 2;
        while true
            endpos = start + N - 1;
            contour_x = C(1, start:endpos);
            contour_y = C(2, start:endpos);
            plot(contour_x, contour_y, 'k:');
            
            if endpos >= length(C)
                break;
            end
            
            start = endpos + 2;
            N = C(2, endpos + 1);
        end
    end
    hold off;   
end