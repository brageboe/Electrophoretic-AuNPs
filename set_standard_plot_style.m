function set_standard_plot_style()
    fig=gcf;
    fig.Position = [100, 300, 600, 400]; %Height 400 normal, 350 colorbar, %Width 600 normal, 800 colorbar
    
    fontsize = 16; %16 normal, 20 colorbar, 28 small

    set(0,'DefaultLineLineWidth',1.8);
    set(0,'DefaultFunctionLineLineWidth',1.8);
    set(0, 'DefaultLegendFontSizeMode','manual');
    set(0, 'DefaultLegendFontSize', fontsize);
    
    set(0,'defaulttextInterpreter','latex');
    set(0,'defaultAxesTickLabelInterpreter','latex');
    set(0,'defaultLegendInterpreter','latex');

    set(0,'defaultAxesXTickLabelRotationMode','manual')
    set(0,'defaultAxesYTickLabelRotationMode','manual')
    set(0,'defaultAxesZTickLabelRotationMode','manual')

    set(gca,'box','on');
    set(gca,'fontname','Arial');
    set(gca,'fontsize',fontsize);
    set(gca,'linewidth',1); % Normal 1, small 1.5
    set(gca,'TickLength',[0.02, 0.01]); % Normal 0.02, small 0.03
    
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    %legend('boxoff');
    colororder(["#004488", "#994455", "#997700", "#6699CC"]);
    %colororder(["#4477AA", "#AA3377", "#228833", "#66CCEE", "#CCBB44"]);
    %colororder(["#4477AA", "#AA3377", "#66CCEE", "#228833", "#CCBB44"])

    C = findall(gcf,'type','ColorBar');
    if ~isempty(C)
        c = C(1);
        c.Label.FontSize = 18; %18 normal, 22 colorbar only
        c.Label.Interpreter = 'latex';
        c.TickLabelInterpreter = 'latex';
        c.LineWidth = 1;
        c.TickLength = 0.02;
        c.Ruler.TickLabelRotation = 0;
    end
end