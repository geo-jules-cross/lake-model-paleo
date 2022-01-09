clear ;
%% Load Data

load('DATA/lake.mat');

%%
figure(2); clf; clear ha; ha = tight_subplot(3,1, [0.03 0.05], [.07 .05], [.06 .08]);
set(gcf,'units','normalized','outerposition',[0.1 0.0 0.8 1])

j = 0; % counter for plotting

datum = 1; % plot relative to datum (2015 is 20th year in series)

name = ["Bonney", "Hoare", "Fryxell"];

for l = 1:3

% plot lake-level history h(t)
% ----------------------------
    j = j+1;
    axes(ha(j)); hold on; box on; grid on;
    
    yyaxis left
    plot(lake(l).t_vec, lake(l).h-lake(l).h(datum), 'linewidth', 1.2, 'linestyle','-', 'color', [0 0 0])
    if(j==2)
        ylabel('Lake-level change [m]')
    end
    %ylim auto
    if(j==1)
        ylim([-15 35]);
    elseif(j==2)
        ylim([-5 15]);
    elseif(j==3)
        ylim([-10 10]);
    end
    yticklabels('auto')
    set(ha(1:3),'YColor', 'k')
    % Zero Line
    line([1990, 2300],[0, 0], 'Color','k', 'LineWidth', 1, 'HandleVisibility','off')
    left_ylim = ylim;
    left_yticks = yticks;
    
    yyaxis right
    plot(lake(l).t_vec, lake(l).h, 'linewidth', 1.2, 'linestyle','--', 'color', [0 0 0])
    if(j==2)
        ylabel('Lake surface height [m asl]')
    end
    right_ylim = left_ylim + lake(l).h(datum);
    ylim(left_ylim + lake(l).h(datum))
    yticks(left_yticks + lake(l).h(datum))
    yticklabels('auto')
    ytickformat('%.0f')

    % general axes parameters
    title(name(j))
    set(ha(1:3),'XColor','k', 'YColor', 'k', 'FontWeight', 'bold', 'LineWidth', 1.25, 'FontSize', 14, 'GridColor', 'k')
    set(ha(1:3), 'XGrid', 'off');
    xlim([1990, 2300])
    xticks((1990:25:2300))
    if(j==3)
        xlabel('Date')
        xticklabels('auto')
    end
    
end