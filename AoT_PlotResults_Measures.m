function [CORR_K_KL] = AoT_PlotResults_Measures(DK_RS,DKL_RS,CM2,n_regions,Title)

    % Comparison between univariate Kurtosis and KL divergence, across regions
    % One plot per paradigm
    figure;
    subplot(2,1,1);
    hold on;
    title(Title);
    h1 = bar(1:n_regions,mean(DK_RS));
    set(h1,'FaceColor',CM2(1,:),'EdgeColor','None');
    hold on;
    h2 = errorbar(1:n_regions,mean(DK_RS),abs(mean(DK_RS)-prctile((DK_RS),2.5)),abs(mean(DK_RS)-prctile((DK_RS),97.5)),'Color',CM2(1,:),'LineWidth',0.1,'LineStyle','None','CapSize',realmin);
    xlim([0,420]);
    set(gca,'Box','off');
    xlabel('Regions');
    ylabel('Regional \Delta kurtosis (F - B)');
    
    alpha = 0.2;   
    % Set transparency (undocumented)
    set([h2.Bar, h2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alpha])

    subplot(2,1,2);
    h3 = bar(1:n_regions,mean(DKL_RS));
    set(h3,'FaceColor',CM2(2,:),'EdgeColor','None');
    hold on;
    h4 = errorbar(1:n_regions,mean(DKL_RS),abs(mean(DKL_RS)-prctile((DKL_RS),2.5)),abs(mean(DKL_RS)-prctile((DKL_RS),97.5)),'Color',CM2(2,:),'LineWidth',0.1,'LineStyle','None','CapSize',realmin);
    xlim([0,420]);
    set(gca,'Box','off');
    xlabel('Regions');
    ylabel('Regional \Delta KL divergence (F - B)');
    
    % Set transparency (undocumented)
    set([h4.Bar, h4.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h4.Line.ColorData(1:3); 255*alpha])
    
    CORR_K_KL = corr(mean(DK_RS)',mean(DKL_RS)');

end