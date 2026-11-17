function formatAxesWithScaleBar(ax, sz, barX1, barX2, barY, barLengthUm)
    % force boundary to original image dimensions
    xlim(ax, [1 sz(2)]);
    ylim(ax, [1 sz(1)]);
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';

    % remove ticks
    ax.XTick = [];
    ax.YTick = [];
    daspect(ax, [1 1 1]);
    hold(ax, 'on'); 

    % add scale bar
    plot(ax, [barX1 barX2], [barY barY], 'w-', 'LineWidth', 4);
    text(ax, mean([barX1 barX2]), barY - 0.03*sz(1), ...
        sprintf('%d \\mum', barLengthUm), ...
        'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 12);
end