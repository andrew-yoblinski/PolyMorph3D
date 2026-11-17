function formatAxes3D(ax, sz)
    % force boundary to original image dimensions
    xlim(ax, [1 sz(2)]);
    ylim(ax, [1 sz(1)]);
    zlim(ax, [1 sz(3)]);
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';

    % remove ticks
    ax.XTick = [];
    ax.YTick = [];
    daspect(ax, [1 1 1]);

    % border around XY volume bounds
    borderZ = 1;
    plot3(ax, [1 sz(2) sz(2) 1 1], [1 1 sz(1) sz(1) 1], [borderZ borderZ borderZ borderZ borderZ], ...
        'w-', 'LineWidth', 1.5);
    view(0, -90);
    daspect(ax, [1 1 1]);
end