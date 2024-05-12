function horizontalScroll(src, event)
    ax = gca;
    scrollSpeed = 1; % Adjust the scroll speed here
    currentXLim = ax.XLim;
    rangeX = diff(currentXLim);
    
    if event.VerticalScrollCount > 0
        % Scrolling down, move right
        ax.XLim = currentXLim + rangeX * scrollSpeed * 0.05;
    elseif event.VerticalScrollCount < 0
        % Scrolling up, move left
        ax.XLim = currentXLim - rangeX * scrollSpeed * 0.05;
    end
end