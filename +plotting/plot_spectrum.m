% figure;
hold on;

% title('Amplitude');
group_symbols = {'v', 'o'};
group_colors = {[83 148 255]/255, [255 117 117]/255};


figure;
set(gcf,'color','w');
    
    
    
    
        ylabel('Amplitude ( microvolt ^2 )', 'FontSize', 14);
        xlabel('Frequency', 'FontSize', 14);
%         ylim([0, 1]);
        xlim([1 20]);
        hold on;
        
        x = spectrum_freqs;
        
        y = nanmean( cat(1, spectrum{:, 1}), 1);
        
        e = nansem( cat(1, spectrum{:, 1}), 1);
        
        plot(x, y, 'LineWidth', 1, 'Color', group_colors{1});
        
        x_fill = [x, fliplr(x)];
        y_fill = [y-e, fliplr(y+e)];
        
        fill( x_fill, y_fill, group_colors{1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none' );
        
    




% Overlay the harmonic frequencies

% set(gca, 'XTick', [28.8, 36]);