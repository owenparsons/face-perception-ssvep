

titles{1} = 'Right Hemisphere';
titles{2} = 'Left Hemisphere';
set(gcf,'color','w');

for hemisphere = 1:2
    
    ylabel('Amplitude (\muV^2)', 'FontSize', 14);
    
    subplot(1, 2, hemisphere);
    hold on;
    title(titles{hemisphere});
    xlim([1 13]);
    xlabel('Frequency', 'FontSize', 14);
    set(currfig, 'XTick', 0:1.2:13);
    
    ymat = cat(1, spectrum{:, hemisphere});
    ymat(ymat==0)=NaN;
    
    y = nanmean(ymat, 1);
    x = fft_data.freq;
    e = nansem(ymat, 1);
    
    plot(x, y, 'LineWidth', 1, 'Color', [83 148 255]/255);
    
    fill_index = ~isnan(y);
    x_fill = [x(fill_index), fliplr(x(fill_index))];
    y_fill = [y(fill_index)-e(fill_index), fliplr(y(fill_index)+e(fill_index))];
    
    fill( x_fill, y_fill, [83 148 255]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none' );
    
    
end

suptitle('Amplitude Spectrum');