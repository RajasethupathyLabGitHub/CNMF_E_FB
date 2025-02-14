%% manually pick neurons
axes(ax_selected);
xrange = get(gca, 'xlim');
yrange = get(gca, 'ylim');
title('click pixels out of the field to end');

while true
    [tmp_x, tmp_y] = ginput(1);
    if tmp_x<xrange(1) || tmp_x>xrange(2) || (tmp_y<yrange(1)) || (tmp_y>yrange(2))
        break;
    else
        
        
        % compute distance to all neurons
        if isempty(IDs)
            return;
        end
        % find the closest neuron
        dist = (tmp_x- ctr(IDs,2)).^2 + (tmp_y - ctr(IDs,1)).^2;
        [~, ind] = min(dist);
        IDs(ind) = [];
        
        %% show spatial shapes of the selected neurons
        img_selected = zeros(size(AA,1), 3);
        col = col0;
        for m=1:3
            img_selected(:, m) = sum(bsxfun(@times, AA(:, IDs), mod(col(:, IDs), 2)), 2);
            col = floor(col/2);
        end
        img_selected = neuron.reshape(img_selected, 2);
        img_selected = img_selected/max(img_selected(:))*(2^16);
        img_selected = uint16(img_selected);
        
        axes(ax_selected);
        imagesc(img_selected);
        axis equal off tight;
        xlim([xmin-gSiz, xmax+gSiz]);
        ylim([ymin-gSiz, ymax+gSiz]);
        
        %% show temporal traces of the selected neurons
        tmp_C = neuron.C_raw(IDs, :)';
        tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
        axes(ax_trace); cla; hold on; 
        for m=1:length(IDs)
            plot(tmp_C(:,m), 'color', color_all(IDs(m),:),  'linewidth', 2);
        end
    end
    
end


