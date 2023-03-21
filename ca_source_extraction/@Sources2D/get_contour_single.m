%% get contours of the all neurons
function Coor = get_contour_single(obj,thr,A_)
    d1 = obj.options.d1;
    d2 = obj.options.d2;
    %             tmp_kernel = strel('square', 3);
    for m=1
        % smooth the image with median filter
        A_temp = obj.reshape(full(A_(:, m)),2);
        % find the threshold for detecting nonzero pixels

        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend');
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        thr_a = A_temp(ind(ff));
        A_temp = obj.reshape(A_temp,2);

        % crop a small region for computing contours
        [tmp1, tmp2, ~] = find(A_temp);
        if isempty(tmp1)
            Coor{m} = zeros(2,1);
            continue;
        end
        rmin = max(1, min(tmp1)-3);
        rmax = min(d1, max(tmp1)+3);
        cmin = max(1, min(tmp2)-3);
        cmax = min(d2, max(tmp2)+3);
        A_temp = A_temp(rmin:rmax, cmin:cmax);

        if nnz(A_temp)>36
            l = bwlabel(medfilt2(A_temp>thr_a));
        else
            l = bwlabel(A_temp>=thr_a); 
        end
        l_most = mode(l(l>0));
        if isnan(l_most)
            Coor{m} = zeros(2, 1); 
            continue; 
        end
        ind = (l==l_most);
        A_temp(ind) =  max(A_temp(ind), thr_a);
        A_temp(~ind) = min(A_temp(~ind), thr_a*0.99);

        pvpairs = { 'LevelList' , thr_a, 'ZData', A_temp};
        h = matlab.graphics.chart.primitive.Contour(pvpairs{:});
        temp = h.ContourMatrix;
        if isempty(temp)
            temp = obj.get_contours((thr+1)/2, ind_show(m));
            Coor{m} = temp{1};
            continue;
        else
            temp(:, 1) = temp(:, 2);
            temp = medfilt1(temp')';
            temp(:, 1) = temp(:, end);
            Coor{m} = bsxfun(@plus, temp, [cmin-1; rmin-1]);
        end

    end
end
        