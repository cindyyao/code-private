function com = ei_com_xy(ei,positions,ROI,mode)


% threshold for significant electrodes, in EI sigmas
sig_elec_thresh = 5;

switch mode
    case 'p2t'
        % identify ROI
        % identify electrode with largest amplitude signal
        [junk,peak] = max(max(ei, [], 2) - min(ei, [], 2));
        % get distance of each electrode from this one
        dists = sqrt(sum((positions - repmat(positions(peak,:),size(positions,1),1)).^2,2));
        % only use electrodes within the radius
        roi = dists < ROI;

        thresh = sig_elec_thresh*robust_std(ei);

        % initialize

        [max_V, maxi] = max(ei, [], 2);
        max_V(max_V'<=thresh(maxi)) = 0;

        [min_V, mini] = min(ei, [], 2);
        min_V(min_V'>=-thresh(mini)) = 0;
        ei_amp = max_V(roi) - min_V(roi);
        elec_roi = find(roi == 1);
        com = [ei_amp'*positions(elec_roi,1)/sum(ei_amp)...
                ei_amp'*positions(elec_roi,2)/sum(ei_amp)];
    case 'neg'
        % identify ROI
        % identify electrode with largest amplitude signal
        [junk,peak] = min(min(ei, [], 2));
        % get distance of each electrode from this one
        dists = sqrt(sum((positions - repmat(positions(peak,:),size(positions,1),1)).^2,2));
        % only use electrodes within the radius
        roi = dists < ROI;
        elec_roi = find(roi == 1);
        % ensure that the roi electrodes are central-symmetric about peak
        % electrod.
        roi_positions = positions(elec_roi, :);
        peak_position = positions(peak, :);
        i = 1;
        while i <= length(elec_roi)
            sym_pos = peak_position*2 - roi_positions(i, :);
            if isempty(find(sum((roi_positions - repmat(sym_pos, length(elec_roi), 1)).^2, 2) == 0))
                elec_roi(i) = [];
                roi_positions(i, :) = [];
            else
                i = i+1;
            end
        end
        
        thresh = sig_elec_thresh*robust_std(ei);

        % initialize

        [min_V, mini] = min(ei, [], 2);
        min_V(min_V'>=-thresh(mini)) = 0;
        ei_amp = -min_V(elec_roi);
        com = [ei_amp'*positions(elec_roi,1)/sum(ei_amp)...
                ei_amp'*positions(elec_roi,2)/sum(ei_amp)];

end

