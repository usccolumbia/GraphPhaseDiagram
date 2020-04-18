function xrd_peak_path = peakfinder(xrd_path)
    set(0, 'DefaultFigureVisible', 'off');
    
    [pathstr,name,ext] = fileparts(xrd_path);
    xrd = csvread(xrd_path);
    twotheta = xrd(1, :);
    sample_num = size(xrd,1) - 1;
    feature_num = size(xrd,2);
    xrd = xrd(2:sample_num + 1, :);
    x = 1:feature_num;
    xrd_peakdata = zeros(sample_num + 1, feature_num);
    xrd_peakdata(1, :) = twotheta;
%     xrd_peakdata_filter = zeros(sample_num, feature_num);
    xrd_peakcenter_location = zeros(sample_num, feature_num);
    start_sample = 1;
    stop_sample = sample_num;

    for sample = start_sample:stop_sample
        close all;
        y = xrd(sample, :);
        
        sorted_y = sort(y);
        Amp = sorted_y((floor(feature_num * 0.5)));
        
        % for loop
        min_peak = inf;
        loop_continue = 1;
        peak_num_list = inf(100, 1);
        peak_width_list = zeros(100, 1);
        peak_area_list = zeros(100, 1);
        peak_height_list = zeros(100, 1);
        final_p = [];
        
        for WP = 42
            Slope = 0.7 * WP ^ -2;
            Smooth = WP/2.0;
            Fit = WP/2.0;
            p = findpeaksb(x, y, Slope, Amp, Smooth, Fit, 3, WP*2, 1, 0, 1, 1);

            height_sum = 0;
            for peak = 1:size(p, 1)
                if p(peak, 2) >= 1 && p(peak, 2) <= feature_num && sum(p(:, 4)) ~= 0
                    height_sum = height_sum + xrd(sample, round(p(peak, 2)));
                end
            end
            
            % for loop
            peak_height_list(WP, 1) =  height_sum/size(p, 1);
            peak_width_list(WP, 1) = sum(p(:, 4))/size(p, 1);
            peak_area_list(WP, 1) = sum(p(:, 5));
            peak_num = size(p, 1);
            peak_num_list(WP, 1) = peak_num; 

            if peak_num < min_peak
                WP_final = WP;
                min_peak = peak_num;
                final_p = p;
            end
            
        end
        
        % loop 
        [~, WP_final] = min(peak_num_list);
        %[~, WP_final] = max(peak_height_list);
        min_peak = peak_num_list(WP_final);
        Slope = 0.7 * WP_final ^ -2;
        Smooth = WP_final/2.0;
        Fit = WP_final/2.0;
        p = findpeaksb(x, y, Slope, Amp, Smooth, Fit, 3, WP_final*2, 1, 0, 3, 1);
        original_p = p;
        
        p = round(p);

        for peak = 1:size(p, 1)
            if p(peak, 2) <= feature_num && p(peak, 2) >= 1
                if p(peak, 2) - p(peak, 4) >= 1 && p(peak, 2) + p(peak, 4) <= feature_num
                    peak_range = (p(peak, 2) - p(peak, 4)):(p(peak, 2) + p(peak, 4));
                elseif p(peak, 2) - p(peak, 4) < 1
                    peak_range = 1:(p(peak, 2) + p(peak, 4));
                elseif p(peak, 2) + p(peak, 4) > feature_num
                    peak_range = (p(peak, 2) - p(peak, 4)):feature_num;
                end
                xrd_peakdata(sample + 1, peak_range) = xrd(sample, peak_range);
                % Another method is do a second-order polynomial approximation
                % poly = polyfit(peak_range, xrd(sample, peak_range), 3);
                % xrd_peakdata(sample, peak_range) = polyval(poly, peak_range);
                [val, idx] = max(xrd_peakdata(sample + 1, peak_range));
                xrd_peakcenter_location(sample, peak_range) = ones(1, length(peak_range)) * peak_range(idx);
            end
        end

%         fig = figure;
%         subplot(2,1,1);
%         plot(x,y);
%         text(original_p(:,2),xrd(sample, round(abs(original_p(:, 2)))),num2str(original_p(:,1)))
%         subplot(2,1,2);
%         y = xrd_peakdata(sample, :);
%         plot(x,y);
%         filename = [pathstr + '/xrd_spectrum/', num2str(sample), '.jpg'];
%     	  saveas(fig, filename);
        
    end
    csvwrite(strcat(pathstr,'/FePdGa_XRD_peak.csv'), xrd_peakdata);
    xrd_peak_path = strcat(pathstr,'/FePdGa_XRD_peak.csv');
end