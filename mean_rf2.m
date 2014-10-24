function [sta_f_weighted, sta_f_n, sta_c] = mean_rf2(datarun0, datarun1, cell_type, cell_id, cell_idx)
% [sta_f_weighted, sta_f_n, sta_c, fit_sta_c, params_c] = mean_rf2(datarun0, datarun1, cell_type, cell_id, cell_idx)
d0 = size(datarun0.stas.stas{1});
d1 = size(datarun1.stas.stas{1});
f0 = 1200/d0(1);
f1 = 1200/d1(1);

n = size(cell_id, 1);
sta_f_weighted = cell(n, 2);
sta_f_n = zeros(n, 2);
for k = 1:n;
    index = cell_idx{k};
    cell_numb = size(cell_id{k}, 1);
    sta_f0 = zeros(f0*d0(1)/2, f0*d0(2)/2, d0(3), d0(4));
    sta_f1 = zeros(f1*d1(1)/2, f1*d1(2)/2, d0(3), d1(4));
    sta_f_n0 = 0;
    sta_f_n1 = 0;
    if isempty(index) == 0
       for i = 1:cell_numb
           sta = datarun0.stas.stas{index(i, 1)};
           if cell_type{k}(2) == 'N' 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d0(1)*d0(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f0*d0(1), f0*d0(2), d0(3), d0(4));
           for j = 1:d0(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f0, f0));
           end
           xy = floor(f0*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f0 = sta_f0 + sta_temp(xy(2)-f0*d0(1)/4+1:xy(2)+f0*d0(1)/4, xy(1)-f0*d0(2)/4+1:xy(1)+f0*d0(2)/4, :, :)/noise;
           sta_f_n0 = sta_f_n0 + 1;
           end
           
    
           sta = datarun1.stas.stas{index(i, 2)};
           if cell_type{k}(2) == 'N' 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d1(1)*d1(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f1*d1(1), f1*d1(2), d1(3), d1(4));
           for j = 1:d1(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f1, f1));
           end
           xy = floor(f1*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f1 = sta_f1 + sta_temp(xy(2)-f1*d1(1)/4+1:xy(2)+f1*d1(1)/4, xy(1)-f1*d1(2)/4+1:xy(1)+f1*d1(2)/4, :, :)/noise;
           sta_f_n1 = sta_f_n1 + 1;
           end
           
       end       
    end
    sta_f_weighted{k, 1} = sta_f0;
    sta_f_weighted{k, 2} = sta_f1;
    sta_f_n(k, 1) = sta_f_n0;
    sta_f_n(k, 2) = sta_f_n1;
end

fprintf('weighted addition done')

sta_c = cell(n, 2);
for i = 1:n
        sta = sta_f_weighted{i, 1};
        sta_c_temp = zeros(d0(1), d0(2), d0(3), d0(4));
        for j = 1:d0(1)
            for k = 1:d0(2)
                for m = 1:d0(4)
                    sta_temp = sta(f0/2*(j-1)+1:f0/2*j, f0/2*(k-1)+1:f0/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 1} = sta_c_temp;
    
    sta = sta_f_weighted{i, 2};
    sta_c_temp = zeros(d1(1), d1(2), d1(3), d1(4));
        for j = 1:d1(1)
            for k = 1:d1(2)
                for m = 1:d1(4)
                    sta_temp = sta(f1/2*(j-1)+1:f1/2*j, f1/2*(k-1)+1:f1/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 2} = sta_c_temp;
end

fprintf('down scale done')


end

