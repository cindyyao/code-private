% cell_type = {'on brisk transient', 'on transient', 'off brisk transient', 'off transient', 'off transient large'};
% % cell_type = {'on', 'off'};
% n = length(cell_type);
% cell_id = cell(n, 1);
% cell_idx = cell(n, 1);
% for i = 1:n
%     idx = get_cell_indices(datarun{1}, cell_type{i});
%     id = datarun{1}.cell_ids(idx);
%     
%     for j = 1:length(idx)
%         a = datarun{2}.cell_ids - id(j);
%         if isempty(find(a == 0)) == 1
%             id(j) = 0;
%             idx(j) = 0;
%         end
%     end
%     id(id == 0) = [];
%     idx(idx == 0) = [];
%     
%     cell_idx{i} = idx;
%     cell_id{i} = id;        
% end
% 


figure
for j = 1:15
    sta = sta_c{1, 1}(:, :, :, 12);
    im = norm_image(sta);
    image(im)
    pause
end




       
%%

NDF = [0 1];
figure
for nn = 2:6
for k = 1:2
    
sta = sta_c{nn, k};
d = size(sta, 1);

if cell_type{nn}(2) == 'N' 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);


if cell_type{nn}(2) == 'N' 
    com = rf_com(sta_1);
else
    com = rf_com(-sta_1);
end


pix = [1:d] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(d, d);
for i = 1:d
    for j = 1:d
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end

dis1 = reshape(dis, 1, d^2);
if k ~= 1
    if cell_type{nn}(2) == 'N'  
        r = max(sta_temp1(:))/max(sta_1(:));
    else
        r = min(sta_temp1(:))/min(sta_1(:));
    end
end
    
sta_temp = reshape(sta_1, 1, d^2);

if k == 1
    sta_temp1 = sta_temp;
end

if k == 1
    subplot(2, 3, nn-1);
    if cell_type{nn}(2) == 'N' 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7)
        hold on 
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    end
else
    if cell_type{nn}(2) == 'N' 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    end
end


end
plot(dis1, zeros(1, d^2), 'm')
xlim([0 20])
legend(['NDF ' num2str(NDF(1))], ['NDF ' num2str(NDF(2))])
title(cell_type{nn});

end

