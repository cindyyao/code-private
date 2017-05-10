%% 2016-11-22-0
for cc = 2:2
    figure(1)
    set(gcf, 'Position', [1 1 1000 500])
    id = ds_id(cc);
    for i = 1:5
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc}(5:end, 5:end, :);
    
            subplot(3,5,i)
            imagesc(sum(rf,3))
            colormap gray
            axis image
            axis off

            subplot(3,5,5+i)
            imagesc(rf(:,:,1))
            colormap gray

            axis image
            axis off

            subplot(3,5,10+i)
            imagesc(rf(:,:,2))
            colormap gray

            axis image
            axis off

        end
    end
    print_close(1,[24 12],num2str(id))
end

%% 2016-08-29-0
for cc = 7:7%length(ds_id)
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
            rf = rf_all{i}{cc};
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            title('on')
            axis image

            subplot(3,3,3*(i-1)+3)
            imagesc(rf(:,:,2))
            colormap gray
            title('off')
            axis image
        end
    end
    print_close(1,[12 12],num2str(id))
end

%% RF area ratio

load('DS161115.mat', 'rf_area_clean_all_mean')
rf_other = rf_area_clean_all_mean;
load('DS161122.mat', 'rf_area_clean_mean')
for onoff = 1:2
    rf_superior{onoff} = rf_area_clean_mean{onoff}(:, 1)';
    ratio{onoff} = rf_superior{onoff}./rf_other{onoff};
end

figure
plot(0:4, ratio{1}, 'b')
hold on
plot(0:4, ratio{2}, 'r')
plot([-1 5], [1 1], 'k--')
ylim([0 9])
ylabel('superior/others RF ratio')
xlim([-1 5])