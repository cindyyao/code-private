%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
load('DS161017.mat')
datawn = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/data001-map/data001-map', opt);
datawn = load_ei(datawn, cell2mat(id_dir));

%% white noise cross correlation
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
ct = 1;
N = 10000;
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
corr_cells_test = [];
dis = [];
area_ccf = [];
area_ccf_min = [];
pos = datawn.ei.position;
mode = 'neg';

for c1 = 12:length(id_dir{ct})-1
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 1980 1080])
    for c2 = c1+1:length(id_dir{ct})
        if c1 == 11 & c2 == 22
            continue
        end
        if c1 ~= c2
            id1 = id_dir{ct}(c1);
            id2 = id_dir{ct}(c2);
            idx1 = get_cell_indices(datawn, id1);
            idx2 = get_cell_indices(datawn, id2);
            spikes1 = datawn.spikes{idx1};
            spikes1_TF= ceil(spikes1/bin_size);
            spikes1 = zeros(duration/bin_size, 1);
            spikes1(spikes1_TF) = 1;
            
            spikes2 = datawn.spikes{idx2};
            spikes2_TF= ceil(spikes2/bin_size);
            spikes2 = zeros(duration/bin_size, 1);
            spikes2(spikes2_TF) = 1;
            
            A = xcorr(spikes1, spikes2, max_lag);
            [maxv(c1, c2), maxi(c1, c2)] = max(A);
            a = round(0.001/bin_size)+max_lag;
            b = conv(A, ones(1, 11), 'valid');
            ratio(c1, c2) = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*22)/(min(b)*2 - min(A)*22);
            [h, filteredA] = find_smallest_h(A);
            [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
            p(c1, c2) = sum(bootstat > h)/N;
            subplot(5, 5, c2)
            if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1
                bar(xx, A, 'r')
                corr_cells_test = [corr_cells_test; id1 id2];
            else
                bar(xx, A, 'b')
            end
            title([num2str(id1) '  ' num2str(id2) '  ' num2str(p(c1, c2))])
            xlim([-0.01 0.01])
        end
    end
    c1
end
