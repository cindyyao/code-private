
%stim
[mov,~,~,dur,refresh] = get_movie('/Volumes/G-DRIVE mobile USB-C/Analysis/BW-15-1-0.48-11111-53x40-60.35.xml',datarun.triggers,216000);
movie_white_noise = squeeze(mov(:,:,1,:));
slen = 215000;
mov1 = movie_white_noise(:,:,1:slen);
Stim0 = reshape(mov1,nky*nkx,[])'; %nky = size(movie_white_noise,1) nkx = size(movie_white_noise,2)
Stim = bsxfun(@minus,Stim0,mean(Stim0)); %subtract mean

%spikes
spikes = dataruns.(wn_runs_names{j}).spikes{ind};
nframes = 216000;
upsamplefactor = 10;
[binned_spikes_coarse,binned_spikes_fine,dtStim,~] = bin_spikes( spikes,dataruns.(wn_runs_names{j}).triggers,nframes,upsamplefactor); %just bin up spikes according to frame rate
rlen = upsamplefactor*slen; 
sps = binned_spikes_fine(1:rlen); 
%rebin coarsely 
sps_coarse = sum(reshape(sps,[],slen),1)';


%get sta
nkt = 30; %number of frames to go back in time
sta = simpleSTC(Stim,sps_coarse,nkt);









