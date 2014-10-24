spec.x_start = 0; spec.x_end = 800;
spec.y_start = 0; spec.y_end = 600;
spec.orientation = 0;
spec.spatial_period = 300;
spec.temporal_period = 480;
[framesin, framecos, tdrift] = calc_drifting_grating_frame_intensities_jf(spec);
figure
colormap gray
% set(gca,'YDir','normal')
im = framesin*cos(0)+framecos*sin(0);
% im = im(:, end:-1:1);
imagesc(im)
pause
im = framesin*cos(0.1)+framecos*sin(0.1);
% im = im(:, end:-1:1);
imagesc(im)
% set(gca,'YDir','normal')