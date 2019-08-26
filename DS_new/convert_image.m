imageNB = imread('/Users/xyao/Dropbox/DS_paper/figures/raw_figures/fig5_new/Untitled-4.jpg');
imageNB = imageNB(3:end, :, :);
figure
image(imageNB);
imageNB_new = ones(size(imageNB));
imageNB_new(:, :, 1) = imageNB(:, :, 2);
imageNB_new = uint8(imageNB_new)*1.2;
% imageNB_new(:, :, 2) = imageNB(:, :, 1);
% imageNB_new(:, :, 3) = imageNB(:, :, 3);
figure
image(imageNB_new);
imwrite(imageNB_new, 'NB', 'TIFF')
%%
imageAI9 = imread('/Users/xyao/Dropbox/DS_paper/figures/raw_figures/fig5_new/Untitled-5.jpg');
imageAI9 = imageAI9(3:end-1, 4:end-1, :);
figure
image(imageAI9);
imageAI9_new = ones(size(imageAI9));
imageAI9_new(:, :, 2) = imageAI9(:, :, 1);
imageAI9_new = uint8(imageAI9_new)*0.8;
% imageNB_new(:, :, 2) = imageNB(:, :, 1);
% imageNB_new(:, :, 3) = imageNB(:, :, 3);
figure
image(imageAI9_new);
imwrite(imageAI9_new, 'AI9', 'TIFF')

merge = ones(size(imageAI9));
merge(:, :, 1) = imageNB_new(:, :, 1);
merge(:, :, 2) = imageAI9_new(:, :, 2);
merge = uint8(merge);
figure
image(merge);
imwrite(merge, 'Merge', 'TIFF')
