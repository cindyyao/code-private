for i = 1:32
    figure
    filter = squeeze(norm_image(layer1(:, :, :, i)));
    for j = 1:20
        imshow(filter(:, :, j), 'InitialMagnification', 'fit')
%         colormap gray
        pause
    end
end

bias_dense = [network{7}; network{8}];
[~, scores, latent] = pca(bias_dense');
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')

figure
plot(latent)
