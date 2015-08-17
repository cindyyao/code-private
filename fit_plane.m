function [normal, meanX, pcaExplained, sse] = fit_plane(data)

[coeff, score, roots] = pca(data);
normal = coeff(:,3);
pcaExplained = roots' ./ sum(roots);

[n,~] = size(data);
meanX = mean(data,1);
Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';
error = abs((data - repmat(meanX,n,1))*normal);
sse = sum(error.^2);

figure
[xgrid,ygrid] = meshgrid(linspace(min(data(:,1)),max(data(:,1)),5), ...
                         linspace(min(data(:,2)),max(data(:,2)),5));
zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);

hold on
above = (data-repmat(meanX,n,1))*normal < 0;
below = ~above;
nabove = sum(above);
X1 = [data(above,1) Xfit(above,1) nan*ones(nabove,1)];
X2 = [data(above,2) Xfit(above,2) nan*ones(nabove,1)];
X3 = [data(above,3) Xfit(above,3) nan*ones(nabove,1)];
plot3(X1',X2',X3','-', data(above,1),data(above,2),data(above,3),'o', 'Color',[0 .7 0]);
nbelow = sum(below);
X1 = [data(below,1) Xfit(below,1) nan*ones(nbelow,1)];
X2 = [data(below,2) Xfit(below,2) nan*ones(nbelow,1)];
X3 = [data(below,3) Xfit(below,3) nan*ones(nbelow,1)];
plot3(X1',X2',X3','-', data(below,1),data(below,2),data(below,3),'o', 'Color',[1 0 0]);
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')

end