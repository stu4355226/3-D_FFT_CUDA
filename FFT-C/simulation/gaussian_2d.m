% Domain of our random variable 
% x = [-10:2:-3]
% y = [-3:0.6:3]
% z = [3:2:10]
% x = [x, y , z]
% y = x
x = [-10:.1:10]; 
y = [-10:.1:10]; 
% Grid domain 
[xx yy] = meshgrid(x,y); 
% Define mean and covariance matrix for Gaussian 
mu = [0 0]; 
covMat = [2 0;0 2]; 
% Vectorize inputs for fast calculation
jointRV = [xx(:) yy(:)]; 
muRep = repmat(mu,size(jointRV,1),1); 
% 2D Gaussian pdf 
gaussianPdf = (1/(sqrt(2*pi*det(covMat)))).*... 
    exp(-.5*(sum((jointRV-muRep)*inv(covMat).*(jointRV-muRep),2))); 
% Reshape to be same as our domain grid 
G = reshape(gaussianPdf,size(xx)); 
% 3D mesh plot 
figure(1); 
%surf(xx,yy,G); 
xx = reshape(xx, size(xx, 1) ^ 2, 1);
yy = reshape(yy, size(yy, 1) ^ 2, 1);
scatter3(xx, yy, gaussianPdf);
%colormap('jet'); 
%shading interp; 
title(['2 Dimensional Gaussian: \mu_x = ' num2str(mu(1)) ', \mu_x = ' num2str(mu(2))... 
    ', \sigma_x^2 = ' num2str(covMat(1,1)) ', \sigma_y^2 = ' num2str(covMat(2,2))]); 
% % Contour projection plot 
% figure(2); 
% [c hc] = contour(xx,yy,G,10,'r-'); 
% title(['2 Dimensional Gaussian Contour: \mu_x = ' num2str(mu(1)) ', \mu_x = ' num2str(mu(2))... 
%     ', \sigma_x^2 = ' num2str(covMat(1,1)) ', \sigma_y^2 = ' num2str(covMat(2,2))]);

dlmwrite('dataset',[xx, yy, gaussianPdf]);
% dataset = dlmread('dataset');
