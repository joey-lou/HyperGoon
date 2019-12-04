Z = dlmread('../496_New_Data/Saturation/saturaton1.txt');

%% POD
fprintf('Analyzing POD modes...\n')
[PSI,S,V] = svd(Z','econ');
% Compute DMD (Phi are eigenvectors)
r = length(find(cumsum(diag(S))/sum(diag(S))<=0.9999)); % dominant modes
fprintf('POD trucated at r = %i\n', r);
[X,Y] = meshgrid(1:60,1:220);
for k = 1:min(4,r)
    figure(4)
    subplot(2,2,k)
    contourf(X,Y,reshape(PSI(:,k),60,220)')
    colorbar
    title(['U POD mode ',num2str(k)])
%     figure(5)
%     subplot(3,3,k)
%     Grid.cplot(PSI(N+1:end,k))
%     title(['V POD mode ',num2str(k)])
end

x = normrnd(0,1,100,1);
y = normrnd(1,2,100,1);
scatter(x,y)
z = [x,y];
[PSI,S,V] = svd(z,'econ');

