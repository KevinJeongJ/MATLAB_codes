clear all
close all
load('EstSnap13.mat');
snapshotsV = GlobalSnap(1:size(GlobalSnap,1)/2,:);
snapshotsW = GlobalSnap(size(GlobalSnap,1)/2+1:size(GlobalSnap,1),:);

%% plot of snapshot to view actual flow

% for i=1:100:1000
% figure(1)
% scatter(visz, visy, [], snapshotsV(:,i), 'filled');
% title('snapshotsV');
% axis equal
% end

%% plot
k= 10;
visv = snapshotsV(:, k:k+179);

%make uniform grid
xdummy = 6.3352:1.8296:346.6375;
ydummy = -147.8552:1.8282:150.1393;


nl = length(xdummy)*length(ydummy);
xnew = zeros(1, nl);
ynew = zeros(1, nl);

count = 0;
for i = 1:length(xdummy)
    for j = 1:length(ydummy)
        count = count+1;
        xnew(count) = xdummy(i);
        ynew(count) = ydummy(j);
    end
end

interpolated_visv = zeros(nl, size(visv,2));
for i = 1:size(visv,2)
    interpolated_visv(:,i) = griddata(visz, visy, visv(:,i), xnew, ynew);
end
X = zeros(nl,size(visv,2));
for i= 1:size(visv,2)
    X(:,i)=griddata(visz,visy,visv(:,i),xnew,ynew,'nearest');
end
for i=1:length(ydummy)
    interpolated_visv(i,:)=X(i,:);
end


%prepare data into form for contour plot
p = length(xdummy);
n = length(ydummy);
xCon = reshape(xnew,[n,p]);
yCon = reshape(ynew,[n,p]);

%% DMD
B = interpolated_visv(:,1:end-1);
A = interpolated_visv(:,2:end);

[L,M] = omd(A,B,10,[],'alternating'); 
[Y,D_nu]=eig(M);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd

Vand=zeros(size(M));
for loop1=1:size(B,2)
    Vand(:,loop1)=diag(D_nu).^(loop1-1); %creating vandermonde with increasing powers
end
U=L;
Psi = U*Y;

%% alpha - OMD

Z(:,1) =(Y^-1)*U'*B(:,1);

for i = 1:size(Vand,1)
    for j = 1:size(Vand,1)
        alpha_omd(j,i) = Z(j,1)*Vand(j,i);
    end
end

%%
dt = 1/720;
freq = 0.1965/15/(2*pi)*imag(log(diag(D_nu))/dt);

%% filter
i = 1;
for level = 8
  
    figure(level)
    set(gcf,'pos',[10,10,1200,500]);
    %
    %
    DLamD = reshape(Psi(:,i),[n,p]);
    DLamD = real(DLamD);
    [c h] = contourf(xCon/200,yCon/200,DLamD,level);
%     clevel = h.LevelList;
%     DLamD(DLamD<clevel(2)) = 50;
%     DLamD(DLamD>clevel(3)) = 100;
%     DLamD = imcomplement(DLamD);
    hold all
    % normal image
    subplot(2,4,1)
    contourf(xCon/200,yCon/200,real(DLamD), level);
    title(['Real Mode'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

    
    % mediain filter
    img = medfilt2(real(DLamD));
    subplot(2,4,2)
    contourf(xCon/200,yCon/200,img,level);
    title(['Median'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

  
    % wiener 2 (2D adaptive noise-removal filtering
    subplot(2,4,3)
    imgw = wiener2(real(DLamD));
    contourf(xCon/200,yCon/200,imgw,level);
    title(['Wiener2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    
    %average
    subplot(2,4,4)
    H = fspecial('average',10);
    img2 = imfilter(real(DLamD),H,'replicate');
    contourf(xCon/200,yCon/200,img2,level);
    title(['Average'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

    
    %circular average
    subplot(2,4,5)
    H = fspecial('disk',10);
    imgd = imfilter(real(DLamD),H,'replicate');
    contourf(xCon/200,yCon/200,imgd,level);
    title(['Circ-Avg'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
  
    %guassian
    subplot(2,4,6)
    H = fspecial('gaussian',10);
    imgg = imfilter(real(DLamD),H,'replicate');
    contourf(xCon/200,yCon/200,imgg,level);
    title(['Guassian'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
   
    %laplacian
    subplot(2,4,7)
    H = fspecial('laplacian',1);
    imgl = imfilter(real(DLamD),H,'replicate');
    contourf(xCon/200,yCon/200,imgl,level);
    title(['Laplacian'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

    %log
    subplot(2,4,8)
    H = fspecial('log',10);
    imglo = imfilter(real(DLamD),H,'replicate');
    contourf(xCon/200,yCon/200,imglo,level);
    title(['Laplacian-Guassian'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

    
    
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray 
    
%     saveas(figure(level), fullfile('/home/kevin/Desktop//', ['figure' num2str(level) '.jpeg']));

end
