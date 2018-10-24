clear all
close all
load('PIV13.mat');
VV3= reshape(PIVdata.vz, size(PIVdata.vz,1)*size(PIVdata.vz,2),size(PIVdata.vz,3));
snapshotsV = VV3;
visz = repelem(dim.z,size(dim.y,2));
visy = repmat(dim.y,size(dim.z));

% snapshotsV = GlobalSnap(1:size(GlobalSnap,1)/2,:);
% snapshotsW = GlobalSnap(size(GlobalSnap,1)/2+1:size(GlobalSnap,1),:);
%%
addpath('/home/kevin/Donwloads/PlotPub-master/lib');

%% plot of snapshot to view actual flow
% fig=figure('visible','off');
% set(gcf,'pos',[10 10 2500 500]);
% subplot(1,5,1);
% scatter(visz,visy,[],snapshotsV(:,1),'filled');
% colormap(gray)
% title(['SnapshotsV: 1'],'interpreter','latex','FontSize',20);
% axis equaler transform

% % saveas(figure(3), fullfile('/home/kevin/Desktop/Data/bullet//', ['figure' num2str(3) '.jpeg']));

% ylim([-150 150]);
% 
% hold aller transform

% % saveas(figure(3), fullfile('/home/kevin/Desktop/Data/bullet//', ['figure' num2str(3) '.jpeg']));

% k=2;
% for i=11:10:41
%     subplot(1,5,k);
%     scatter(visz, visy, [], snapshotsV(:,i), 'filled');
%     title([num2str(i)],'interpreter','latex','FontSize',20);
%     axis equal
%     ylim([-150 150]);
%     colormap(gray);
%     k= k+1;
% end
% 
%     saveas(fig, fullfile('/home/kevin/Desktop//', ['flow.jpeg']));


%% plot
k= 10;
visv = snapshotsV(:, k:k+500);

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
for i=1:size(length(ydummy))
    interpolated_visv(i,:)=X(i,:);
end


%prepare data into form for contour plot
p = length(xdummy);
n = length(ydummy);
xCon = reshape(xnew,[n,p]);
yCon = reshape(ynew,[n,p]);

%%
% fig = figure('visible','off');
% set(gcf,'pos',[10 10 2500 500]);
% subplot(1,5,1);
% scatter(xnew/200, ynew/200, [], interpolated_visv(:,1), 'filled');
% 
% axis equal
% axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
% colormap gray
%     title(['smooth: 1'],'interpreter','latex','FontSize',20);
% hold all
% 
% k=2;
% for i=11:10:41
%     subplot(1,5,k);
% scatter(xnew/200, ynew/200, [], interpolated_visv(:,i), 'filled');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     colormap gray
%         title([num2str(i)],'interpreter','latex','FontSize',20);
% 
%     k=k+1;
% end
%     saveas(fig, fullfile('/home/kevin/Desktop//', ['flow_smooth.jpeg']));


%%
% fig=figure('visible','on');
% set(gcf,'pos',[10 10 1200 200]);
% subplot(1,5,1);
% scatter(visz,visy,[],snapshotsV(:,11),'filled');
% colormap(gray)
% axis equal
% ylim([-150 150]);
% 
% hold all
% k=2;
% for i=11:10:41
%     subplot(1,5,k);
%     scatter(visz, visy, [], snapshotsV(:,i+10), 'filled');
%     axis equal
%     ylim([-150 150]);
%     colormap(gray);
%     k= k+1;
% end
% 
% %% all flow pic
% 
% fig=figure('visible','on');
% set(gcf,'pos',[10 10 1200 600]);
% subplot(3,5,1);
% scatter(visz,visy,[],snapshotsV(:,11),'filled');
% colormap(gray)
% title(['a)'],'interpreter','latex','FontSize',15);
% axis equal
% ylim([-150 150]);
% 
% hold all
% k=2;
% for i=11:10:41
%     subplot(3,5,k);
%     scatter(visz, visy, [], snapshotsV(:,i+10), 'filled');
%     axis equal
%     ylim([-150 150]);
%     colormap(gray);
%     k= k+1;
% end
% 
% subplot(3,5,6);
% dLamC = reshape(interpolated_visv(:,1), [n,p]); %reshape inter_visu to [n,p]
% contourf(xCon/200, yCon/200, dLamC, 8);
% axis equal
% axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
% colormap gray
% title(['b)'],'interpreter','latex','FontSize',15);
% hold all
% 
% k=2;
% for i=11:10:41
%     subplot(3,5,k+5);
%     dLamC = reshape(interpolated_visv(:,i), [n,p]); %reshape inter_visu to [n,p]
%     contourf(xCon/200, yCon/200, dLamC, 8);
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     colormap gray
%     
%     k=k+1;
% end
% 
% subplot(3,5,11);
% scatter(xnew/200, ynew/200, [], interpolated_visv(:,1), 'filled');
% axis equal
% axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
% title(['c)'],'interpreter','latex','FontSize',15);
% hold all
% 
% k=2;
% for i=11:10:41
%     subplot(3,5,10+k);
%     scatter(xnew/200, ynew/200, [], interpolated_visv(:,i), 'filled');
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     colormap gray    
%     k=k+1;
% end

%% DMD
B = interpolated_visv(:,1:end-1);
A = interpolated_visv(:,2:end);

[U,S,V] = svd(B,'econ');
F_dmd = U'*A*V/S;
[Y,D_nu]=eig(F_dmd);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd

Vand=zeros(size(F_dmd));
for loop1=1:size(B,2)
    Vand(:,loop1)=diag(D_nu).^(loop1-1); %creating vandermonde with increasing powers
end

Psi = U*Y;
%% omd 
[L,M] = omd(A,B,179,[],'alternating'); 
[Q,D_omd]=eig(M);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd
Psi_omd = L*Q;

Vand_o=zeros(size(M));
for loop1=1:size(B,2)
    Vand_o(:,loop1)=diag(D_omd).^(loop1-1); %creating vandermonde with increasing powers
end

Z(:,1) =(Q^-1)*L'*B(:,1);

for i = 1:size(Vand_o,1)
    for j = 1:size(Vand_o,1)
        alpha_omd(j,i) = Z(j,1)*Vand_o(j,i);
    end
end

%%


%% alpha - POD
alpha = U'*B;
% 
% for i=1:20 %size(alpha,2)
%     figure
%     plot(alpha(i,:));
% end

%mode energy plot
%%
% fig=figure;
% set(gcf,'pos',[10 10 2500 500]);
% bar(diag(M(1:60,1:60)/norm(M)));
% set(gca,'Yscale','log');
% title(['OMD'],'interpreter','latex','FontSize',22);
% ylabel(['normalized mode energy'],'interpreter','latex','fontsize',22);
% xlabel(['mode'],'interpreter','latex','fontsize',22);

%% alpha - DMD

Z(:,1) =(Y^-1)*U'*B(:,1);

for i = 1:size(Vand,1)
    for j = 1:size(Vand,1)
        alpha_dmd(j,i) = Z(j,1)*Vand(j,i);
    end
end

[U_dmd,S_dmd,V_dmd]=svd(U'*B,'econ'); %singular value of U'B, like POD with B

%% plot - POD
fig = figure;
set(gcf,'pos',[10 10 1000/3*5 250]);
subplot(1,5,1);
DLamD = reshape(U(:,1),[n,p]);
contourf(xCon/200,yCon/200,real(DLamD),8);
axis equal
axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
colormap gray
title(['POD: 1'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


k=2;
for i =2:5
    DLamD = reshape(U(:,i),[n,p]);
    subplot(1,5,k)
    contourf(xCon/200,yCon/200,real(DLamD),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
    title([num2str(i)],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


    k=k+1;
%     saveas(figure(i), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_v/bullet_v_pod//', ['figure' num2str(i) '.jpeg']));
end
%{
%% if r</< n DMD

phi_w = A*V*(S^-1)*Y;

% fig =figure;
% set(gcf,'pos',[10 10 2500 500]);
% subplot(1,5,1);
% DLamD = reshape(Psi(:,1),[n,p]);
% contourf(xCon/200,yCon/200,real(DLamD),8);
% title(['DMD: #1'],'interpreter','latex','FontSize',20);
% 
% hold all
% 
% for i = 2:5 %size(Psi,2) %number of modes
%     DLamD = reshape(Psi(:,i),[n,p]);
%     subplot(1,5,i)
%     contourf(xCon/200,yCon/200,real(DLamD), 8);
%     title(['#',num2str(i)]);
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     colormap gray
% end

%% plot - via maximum alpha_dmd
% 
alpha_dmd_transpose = alpha_dmd.';
% row and columns of alpha_dmd are switched
alpha_dmd_abs = abs(real(alpha_dmd_transpose));

[M,J] = max(alpha_dmd_abs);

[H, J] = sort(M,'descend');
% 
% for i=1:30 %size(alpha_dmd,2)
%     k = I(i);
%     DLamD = reshape(Psi(:,k),[n,p]);
%     figure('name','DMD','pos',[10,10,600,900]);
%     subplot(3,3,[1,2,3,4,5,6])
%     contourf(xCon/200,yCon/200,real(DLamD), 8);
%     title(['Real: Mode',num2str(k-1)]);
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     colormap gray
%     
% %         subplot(2,1,2)
% %         contourf(xCon/200,yCon/200,imag(DLamC), 8);
% %         title(['Imag: Mode',num2str(i-1)]);
% %         axis equal
% %         axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
% %         colormap gray
% %     
%     subplot(3,3,[7,8,9])
%     plot(real(alpha_dmd(k,:)));
%     xlim([0, inf]);
%     
%     saveas(figure(i), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_v/bullet_v_dmd_alpha//', ['figure' num2str(i) '.jpeg']));
% end

%% I method DMD
% 
alpha_dmd_abs_I = abs(real(alpha_dmd));
alpha_dmd_abs_I_t = alpha_dmd_abs_I.';

dt = 1/720;
freq = 0.1965/15/(2*pi)*imag(log(diag(D_nu))/dt);


I = sum(alpha_dmd_abs_I_t);
[C,T] = sort(I,'descend');


    
fig = figure;
set(gcf,'pos',[10 10 1000 500]);

subplot(2,3,1);
DLamD = reshape(Psi(:,T(1)),[n,p]);
contourf(xCon/200,yCon/200,real(DLamD),8);
axis equal
axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
colormap gray
title(['I-DMD Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

subplot(2,3,4);
contourf(xCon/200,yCon/200,imag(DLamD),8);
axis equal
axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
colormap gray
title(['Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


k=2;
for i =3
    j = T(i);
    DLamD = reshape(Psi(:,j),[n,p]);
    subplot(2,3,k)
    contourf(xCon/200,yCon/200,real(DLamD),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
    title(['St =',num2str(abs(freq(j)))],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

    subplot(2,3,3+k)
    contourf(xCon/200,yCon/200,imag(DLamD),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray


    k=k+1;
%     saveas(figure(i), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_v/bullet_v_pod//', ['figure' num2str(i) '.jpeg']));
end
for i =5
    j = T(i);
    DLamD = reshape(Psi(:,j),[n,p]);
    subplot(2,3,k)
    contourf(xCon/200,yCon/200,real(DLamD),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
    title(['St =',num2str(abs(freq(j)))],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

    subplot(2,3,3+k)
    contourf(xCon/200,yCon/200,imag(DLamD),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray


    k=k+1;
%     saveas(figure(i), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_v/bullet_v_pod//', ['figure' num2str(i) '.jpeg']));
end

%% misc
k = 60;
% %comparing S_dmd to S
figure
set(gcf,'pos',[10 10 1000 500]);


bar(diag(S(1:k,1:k)));
title(['POD: Energy'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
set(gca,'Yscale','log');
ylabel(['Normalized Mode Energy'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
xlabel(['Mode'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

% 
% figure
% bar(diag(S_dmd(1:k,1:k)));
% title('energy mode DMD');
% set(gca,'Yscale','log');
% 
% % % comparing alpha_dmd against alpha_pod
% % for i = 1:k
% %     figure
% %     subplot(2,1,1)
% %     plot (alpha(i,:));
% %     title('\alpha POD');
% %     subplot(2,1,2)
% %     plot(real(alpha_dmd(i,:)));
% %     title('\alpha DMD');
% % end
% 

% subplot(1,3,2);
% plot(D,'o');
% circle implementation
%     hold on
%     th=0:pi/50:2*pi;
%     plot(cos(th),sin(th));
%     hold off
% xlim([0,1.5]);
% ylim([-1.5,1.5]);
% xlabel('Re(\lambda_i)');
% ylabel('Im(\lambda_i)');
% title('DMD eigenvalues)');



%% freq
dt = 1/720;
freq = 0.1965/15/(2*pi)*imag(log(diag(D_nu))/dt);

for i = 1:size(alpha_dmd,2)
    alpha_dmd_norm(i) = norm(real(alpha_dmd(i,:)));
end
alpha_dmd_norm = alpha_dmd_norm/norm(alpha_dmd_norm);

for i = 1:size(alpha,2)
    alpha_pod_norm(i) = norm(alpha(i,:));
end
t = alpha_pod_norm/norm(alpha_pod_norm);

figure
set(gcf,'pos',[10 10 1500 400]);
subplot(1,2,2);
plot(freq,I,'o');
title(['I-method'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
xlabel(['St'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
ylabel(['||\alpha_{DMD}||'],'fontweight','normal','FontSize',22,'fontname','arial');
set(gca,'Yscale','log');
ylim([0, 5*10^4]);

subplot(1,2,1);
plot(freq,J,'o');
title(['||\alpha_{max}||'],'fontweight','normal','FontSize',22,'fontname','arial');
xlabel(['St'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
ylabel(['||\alpha_{DMD}||'],'fontweight','normal','FontSize',22,'fontname','arial');

ylim([0, 5*10^4]);
set(gca,'Yscale','log');

%%

% subplot(1,3,2)
% for i=1:178
% Kst(i) = max(abs((real(alpha_omd(i,:)))));
% freq_2(i) = freq(i);
% end
% 
% plot(freq_2,Kst,'o');
% set(gca,'Yscale','log')
% title(['\alpha_{OMD}'],'fontsize',13);
% xlabel(['St'],'interpreter','latex','fontsize',15);
% ylabel('||\alpha_{OMD}||');
% ylim([0, 5*10^4]);

% subplot(1,3,3)
% 
clear Ts alpha_omd_5
[L,M] = omd(A,B,5,[],'alternating'); 
[Q,D_omd]=eig(M);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd
Psi_omd = L*Q;

Vand_o=zeros(size(M));
for loop1=1:size(B,2)
    Vand_o(:,loop1)=diag(D_omd).^(loop1-1); %creating vandermonde with increasing powers
end

Ts(:,1) =(Q^-1)*L'*B(:,1);

for i = 1:size(Vand_o,1)
    for j = 1:size(Vand_o,1)
        alpha_omd_5(j,i) = Ts(j,1)*Vand_o(j,i);
    end
end
freq_omd = 0.1965/15/(2*pi)*imag(log(diag(D_omd))/dt);
% for i= 1:5
% Ks(i) = abs(max((real(alpha_omd_5(i,:)))));
% end
% plot(freq_omd,Ks,'o');
% set(gca,'Yscale','log')
% title(['\alpha_{OMD}: 10 modes'],'fontsize',13);
% xlabel(['St'],'interpreter','latex','fontsize',15);
% ylabel('||\alpha_{OMD}||');
% ylim([0, 5*10^4]);



% saveas(figure(1), fullfile('/home/kevin/Desktop/Data/bullet//', ['figure' num2str(1) '.jpeg']));
% 
% 
% figure (2)
% plot(alpha_dmd_norm);
% xlabel('n');
% ylabel('||\phi_{DMD}||');
% saveas(figure(2), fullfile('/home/kevin/Desktop/Data/bullet//', ['figure' num2str(2) '.jpeg']));

%%
fig = figure;
set(gcf,'pos',[10 10 1000 500]);
subplot(2,3,1);
DLamE = reshape(Psi_omd(:,3),[n,p]);
contourf(xCon/200,yCon/200,real(DLamE),8);
axis equal
axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
colormap gray
title(['OMD Real: 1'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');

subplot(2,3,4);
contourf(xCon/200,yCon/200,imag(DLamE),8);
axis equal
axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
colormap gray
title(['Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


hold all
k=2;
for i =2
    DLamE = reshape(Psi_omd(:,i),[n,p]);
    subplot(2,3,k)
    contourf(xCon/200,yCon/200,real(DLamE),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
    title(['St =',num2str(abs(freq_omd(i)))],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    subplot(2,3,k+3);
    contourf(xCon/200,yCon/200,imag(DLamE),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
end

k=3;
for i =4
    DLamE = reshape(Psi_omd(:,i),[n,p]);
    subplot(2,3,k)
    contourf(xCon/200,yCon/200,real(DLamE),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
    title(['St =',num2str(abs(freq_omd(i)))],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    
    subplot(2,3,k+3);
    contourf(xCon/200,yCon/200,imag(DLamE),8);
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    colormap gray
end


%     saveas(figure(i), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_v/bullet_v_pod//', ['figure' num2str(i) '.jpeg']));


%% fourier transform
fig = figure;
set(gcf,'pos',[10 10 1000 300]);
for c = 2:5
    a = V(c,:);
%  a = (((Psi_omd(:,T(c)).')*Psi_omd(:,T(c))).^-1)*Psi_omd(:,T(c)).'*B;


Fs = 720*0.1965/15; %sampling frequency
period = 0:1/Fs:Fs*500;
N = length(a);
xdft = fft(a);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N))*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq_ = 0:Fs/length(a):Fs/2;

hold all
plot(freq_,psdx);
set(gca,'Xscale','log');
end

xlim([-inf,0.5]);
legend('mode 2','mode 3', 'mode 4', 'mode 5');
xlabel(['St'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
ylabel(['V spectra'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
box on
xticks([0.04 0.06 0.2 0.4]);
% % % saveas(figure(3), fullfile('/home/kevin/Desktop/Data/bullet//', ['figure' num2str(3) '.jpeg']));
% 
% set(gca,'Xscale','log');
% xlim([-inf,0.5]);
% end
% legend('mode 2','mode 3', 'mode 4', 'mode 5');
% xlabel(['St'],'interpreter','latex','FontSize',15);
% ylabel(['V spectra'],'interpreter','latex','FontSize',15);

% % saveas(figure(3), fullfile('/home/kevin/Desktop/Data/bullet//', ['figure' num2str(3) '.jpeg']));

%}
%% omd function
function [L,M] = omd(A,B,k,L0,method,userOpts)

%OMD : solves an optimal low rank approximation problem 
% in the form:
%
% min \|A - XB\|^2
% s.t. X = LML', 
%
% where L is an orthonormal matrix with k columns, and M is 
% a square matrix of size k x k.
%
% Usage : [L,M] = omd(A,B,k,L0,method)
%
% Inputs:
%
% A,B    : data matrices of size p x n
% k      : integer specifying the desired output size.  
% L0     : initial condition for L.  If L0 is not specified,
%          then the initial iterate will be based on the 
%          first k singular vectors of the matrix [A B]
%
% method : Specifies which optimization method to use.
%          It must be one of the following:
%          
%          'alternating' : Uses an approximate method of
%                          alternating directions.  Good for
%                          obtaining a rough solution quickly
%
%          'gradient'    : Uses the gradient descent method
%
%          'conjgrad'    : Uses the conjugate gradient method
%
%          'hybrid'      : Uses the alternating method first
%                          to get a rough solution, followed
%                          by the conjugate gradient method
%                          to achieve local optimality.  This
%                          is usually the best way to get a 
%                          high-precision solution quickly.
%
%          'dmd'         : Returns the dynamic mode decomposition
%
%          Note that if the 'dmd' method is chosen, then no
%          optimization is performed over the input basis L0. 
%          If L0 is not specified, then function returns the 
%          standard dynamic mode decomposition.
%                   
%
% Outputs:
%
% L      : optimal solution basis for the problem above,  
%          with L \in R^{p x k} 
% M      : optimal solution matrix for the problem above,
%          with M \in R^{k x k}
%
% The input data matrices (A,B) must both be size
% p x n, and full column rank, with k \le n \le p
%
%  

% Author : P. Goulart - 20 June 2012
%
% This function is an implementation of the algorithms described
% in the following publications: 
%
% Goulart, Wynn & Pearson, 'Optimal mode decomposition for 
% high-dimensional systems. In 51st IEEE Conference on Decision 
% and Control. Maui, Hawaii, Dec. 2012.  Available at 
% http:\\control.ee.ethz.ch\~goularpa\.
%
% Wynn, Pearson, Ganapathisubramani & Goulart, 'Optimal mode 
% decomposition for unsteady and turbulent flows', June 2012.
% Submitted to Journal of Fluid Mechanics. Available at 
% http:\\control.ee.ethz.ch\~goularpa\  
%
% Edelman et. al. 'The geometry of algorithms with orthogonality 
% constraints' SIAM J. Matrix Anal. Appl. 20(2) 303-353).


%set the default solver options, and merge
%in the user ones
opts.relTol  = 1e-4;
opts.maxIter = 100;
opts.maxStep = 0.1;
opts.lineSearch = [];  %leaves line search defaults to linesearch function
if(nargin == 6)
    opts = setstructfields(opts,userOpts);
end


%Check that a solver was specified
if(nargin < 5 || isempty(method))
    method = 'alternating';
end

%Find an initial L if needed
if(nargin < 4 || isempty(L0)) 
    
    %Use the initial bases of the POD modes
    %Note that in some problems svd([A B],0)
    %might actually be better, but we use
    %B only for consistency with the DMD method
	[L0,~] = svd(B,0);
    L0 = L0(:,1:k);
else
    L0 = orth(L0); %use what was provided
end



switch method
    
    case 'hybrid'
    %----------------- 
%         fprintf('\n\nInitiating pre-solve step\n');
        L = lmlOpt_altProject(A,B,L0,opts);
%         fprintf('\n\nInitiating refinement step\n');
        [L,M] = omd(A,B,k,L,'conjgrad',opts);
        return;
    
    
    case 'alternating'
    %-----------------    
                
        %call the solver;
        L = lmlOpt_altProject(A,B,L0,opts);
        
    case 'conjgrad'
    %-----------------    
    
        %this method does not work in the special
        %case of 1-D data and k = 1.  In that case,
        %call the gradient method
        if(k == 1 && size(A,1) == 1)
            fprintf('\nWarning: 1-d problem.  Calling gradient method');
            [L,M] = omd(A,B,1,L0,'gradient');
            return;
        end
                
        %objective function and its gradient
        nA = norm(A,'fro')^2;
        g  = @(L)(sqrt(nA-get_g(L,A,B)^2));
        dg = @(L)(-get_dgdL(L,A,B));
       
        %call the solver
        L = grassOpt_cg(g,dg,L0,opts);
               
        
    case 'gradient'
    %----------------- 
                
        %objective function and its gradient
        nA = norm(A,'fro')^2;
        g  = @(L)(sqrt(nA-get_g(L,A,B)^2));
        dg = @(L)(-get_dgdL(L,A,B));
        
        %call the solver
        L = grassOpt_gradient(g,dg,L0,opts);    
        
        
    case 'dmd'
    %-----------------
        
        %calculates the 'dynamic mode decomposition',
        %i.e. fix L to the first k vectors in L0 and
        %then find the corresponding M
        L = L0(:,1:k);
        
    otherwise
    %-----------------
        error('Unknown solver method')
end


%Compute the optimal M and return;
M = getMfromL(A,B,L);

end


    

%----------------------------------------------------
%----------------------------------------------------

function n = get_g(L,A,B)

%computes the function g(L), where
%g(L) = \|L'*A*C\|^2, where C = orth(B'*L)

C = orth(B'*L);
n = norm((L'*A)*C,'fro');
end



%----------------------------------------------------
%----------------------------------------------------

function dgdL = get_dgdL(L,A,B)

%computes the gradient of function g(L), where
%g(L) = \|L'*A*C\|^2, where C = orth(B'*L)

BtL = B'*L;
AtL = A'*L;
U = AtL'*BtL;
V = BtL'*BtL;
iVU = V\(U');

%Find the derivative 
dgdL = + 2*(A*(BtL*iVU) + B*(AtL*iVU'))  ...
       - 2*B*(BtL*iVU*iVU');       
end
%----------------------------------------------------
%----------------------------------------------------

function M = getMfromL(A,B,L)

% getMfromL : Given input data (A,B,L), computes an 
% optimal M as
%
%   M = (L'AB'L)(L'BB'L)^{-1}
%
% Usage: M = getMfromL(A,B,L)


LtB = L'*B;
LtA = L'*A;
M   = (LtA * LtB') / (LtB * LtB');
end

%
%
%

%% GrassOpt_gradient
function L = grassOpt_gradient(g,dg,L0,opts)

% grassOpt_gradient : solves an optimization problem on the 
% Grassman manifold using a gradient method.
%
% This function solves the problem:
%
% min g(L)
% s.t. L'L = I (not necessarily square)
%
% where the function g() is invariant with respect to
% an orthonormal transformation of L, i.e. when 
% g(L) = g(LR) for any R'R = I, R square.
%
% Usage : L = grassOpt_cg(g,dg,L0,opts)
%
% Inputs:
%
% g     : function handle for evaluating g(L)
% dg    : function handle for evaluating \nabla g(L)  
% L0    : initial condition for L.  This must be
%         a matrix satisfying the constraint L'L=I.
% opts  : options structure with the following fields:
%         'maxIter', 'relTol' and 'maxStep'.
%
% Outputs:
%
% L     : solution for the above problem
%
% This function is an implementation of the algorithm described
% in the following publications: 
%
% Edelman et. al. 'The geometry of algorithms with orthogonality 
% constraints' SIAM J. Matrix Anal. Appl. 20(2) 303-353).

% Author : P. Goulart - 20 June 2012
%


%intialize convergence test variables
iterCount = 0;
relError  = 1;
valueOld  = g(L0); 

%initialize algorithm parameters
L     = L0;
Llast = L0;

%print the reporting header
% printInfo();

%get an initial step-size
tstep = opts.maxStep;

while( abs(relError) > opts.relTol && iterCount < opts.maxIter)
    
    %update counter
    iterCount = iterCount + 1;
            
    %find the component of the gradient tangent to 
    %the Grassman manifold
    %following the method of \S2.5.3 in Edelman
    dgdL   = dg(L);
    G      = dgdL - L*(L'*dgdL);
    
    %search direction is negative gradient
    D = -G;
    
    %compact SVD decomp for D
    [U,S,V] = svd(D,0);
    
    %do a line search to minimize g(t).  
    k = pi/max(diag(S));  %normalization for geodesic search
    fline = @(t)(g(L*V*dcos(S*(k*t))*V' + U*dsin(S*(k*t))*V'));
    [tstep,valueNew] = linesearch(fline,0,tstep,opts.lineSearch,valueOld);
    
    %update L along a geodesic
    L = L*V*dcos(S.*(k*tstep))*V' + U*dsin(S.*(k*tstep))*V';
    
    %renormalize, just in case...
    L = orth(L);
            
    %calculate relative improvement
    relError = (valueOld-valueNew)/abs(valueOld);
    valueOld = valueNew; %for next pass
    
    %print reporting information
%     printInfo(iterCount,valueNew,relError);
        
    if(relError <= 0 || tstep == 0)
        %not improving.  Force stop.
        relError = 0;
        L = Llast; %previous iterate
    else
        Llast = L;
    end

end

if(iterCount == opts.maxIter)
    warning('Didn''t converge')
end
end


%----------------------------------------------------
%----------------------------------------------------

% function printInfo(iterCount,val,relerr)
% 
% %print reporting information
% 
% if(nargin < 1)
%     %print the header
%     fprintf(1,'\n\nStarting solver: %s\n\n',mfilename);
%     fprintf(1,'iterate  |  objective value  |  relative improvement\n');
%     fprintf(1,'----------------------------------------------------\n');
% else
%     fprintf(1,'%4i     |  %0.6e    |  %0.6e \n',iterCount,val,relerr);
% end
% end

   
%----------------------------------------------------
%----------------------------------------------------

function D = dsin(S)

%Computes a diagonal matrix by taking 
%sines of the diagonal elements of the
%input.  This is 'sin' in the sense of 
%Edelman (2.65)

D = diag(sin(diag(S)));

end



%----------------------------------------------------
%----------------------------------------------------

function D = dcos(S)

%Computes a diagonal matrix by taking 
%sines of the diagonal elements of the
%input.  This is 'cos' in the sense of 
%Edelman (2.65)

D = diag(cos(diag(S)));
end

%
%
%

%% grassOpt_cg
function L = grassOpt_cg(g,dg,L0,opts)

% grassOpt_cg : solves an optimization problem on the 
% Grassman manifold using a conjugate gradient method.
%
% This function solves the problem:
%
% min g(L)
% s.t. L'L = I (not necessarily square)
%
% where the function g() is invariant with respect to
% an orthonormal transformation of L, i.e. when 
% g(L) = g(LR) for any R'R = I, R square.
%
% Usage : L = grassOpt_cg(g,dg,L0,opts)
%
% Inputs:
%
% g     : function handle for evaluating g(L)
% dg    : function handle for evaluating \nabla g(L)  
% L0    : initial condition for L.  This must be
%         a matrix satisfying the constraint L'L=I.
% opts  : options structure with the following fields:
%         'maxIter', 'relTol' and 'maxStep'.
%
% Outputs:
%
% L     : solution for the above problem
%
% This function is an implementation of the algorithm described
% in the following publications: 
%
% Edelman et. al. 'The geometry of algorithms with orthogonality 
% constraints' SIAM J. Matrix Anal. Appl. 20(2) 303-353).

% Author : P. Goulart - 20 June 2012
%


%intialize convergence test variables
iterCount = 0;
relError  = 1;
valueOld  = g(L0); 

%initialize algorithm parameters
Yk = L0;
dgdY = dg(Yk); 
Gk   = dgdY - Yk*(Yk'*dgdY);
Hk   = -Gk;  

%print the reporting header
% printInfo();

%get an initial step-size
tstep = opts.maxStep;

while( abs(relError) > opts.relTol && iterCount < opts.maxIter)
    
    %update counter
    iterCount = iterCount + 1;
            
    %minimize g(L_k(t)) over t
    [U,S,V] = svd(Hk,0);
    
    %do a line search to minimize g(t).  
    k = pi/max(diag(S));  %normalization for geodesic search
    fline = @(t)(g(Yk*V*dcos(S*(k*t))*V' + U*dsin(S*(k*t))*V'));
    [tstep,valueNew] = linesearch(fline,0,tstep,opts.lineSearch,valueOld);
   
    %update L along a geodesic
    Ykp1 = Yk*V*dcos(S.*(k*tstep))*V' + U*dsin(S.*(k*tstep))*V';
        
    %renormalize, just in case...
    Ykp1 = orth(Ykp1);
    
    %find gradient for new L.
    dgdY = dg(Ykp1);
    
    %update G
    Gkp1 = dgdY - Ykp1*(Ykp1'*dgdY);
    
    %parallel transport H and G
    tHk = (-Yk*V*dsin(S.*(tstep*k))+U*dcos(S.*(tstep*k)))*S*V';
    tGk = Gk - (Yk*V*dsin(S.*(tstep*k))+U*(eye(size(S)) - dcos(S.*(tstep*k))))*(U'*Gk);
    
    %new search direction
    normGk2 = norm(Gk,'fro')^2;
    gammak  = sum(sum(Gkp1.*(Gkp1-tGk))) ./ normGk2;
    Hkp1    = -Gkp1 + gammak.*tHk; 
            
    %calculate relative improvement
    relError = (valueOld-valueNew)/abs(valueOld);
    valueOld = valueNew; %for next pass
    
    %print reporting information
%     printInfo(iterCount,valueNew,relError);
    
    %If the search direction and gradient are too far
    %apart (arbitarily defined...), then restart
    if(trace(Gkp1'*Hkp1) / trace(Gkp1'*Gkp1) > -0.5)
        Hkp1 = -Gkp1;
        tstep = opts.maxStep;
        fprintf('Reinitializing conjugate gradient method\n');
    end
        
    if(relError <= 0 || tstep == 0)
        %not improving.  Force halt.
        relError = 0;
        L = Yk; %from the previous iteration
    else
        %update L, G, H iterates
        Yk = Ykp1;
        Gk = Gkp1;
        Hk = Hkp1;
        %update the output matrix
        L  = Yk;
    end
    
    %NB: p(n-p) can be really huge, so don't 
    %reset Hk as in Edelman paper.
    
end


if(iterCount == opts.maxIter)
    warning('Didn''t converge')
end


end






%% lmlOpt_alt

function L = lmlOpt_altProject(A,B,L0,opts)

% grassOpt_alternate : finds an approximate solution to the low rank
% matrix approximation problem in the form:
%
% min \|A - XB\|^2
% s.t. X = LML', 
%
% where L is an orthonormal matrix with k columns, and M is 
% a square matrix of size k x k.
%
% Usage : grassOpt_alternate(A,B,L0,opts)
%
% Inputs:
%
% A,B    : data matrices of size p x n
% L0     : initial condition for L.  
% opts   : options structure with the following fields:
%         'maxIter', 'relTol' and 'maxStep'.
%
% Outputs:
%
% L      : optimal solution basis for the problem above,  
%          with L \in R^{p x k} 
%
% The input data matrices (A,B) must both be size
% p x n, and full column rank, with k \le n \le p
%
%  

% Author : P. Goulart - 20 June 2012
%
% This function is an implementation of the alternating 
% projection method described in Algorithm 1 of:
%
% Goulart, Wynn & Pearson, 'Optimal mode decomposition for 
% high-dimensional systems. In 51st IEEE Conference on Decision 
% and Control. Maui, Hawaii, Dec. 2012.  Available at 
% http:\\control.ee.ethz.ch\~goularpa\.


%initialize convergence test variables
iterCount = 0;
relError  = 1;
valueOld  = realmax; 

%initialize algorithm parameters
L = L0;
Llast = L;
k = size(L,2);
%C = orth(A'*L);
C = eye(size(A,2),size(L,2));

%print the reporting header
% printInfo();

%norm squared of A.  Calculated so that
%reporting of values will be consistent
%with the overall objective function
nA = trace(A'*A);


while( abs(relError) > opts.relTol && iterCount < opts.maxIter)
    
    %update counter
    iterCount = iterCount + 1;
            
    %Find optimal L' ignoring constraints on C
    L = solve_inner((A*C)',k);
    
    %compute the image of B'*L
    C = orth(B'*L);
    
    %calculate relative improvement
    valueNew = sqrt(nA - getNorm(L,A,C)^2);
    relError = (valueOld - valueNew)/abs(valueOld);
    valueOld = valueNew; %for next pass
    
    %print reporting information
%     printInfo(iterCount,valueNew,relError);
    
    if(relError < 0)
        %going backwards.  Force halt.
        relError = 0;
        L = Llast; %from the previous iteration
    else
        Llast = L;
    end

end

if(iterCount == opts.maxIter)
    warning('Didn''t converge')
end
end

%----------------------------------------------------
%----------------------------------------------------

function X = solve_inner(D,k)

% solver for the unconstrained step
%

%Compute a basis for D'
V = orth(D');

%X coincides with the first k right singular vectors
X = V(:,1:k);

end


%----------------------------------------------------
%----------------------------------------------------

function n = getNorm(L,A,C)

n = norm((L'*A)*C,'fro');

end

%
%
%

%% linesearch

function [tout,fout] = linesearch(f,t1,t2,userOpts,f1)

% LINESEARCH : minimizes a 1-d function over an interval
% using a backstepping procedure
%
% Usage : t = linesearch(f,tmin,tmax)
%
% where f : R -> R is a 1 dimensional function, and 
% t is argmin_t f(t) over the interval.
%

% This function implements an approximate forward/backstep procedure

%set the default search options, and merge
%in the user ones
opts.tol              = 1e-5;
opts.backStepScaling  = 0.5;  %backtracking multiplier
if(nargin == 4)
    opts = setstructfields(opts,userOpts);
end

%small sanity check on options
if(opts.backStepScaling <= 0 || opts.backStepScaling >= 1)
    error('Line search scale factor should be in interval (0,1)');
end

%compute the value at t1 if not provided
if(nargin < 5)
    f1 = f(t1);
end
f2 = f(t2);  %value at nominal step size

%initial search window size
tdiff = t2-t1;

%create a list of points evaluated so far
tList = [t1 t2];
fList = [f1 f2];

%Try a forward scaling search
didForwardSearch = 0;
if(f2 < f1)
    didForwardSearch = true;
    [tF,fF,tFex,fFex] = scalingSearch(t1,t2,f,1./opts.backStepScaling,f2);
    tList = [tList tF tFex];
    fList = [fList fF fFex];
end

%if forward search wasn't used,
%try a backwards scaling search
if(~didForwardSearch)
    while((t2-t1)/tdiff > opts.tol)
        if(f2 < f1)
            %it worked!  Continue for best improvement
            [tB,fB,tBex,fBex] = scalingSearch(t1,t2,f,opts.backStepScaling,f2);
            tList = [tList tB tBex];
            fList = [fList fB fBex];
            break;
        else
            %reduce by factor c and try again
            t2 = t1 + (t2-t1)*opts.backStepScaling;
            f2 = f(t2);
        end
    end
end

%One last check - try to fit a quadratic
%function to the points found so far
tQ = quadraticSearch(tList,fList);
fQ = f(tQ);

%assemble everything and find the approximate minimizer
tList = [tList tQ];
fList = [fList fQ];
    
[fout,idx] = min(fList);
tout = tList(idx);

end
%----------------------------------------
%----------------------------------------

function [tout,fout,tExtreme,fExtreme] = scalingSearch(t1,t2,f,c,f2Init)

%helper function that keeps scaling until 
%the improvement stops

f2 = f2Init; f2last = realmax;

while(f2 <= f2last)
    f2last = f2; t2last = t2;
    %scale by a factor c and try again
    t2 = t1 + (t2-t1).*c;
    f2 = f(t2);
end

%the approximate optimizer
tout = t2last;
fout = f2last;

%the last point evaluated
tExtreme = t2;
fExtreme = f2;

end
%----------------------------------------
%----------------------------------------

function tOpt = quadraticSearch(t,f)

%Try to fit a quadratic function and
%find the minimizer
A = [t(:).^2 t(:) ones(size(t(:)))];

if(rank(A)<3)
    tOpt = t(1);
else
    c = A\f(:);
    tOpt = -c(2)/c(1)/2;
end

%further protection against ill conditioning
%and going backwards
if(isinf(tOpt) || isnan(tOpt) || tOpt < min(t))
    tOpt = t(1);
end
end

%
%
%

%% Orth
function Q = orth(A)

%ORTH   Orthogonalization.
%   Q = ORTH(A) is an orthonormal basis for the range of A.
%   That is, Q'*Q = I, the columns of Q span the same space as 
%   the columns of A, and the number of columns of Q is the 
%   rank of A.

% This code is copied from the standard
% matlab function of the same name, but with
% a different tolerance setting.  The matlab
% code does not seem able to identify which 
% columns of the SVD are a basis for the image,
% because there are a large number of singular
% values just over the canned matlab threshold
%
% This seems to be a particular problem when
% finding the image for a projected matrix, i.e.
% something like (I-UU')A, where U is orthonormal

[U,S] = svd(A,0);
[m,n] = size(A);
if m > 1, s = diag(S);
   elseif m == 1, s = S(1);
   else s = 0;
end

%mathworks value
%tol = max(m,n) * max(s) * eps(class(A));

%new value
tol = sqrt(eps);

r = sum(s > tol);
Q = U(:,1:r);
end
%
%}