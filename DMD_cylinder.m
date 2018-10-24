clear all
close all
load('EstimationSnapRed','GlobalSnap','visx','visy');
%load estimateionsnapred with variables: global snap, visx, visy
snapshotsU = GlobalSnap(1:size(GlobalSnap,1)/2,:);
snapshotsV = GlobalSnap(size(GlobalSnap,1)/2+1:size(GlobalSnap,1),:);
%u and v components o GlobalSnap

clear GlobalSnap
% plot of snapshot to view actual flow
% for i=1:100:1000
% 
% figure(i)
% scatter(visx, visy, [], snapshotsU(:,i), 'filled');
% title('snapshotsU');
% axis equal
% end

% figure(2)
% scatter(visx,visy,[],snapshotsV(:,1),'filled');
% title('snapshotsV');
% axis equal
%this are scatter of snapshotU and V at t=1. 
%ex: snapshotsU(:,1) = snapshot of U at t=1


%% plots
for k=10 %change k to change "t" of the snapshot evaulated
%select snapshots to plot
visu = snapshotsU(:, k:k+100);%just so you can select a portion of the data to analyse
visv = snapshotsV(:, k:k+100); %snapshotsV from t= k to k+10(dt)

%visu = snapshotsU(:, k:10:k+100);

%make uniform grid
xdummy = -0.98:0.05:14.9;
ydummy = -2.9:0.025:2.9;

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


%interpolate data onto new coordinates
% size(visu,2) == at t=2
interpolated_visu = zeros(nl, size(visu,2));
for i = 1:size(visu,2)
    interpolated_visu(:,i) = griddata(visx, visy, visu(:,i), xnew, ynew);
end
% 
interploated_visv=zeros(nl,size(visv,2));
for i = 1:size(visv,2)
    interpolated_visv(:,i) = griddata(visx,visy,visv(:,i),xnew,ynew);
end


%prepare data into form for contour plot
p = length(xdummy);
n = length(ydummy);
xCon = reshape(xnew,[n,p]);
yCon = reshape(ynew,[n,p]);


%actual plots
for i = 1:size(visu,2)
    dLamC = reshape(interpolated_visu(:,i), [n,p]); %reshape inter_visu to [n,p]
    figure (1)
    contourf(xCon,yCon,dLamC, 8);
    hold all
    pos = [-0.5 -0.5 1 1];
    rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'w');
    %position of the cylinder
    axis equal
    axis([-0.94 14.9 -2.9 2.9]);
    colormap gray
end
% 
% for i = 1:size(visv,2)    
%     dLamD = reshape(interpolated_visv(:,i), [n,p]);
%     figure (2)
%     contourf(xCon,yCon,dLamD, 8);
%     hold all
%     pos = [-0.5 -0.5 1 1];
%     rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'w');
%     axis equal
%     axis([-0.94 14.9 -2.9 2.9]);
%     colormap gray
% end

end

%% smooth plots - hig res

for i = 1:size(visu,2)
    
    dLamC = reshape(interpolated_visu(:,i), [n,p]);
    figure
    scatter(xnew, ynew, [], interpolated_visu(:,i), 'filled');
    hold all
    pos = [-0.5 -0.5 1 1];
    rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'w');%this line just plots the cylinder
    axis equal
    axis([-0.94 14.9 -2.9 2.9]);
end

for i = 1:size(visv,2)
    
    dLamD = reshape(interpolated_visv(:,i), [n,p]);
    figure
    scatter(xnew, ynew, [], interpolated_visv(:,i), 'filled');
    hold all
    pos = [-0.5 -0.5 1 1];
    rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'w');%this line just plots the cylinder
    axis equal
    axis([-0.94 14.9 -2.9 2.9]);
end

%% DMD - u
B = interpolated_visu(:,1:end-1);
A = interpolated_visu(:,2:end);
% 'Before' and 'After' snapshots

[U,S,V] = svd(B,'econ'); 
F_dmd = U'*A*V/S;
[Y,D_nu]=eig(F_dmd);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd

Vand=zeros(size(F_dmd));
for loop1=1:size(B,2)
    Vand(:,loop1)=diag(D_nu).^(loop1-1); %creating vandermonde with increasing powers
end

Psi = U*Y;

%% alpha - POD
alpha = U'*B;

k = 20 %number of modes

plot alpha

for i=1:20 %size(alpha,2)
    figure
    plot(alpha(i,:));
end

%mode energy plot
bar(diag(S(1:k,1:k)));
set(gca,'Yscale','log');


%% alpha - DMD

k= 20; %number of modes

X(:,1) =(Y^-1)*U'*B(:,1);


for i = 1:size(Vand,1)
    for j = 1:size(Vand,1)
        alpha_dmd(j,i) = X(j,1)*Vand(j,i);
    end
end

[U_dmd,S_dmd,V_dmd]=svd(U'*B,'econ'); %singular value of U'B, like POD with B

%comparing S_dmd to S
figure
subplot(2,1,1)
bar(diag(S(1:k,1:k)));
title('energy mode POD');
set(gca,'Yscale','log');
subplot(2,1,2)
bar(diag(S_dmd(1:k,1:k)));
title('energy mode DMD');
set(gca,'Yscale','log');

%comparing alpha_dmd against alpha_pod
for i = 1:k 
    figure 
    subplot(2,1,1)
    plot (alpha(i,:));
    title('\alpha POD');
    subplot(2,1,2)
    plot(real(alpha_dmd(i,:)));
    title('\alpha DMD');
end

%% eigenvalues

D = diag(D_nu);
plot(D,'o');
%circle implementation
    hold on 
    th=0:pi/50:2*pi;
    plot(cos(th),sin(th));
    hold off
xlim([0,1.5]);
ylim([-1.5,1.5]);
xlabel('Re(\lambda_i)');
ylabel('Im(\lambda_i)');
title('DMD eigenvalues)');

%% plot
% 
% for i = 1:k %size(Psi,2) %number of modes
%     DLamC = reshape(Psi(:,i),[n,p]);
%     figure('name','DMD-u');
%     subplot(3,1,1)
%     contourf(xCon,yCon,real(DLamC), 8);
%     title(['Real: Mode',num2str(i-1)]);
%     hold all
%     pos = [-0.5 -0.5 1 1];
%     rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'w');
%     %position of the cylinder
%     axis equal
%     axis([-0.94 14.9 -2.9 2.9]);
%     colormap gray
%     
%     subplot(3,1,2)
%     contourf(xCon,yCon,imag(DLamC), 8);
%     title(['Imag: Mode',num2str(i-1)]);
%     hold all
%     pos = [-0.5 -0.5 1 1];
%     rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'w');
%     %position of the cylinder
%     axis equal
%     axis([-0.94 14.9 -2.9 2.9]);
%     colormap gray
%     
%     subplot(3,1,3)
%     plot(real(alpha_dmd(i,:)));
% end

%% freq
dt = 1/720;
freq = 2*pi*imag(log(diag(D_nu)))/dt;

for i = 1:size(alpha_dmd,2)
    alpha_dmd_norm(i) = norm(real(alpha_dmd(i,:)));
end
alpha_dmd_norm = alpha_dmd_norm/norm(alpha_dmd_norm);

for i = 1:size(alpha,2)
    alpha_pod_norm(i) = norm(alpha(i,:));
end
t = alpha_pod_norm/norm(alpha_pod_norm);

figure(1)
plot(freq,alpha_dmd_norm,'o');
xlabel('frequency');
ylabel('||\alpha_{DMD}||');

figure (2)
plot(alpha_dmd_norm);
xlabel('n');
ylabel('||\phi_{DMD}||');

%% fourier transform
for i = 1:k
a = (((Psi(:,i).')*Psi(:,i)).^-1)*Psi(:,i).'*B;

Fs = 720*0.1965/15; %sampling frequency
period = 0:1/Fs:Fs*500;
N = length(a);
xdft = fft(a);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N))*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq_ = 0:Fs/length(a):Fs/2;

hold all
figure(3)
plot(freq_,psdx);
set(gca,'Xscale','log');
end



