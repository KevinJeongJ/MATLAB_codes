clear all
close all
%% load data
load('EstSnap13.mat');

%% flow plot
fig = figure;
set(gcf,'pos',[10 10 5*500 500]);

k=1;
for i = 1:10:51
    subplot(1,5,k);
    [p_e(:,k),r_e,th_e] = interpP(p(:,i),80,90);
    [R(k),Theta(k)] = calc_cop_t(p(:,i));
    th_e = rescale(th_e,0,2*pi);
    dLamC = reshape(p_e(:,k), size(r_e,2),size(th_e,2));
    polarcont(r_e,th_e,dLamC,8);
    hold on
    [X(k),Y(k)] = pol2cart(Theta(k),R(k));
    plot(X(k),Y(k),'r*');  
%     saveas(figure(k), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_pressure/flow_w_cop/', ['figure' num2str(k) '.jpeg']));
    colormap(gray)
    k=k+1;
end
colorbar('position',[1150 10 10 200]);
%{

%% DMD
for i = 1:size(p,2)/20
    [p_e(:,i),r_e,th_e] = interpP(p(:,i),80,90);
    [R(i),Theta(i)] = calc_cop_t(p_e(:,i));
end

B = p_e(:,1:end-1);
A = p_e(:,2:end);

[U,S,V] = svd(B,'econ');
F_dmd = U'*A*V/S;
[Y,D_nu]=eig(F_dmd);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd

Vand=zeros(size(F_dmd));
for loop1=1:size(B,2)
    Vand(:,loop1)=diag(D_nu).^(loop1-1); %creating vandermonde with increasing powers
end

Psi = U*Y;
Z(:,1) =(Y^-1)*U'*B(:,1);

for i = 1:size(Vand,1)
    for j = 1:size(Vand,1)
        alpha_dmd(j,i) = Z(j,1)*Vand(j,i);
    end
end

[U_dmd,S_dmd,V_dmd]=svd(U'*B,'econ'); %singular value of U'B, like POD with B


alpha_dmd_abs_I = abs(real(alpha_dmd));
alpha_dmd_abs_I_t = alpha_dmd_abs_I.';

I = sum(alpha_dmd_abs_I_t);
[C,T] = sort(I,'descend');

%% DMD plot
%  for i=1:30 %size(alpha_dmd,2)
%     k = T(i);
%     DLamD = reshape(Psi(:,k), size(r_e,2),size(th_e,2));
%     th_e = rescale(th_e,0,2*pi);
%     figure('name','DMD','pos',[10,10,600,900])
%     subplot(3,3,[1,2,3,4,5,6])
%     polarcont(r_e,th_e,real(DLamD),10);
%     colormap gray
% 
%     subplot(3,3,[7,8,9])
%     plot(real(alpha_dmd(k,:)));
%     xlim([0,size(Psi,2)]);
%     
%     saveas(figure(i), fullfile('/home/kevin/Desktop/Data/bullet/bullet_figure/bullet_pressure/bullet_pressure_dmd//', ['figure' num2str(i) '.jpeg']));
% end


%}