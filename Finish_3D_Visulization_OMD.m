clear all;
load('02_5omd.mat');

clearvars -except DLamSt St Theta iptx_f xCon yCon

%% if Theta is lost from the varaibles 
% maxchanges = 80;
% ipt = findchangepts(Theta,'MaxNumChanges',maxchanges);
% iptx = [1; ipt(:); size(Theta,2)];
% clear ipt
% 
% Theta_f1 = zeros(size(Theta));
% for i = 1:maxchanges
%     if size(Theta(iptx(i):iptx(i+1)),2) < 500 
%         Theta_f1(iptx(i):iptx(i+1)) = 0; 
%     else
%         Theta_f1(iptx(i):iptx(i+1)) = Theta(iptx(i):iptx(i+1));
%     end
% end
% 
% k=1;
% for i=1:size(iptx,1)-1
%     if (iptx(i+1)-iptx(i))>=500
%         iptx_f(k) = iptx(i);
%         iptx_f(k+1)= iptx(i+1);
%         k = k+2;
%     end
% end
% clear Theta iptx
% 
% c= 1;
% for i = 1:length(iptx_f)/2
%     a = iptx_f((2*i)-1);
%     b = iptx_f((2*i))-1; 
%     d = c+(b-a);
%     Theta(i) = mean(Theta_f1(c:d));
%     c = c+(b-a);
% end




%% filter - only top flow 
%this section effectively only considers the top half of the mode (2D mode) 
disp('filter');
level =4;
for i = 1:length(iptx_f)/2
    %%%disk filter
    H = fspecial('disk',1);
    dLamC = imfilter(real(DLamSt(:,:,i)),H,'replicate');
    dLamD = imfilter(imag(DLamSt(:,:,i)),H,'replicate'); 
    DLamC(:,:,i) = dLamC;
    DLamD(:,:,i) = dLamD;
end

%%
%percentage filter, only consider chosen percentages of 
%data from the max and min are considered

per = 0.5; %percentile
disp('10%filter')

DLamCmax = max(max(max(DLamC)));
DLamCmin = min(min(min(DLamC)));
DLamDmax = max(max(max(DLamD)));
DLamDmin = min(min(min(DLamD)));

DLamC(per*DLamCmax>DLamC & DLamC>per*DLamCmin) = 0;
DLamD(per*DLamDmax>DLamD & DLamD>per*DLamDmin) = 0;

clear DLamCmax DLamCmin DLamDmax DLamDmin

%% color filter
% to adjust for the prevention of "color-flipping" mechanism 
for i = 1:size(DLamC,3)
    format longE
    if mean(mean(DLamC(:,:,i)))<0
        DLamC_f(:,:,i) = -DLamC(:,:,i);
    else
        DLamC_f(:,:,i) = DLamC(:,:,i);
    end
    
    if mean(mean(DLamD(:,:,i)))<0
        DLamD_f(:,:,i) = -DLamD(:,:,i);
    else
        DLamD_f(:,:,i) = DLamD(:,:,i);
    end
end

%%
%{
for i =14%1:length(iptx_f)/2
    fig = figure('Visible','on');
    set(gcf,'pos',[10,10,3600,400]);
   
    hold all 
    cx1 = subplot(1,4,1)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['Theta:',num2str(Theta(i)),' a) Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(1,4,2)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['b) Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    colormap gray
%     saveas(fig, fullfile('/home/clear all;
load('02_5omd.mat');

clearvars -except DLamSt St Theta iptx_f xCon yCon


% maxchanges = 80;
% ipt = findchangepts(Theta,'MaxNumChanges',maxchanges);
% iptx = [1; ipt(:); size(Theta,2)];
% clear ipt
% 
% Theta_f1 = zeros(size(Theta));
% for i = 1:maxchanges
%     if size(Theta(iptx(i):iptx(i+1)),2) < 500 
%         Theta_f1(iptx(i):iptx(i+1)) = 0; 
%     else
%         Theta_f1(iptx(i):iptx(i+1)) = Theta(iptx(i):iptx(i+1));
%     end
% end
% 
% k=1;
% for i=1:size(iptx,1)-1
%     if (iptx(i+1)-iptx(i))>=500
%         iptx_f(k) = iptx(i);
%         iptx_f(k+1)= iptx(i+1);
%         k = k+2;
%     end
% end
% clear Theta iptx
% 
% c= 1;
% for i = 1:length(iptx_f)/2
%     a = iptx_f((2*i)-1);
%     b = iptx_f((2*i))-1; 
%     d = c+(b-a);
%     Theta(i) = mean(Theta_f1(c:d));
%     c = c+(b-a);
% end




%% filter - only top flow 
disp('filter');
level =4;
for i = 1:length(iptx_f)/2
    %%%disk filter
    H = fspecial('disk',1);
    dLamC = imfilter(real(DLamSt(:,:,i)),H,'replicate');
    dLamD = imfilter(imag(DLamSt(:,:,i)),H,'replicate'); 
    DLamC(:,:,i) = dLamC;
    DLamD(:,:,i) = dLamD;
end

%%
per = 0.5; %percentile
disp('10%filter')

DLamCmax = max(max(max(DLamC)));
DLamCmin = min(min(min(DLamC)));
DLamDmax = max(max(max(DLamD)));
DLamDmin = min(min(min(DLamD)));

DLamC(per*DLamCmax>DLamC & DLamC>per*DLamCmin) = 0;
DLamD(per*DLamDmax>DLamD & DLamD>per*DLamDmin) = 0;

clear DLamCmax DLamCmin DLamDmax DLamDmin
%% color filter
for i = 1:size(DLamC,3)
    format longE
    if mean(mean(DLamC(:,:,i)))<0
        DLamC_f(:,:,i) = -DLamC(:,:,i);
    else
        DLamC_f(:,:,i) = DLamC(:,:,i);
    end
    
    if mean(mean(DLamD(:,:,i)))<0
        DLamD_f(:,:,i) = -DLamD(:,:,i);
    else
        DLamD_f(:,:,i) = DLamD(:,:,i);
    end
end

%% Mode visualization
%{
for i =14%1:length(iptx_f)/2
    fig = figure('Visible','on');
    set(gcf,'pos',[10,10,3600,400]);
   
    hold all 
    cx1 = subplot(1,4,1)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['Theta:',num2str(Theta(i)),' a) Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(1,4,2)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['b) Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end


for i =3%1:length(iptx_f)/2
    cx1 = subplot(1,4,3)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['Theta:',num2str(Theta(i)),' c) Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(1,4,4)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['d) Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end



hold off
clear cx1 cx2 cx3 cx4 cx5 cx6 cx7 cx8




%%
for i =14%1:length(iptx_f)/2
    fig = figure('Visible','on');
    set(gcf,'pos',[10,10,3600,800]);
   
    hold all 
    cx1 = subplot(2,4,1)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['a) Real, level 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(2,4,5)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['e) Imag, level = 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % lower contour level
    cx3 = subplot(2,4,2)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level);
    title(['b) Real, level 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx4 = subplot(2,4,6)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level);
    title(['f) Imag, level = 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % filter data 
    cx5 = subplot(2,4,3)
    contourf(xCon/200,yCon/200,DLamC(:,:,i), level);
    title(['c) Real, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%       ylim([0.1,inf]);    
    caxis([-20*10^-3 20*10^-3]);


    cx6 = subplot(2,4,7)
    contourf(xCon/200,yCon/200,DLamD(:,:,i),level);
    title(['g) Imag, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    colormap gray
    caxis([-20*10^-3 20*10^-3]);

    
    % corrected filter data
    cx7 = subplot(2,4,4)
    contourf(xCon/200,yCon/200,DLamC_f(:,:,i), level);
    title(['d) Real, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);

    cx8 = subplot(2,4,8)
    contourf(xCon/200,yCon/200,DLamD_f(:,:,i),level);
    title(['h) Imag, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);


%       if Theta(i)>pi
%           ylim([cx5,cx6,cx7,cx8],[-inf, 0]);
%       else
%           ylim([cx5,cx6,cx7,cx8],[0, inf]);
%       end


    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end
hold off
clear cx1 cx2 cx3 cx4 cx5 cx6 cx7 cx8

%%
k =1;
fig = figure('Visible','on');
set(gcf,'pos',[10,10,1000,150]);

for i =13;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    colormap gray

end

k=2; 

for i =4;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    

    colormap gray

end


%{
%% slice images
disp('slices');
for i =14%1:length(iptx_f)/2
    fig = figure('Visible','on');
    set(gcf,'pos',[10,10,3600,800]);
   
    hold all 
    cx1 = subplot(2,4,1)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['a) Real, level 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(2,4,5)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['e) Imag, level = 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % lower contour level
    cx3 = subplot(2,4,2)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level);
    title(['b) Real, level 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx4 = subplot(2,4,6)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level);
    title(['f) Imag, level = 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % filter data 
    cx5 = subplot(2,4,3)
    contourf(xCon/200,yCon/200,DLamC(:,:,i), level);
    title(['c) Real, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%       ylim([0.1,inf]);    
    caxis([-20*10^-3 20*10^-3]);


    cx6 = subplot(2,4,7)
    contourf(xCon/200,yCon/200,DLamD(:,:,i),level);
    title(['g) Imag, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    colormap gray
    caxis([-20*10^-3 20*10^-3]);

    
    % corrected filter data
    cx7 = subplot(2,4,4)
    contourf(xCon/200,yCon/200,DLamC_f(:,:,i), level);
    title(['d) Real, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);

    cx8 = subplot(2,4,8)
    contourf(xCon/200,yCon/200,DLamD_f(:,:,i),level);
    title(['h) Imag, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);


%       if Theta(i)>pi
%           ylim([cx5,cx6,cx7,cx8],[-inf, 0]);
%       else
%           ylim([cx5,cx6,cx7,cx8],[0, inf]);
%       end


    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end
hold off
clear cx1 cx2 cx3 cx4 cx5 cx6 cx7 cx8

%%
k =1;
fig = figure('Visible','on');
set(gcf,'pos',[10,10,1000,300]);

for i =13;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    colormap gray

end

k=2; 

for i =4;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    

    colormap gray

end

%}
%}

%% image matching between real and imag
disp('image/color matching');
dummy1 = [St.',Theta.',[1:size(DLamSt,3)].'];
dummy1 = dummy1(dummy1(:,2)<2.3,:); %either >3.3 or <1.8 modes (switch ylim accordingly)
dummy1 = dummy1(dummy1(:,1)>.15 & dummy1(:,1)<.25,:);
dummy1 = sortrows(dummy1,2);
%0.03 0.08 for 0.06
%0.15 0.25 for 0.2
%0.3 0.5 for 0.4
dummy2 = [St.',Theta.',[1:size(DLamSt,3)].'];
dummy2 = dummy2(dummy2(:,2)>3.3,:); %either >3.3 or <1.8 modes (switch ylim accordingly)
dummy2 = dummy2(dummy2(:,1)>.15 & dummy2(:,1)<.25,:);
dummy2 = sortrows(dummy2,2);

DLamF1(:,:,1) = DLamC_f((end+1)/2:end,:,dummy1(1,3)); %starting mode
DLamG1(:,:,1) = DLamC_f(1:(end+1)/2,:,dummy1(1,3)); %starting mode


ImDLamF1(:,:,1) = DLamD_f((end+1)/2:end,:,dummy1(1,3)); %starting mode
ImDLamG1(:,:,1) = DLamD_f(1:(end+1)/2,:,dummy1(1,3)); %starting mode

% switch the focus to DLamF1 if need be and vice versa
for i = 2:size(dummy1,1)
    k = dummy1(i,3);    
    
        DLamCdummy = DLamC_f(1:((end+1)/2),:,k);
        DLamDdummy = DLamD_f(1:((end+1)/2),:,k);
    
    DLamCdummy2 = DLamCdummy/norm(DLamCdummy);
    DLamDdummy2 = DLamDdummy/norm(DLamDdummy);
    
    DLamFdummy = DLamG1(:,:,i-1)/norm(DLamG1(:,:,i-1));
    
    if sum(max(abs(DLamFdummy-DLamCdummy2)))<= sum(max(abs(DLamFdummy-DLamDdummy2)))
        DLamG1(:,:,i) = DLamCdummy;
        DLamF1(:,:,i) = DLamC_f((end+1)/2:end,:,k);
        ImDLamG1(:,:,i) = DLamDdummy;
        ImDLamG1(:,:,i) = DLamD_f((end+1)/2:end,:,k);
    else
        DLamG1(:,:,i) = DLamDdummy;
        DLamF1(:,:,i) = DLamD_f((end+1)/2:end,:,k);
        ImDLamG1(:,:,i) = DLamCdummy;
        ImDLamF1(:,:,i) = DLamC_f((end+1)/2:end,:,k);
    end
end
clear DLamCdummy DLamCdummy2 DLamDdummy DLamDdummy2 DLamFdummy


%
DLamF2(:,:,1) = DLamC_f((end+1)/2:end,:,dummy2(1,3)); %starting mode
DLamG2(:,:,1) = DLamC_f(1:(end+1)/2,:,dummy2(1,3)); %starting mode

ImDLamF2(:,:,1) = DLamD_f((end+1)/2:end,:,dummy2(1,3)); %starting mode
ImDLamG2(:,:,1) = DLamD_f(1:(end+1)/2,:,dummy2(1,3)); %starting mode

for i = 2:size(dummy2,1)
    k = dummy2(i,3);    
    
        DLamCdummy = DLamC_f(1:(end+1)/2,:,k);
        DLamDdummy = DLamD_f(1:(end+1)/2,:,k);
    
    DLamCdummy2 = DLamCdummy/norm(DLamCdummy);
    DLamDdummy2 = DLamDdummy/norm(DLamDdummy);
    
    DLamFdummy = DLamG2(:,:,i-1)/norm(DLamG2(:,:,i-1));
    
    if sum(max(abs(DLamFdummy-DLamCdummy2)))<= sum(max(abs(DLamFdummy-DLamDdummy2)))
        DLamG2(:,:,i) = DLamCdummy;
        DLamF2(:,:,i) = DLamC_f((end+1)/2:end,:,k);
        ImDLamG2(:,:,i) = DLamDdummy;
        ImDLamF2(:,:,i) = DLamD_f((end+1)/2:end,:,k);
    else
        DLamG2(:,:,i) = DLamDdummy;
        DLamF2(:,:,i) = DLamD_f((end+1)/2:end,:,k);
        ImDLamG2(:,:,i) = DLamCdummy;
        ImDLamF2(:,:,i) = DLamC_f((end+1)/2:end,:,k);
    end
end
clear DLamCdummy DLamCdummy2 DLamDdummy DLamDdummy2 DLamFdummy

clear DLamCdummy DLamCdummy2 DLamDdummy DLamDdummy2 DLamFdummy



%%
fig = figure('visible','on');
set(gcf,'pos',[10,10,900*size(DLamF1,3),800]);

DLamMix = zeros(size(DLamC_f,1),size(DLamC_f,2),size(DLamG2,3));
for i= 1:size(DLamG2,3)
    DLamMix(1:(end+1)/2,:,i) = DLamG2(:,:,i);
    DLamMix((end+1)/2:end,:,i) = DLamF2(:,:,i);
end

for i = 1:size(DLamG2,3)
    subplot(2,((size(DLamF2,3))+1)/2,i)
    contourf(xCon/200,yCon/200,DLamMix(:,:,i),level);
    title(['\theta:' num2str(dummy2(i,2))],'fontweight','normal','FontSize',22,'fontname','arial');
    caxis([-20*10^-3 20*10^-3]);

    colormap(gray)
end
% 
% %  saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd//', ['filtered.jpeg']));



%% contourslice stack
% xNew = xCon((end+1)/2:end,:);
% yNew = yCon((end+1)/2:end,:);
% 
% for i=1:size(dummy2,1)-1
%     fig(i) = figure('visible','on');
%     hold all
%     contour(xNew/200,yNew/200,DLamG2(:,:,i),1)
%     saveas(fig(i), fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['Stack: Theta:',num2str(dummy2(i,2)),' St=',num2str(dummy2(i),1),'.jpeg']));
%     
%     contour(xNew/200,yNew/200,DLamG2(:,:,i+1),1)
%     saveas(fig(i), fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['Stack: Theta:',num2str(dummy2(i,2)),' St=',num2str(dummy2(i),1),'.jpeg']));
% 
%     close all
% end
% 
% 
% fig = figure('visible','on');
% for i=1:size(dummy2,1)
%     hold all
%     contour(xNew/200,yNew/200,DLamG2(:,:,i),1)
% end
% saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['ALL Stack: Theta:',num2str(dummy2(1,2)),'~',num2str(dummy2(end,2)),'.jpeg']));
% 
% fig = figure('visible','on');
% 
% for i=1:size(dummy2,1)
%     hold all
%     contour(xNew/200,yNew/200,DLamG2(:,:,i),1)
% end
% saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['180 ALL Stack: Theta:',num2str(dummy2(1,2)),'~',num2str(dummy2(end,2)),'.jpeg']));

%% setting up polar 
disp('polar');
N = size(DLamSt,3);

k = size(DLamSt,1);
R = zeros(k,1);
for i = 1:k
    R(i) = -19.65 +19.65*2/(k-1)*(i-1);
end 
R=R.';

%%
disp('no filter');

Theta_f1 = dummy1(:,2);
dummyG1 = [Theta_f1,[1:size(DLamF1,3)].'];
dummyG1 = sortrows(dummyG1,1);
R = R(R>=0);
Theta_f1 = dummyG1(:,1);

Theta_f2 = dummy2(:,2)-pi;
dummyG2 = [Theta_f2,[1:size(DLamF2,3)].'];
dummyG2 = sortrows(dummyG2,1);
R = R(R>=0);
Theta_f2 = dummyG2(:,1);

for i = 1:size(DLamF1,3)
    DLamL1(:,:,i) = DLamF1(:,:,dummyG1(i,2));
    ImDLamL1(:,:,i) = ImDLamF1(:,:,dummyG1(i,2));
end

for i = 1:size(DLamF2,3)
    DLamL2(:,:,i) = DLamF2(:,:,dummyG2(i,2));
    ImDLamL2(:,:,i) = ImDLamF2(:,:,dummyG2(i,2));

end

% for the other 180 degrees
Theta_fr1 = dummy1(:,2)+pi;
dummyGr1 = [Theta_fr1,[1:size(DLamG1,3)].'];
dummyGr1 = sortrows(dummyGr1,1);
Theta_fr1 = dummyGr1(:,1);

for i = 1:size(DLamG1,3)
    DLamLr1(:,:,i)=DLamG1(:,:,dummyGr1(i,2));
    ImDLamLr1(:,:,i)=ImDLamG1(:,:,dummyGr1(i,2));

end

Theta_fr2 = dummy2(:,2);
dummyGr2 = [Theta_fr2,[1:size(DLamG2,3)].'];
dummyGr2 = sortrows(dummyGr2,1);
Theta_fr2 = dummyGr2(:,1);

for i = 1:size(DLamG2,3)
    DLamLr2(:,:,i)=DLamG2(:,:,dummyGr2(i,2));
    ImDLamLr2(:,:,i)= ImDLamG2(:,:,dummyGr2(i,2));
end

% for i =1:size(DLamF1,2)
% %       not displaying figures
%     fig = figure('Visible','off');
%     set(gcf,'pos',[10,10,900,900]);
%     
%     hold all 
%     dLamL = reshape(DLamL1(:,i,:),size(DLamL1,1),size(DLamL1,3));
%     polarcont(R,Theta_f1,dLamL,level);    
%     colormap(gray)
%     saveas(fig, fullfile('/home/kevin/Desktop/meh/polar<3/<pi/', ['polarcont',num2str(i),'.jpeg']));
% 
% end
clear dLamF dLamG 
%% data shuffle (cart->polar
disp('cart->polar');
level = 4;
for i = 1:size(DLamF1,2)
    dLamL1 = reshape(DLamL1(:,i,:),size(DLamL1,1),size(DLamL1,3));
    dLamH3_1(:,i,:) = dLamL1;
    ImdLamL1 = reshape(ImDLamL1(:,i,:),size(ImDLamL1,1),size(ImDLamL1,3));
    ImdLamH3_1(:,i,:) = ImdLamL1;
    
end
for i = 1:size(DLamF2,2)
    dLamL2 = reshape(DLamL2(:,i,:),size(DLamL2,1),size(DLamL2,3));
    dLamH3_2(:,i,:) = dLamL2;
    ImdLamL2 = reshape(ImDLamL2(:,i,:),size(ImDLamL2,1),size(ImDLamL2,3));
    ImdLamH3_2(:,i,:) = ImdLamL2;
end
% for the other 180 degrees
for i = 1:size(DLamG1,2)
    dLamLr1 = reshape(DLamLr1(:,i,:),size(DLamLr1,1),size(DLamLr1,3));
    dLamH3r_1(:,i,:)=dLamLr1;
    ImdLamLr1 = reshape(ImDLamLr1(:,i,:),size(ImDLamLr1,1),size(ImDLamLr1,3));
    ImdLamH3r_1(:,i,:)=ImdLamLr1;
end
for i = 1:size(DLamG2,2)
    dLamLr2 = reshape(DLamLr2(:,i,:),size(DLamLr2,1),size(DLamLr2,3));
    dLamH3r_2(:,i,:)=dLamLr2;
    ImdLamLr2 = reshape(ImDLamLr2(:,i,:),size(ImDLamLr2,1),size(ImDLamLr2,3));
    ImdLamH3r_2(:,i,:)=ImdLamLr2;
end
close all
clear  dLamH dLamL1 dLamL2 dLamLr1 dLamLr2
clear  dLamH ImdLamL1 ImdLamL2 ImdLamLr1 ImdLamLr2

%% 3D - poly_int
disp('poly_int');

%plevel = polyfit level
plevel = 2;
buffer = 10/180*pi;

for i =1:size(DLamF1,2)
    dLamL1 = reshape(dLamH3_1(:,i,:),size(dLamH3_1,1),size(dLamH3_1,3));
    
    Theta_f_0 = linspace(-pi/4,min(Theta_f1)-buffer,20);
    
    Theta_ft = linspace(min(Theta_f1)-buffer,max(Theta_f1)+buffer,25*size(Theta_f1,1));
    
    ddLamG_r1 = zeros(length(dLamL1),25*size(Theta_f1,1));
    for k = 1:size(Theta_f1,1)
        [O,T] = sort(abs(Theta_ft-Theta_f1(k)));
        ddLamG_r1(:,T(1)) = dLamL1(:,k);
    end
    for l = 1:size(ddLamG_r1,1)
        for j = 1:size(Theta_f1)
            poly_r1(l,j+1)= dLamL1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f1)+2,1);
        Thetadummy(2:end-1) = Theta_f1.';
        Thetadummy(1) = min(Theta_f1)-buffer;
        Thetadummy(end) = max(Theta_f1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft);
        ddLamG_r1(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f1)+buffer,pi-pi/4,20);
    
    Theta_final1t = cat(2,Theta_f_0, Theta_ft,Theta_f_180);
    Filler = zeros(size(dLamL1,1),20);
    ddLamR1 = cat(2,Filler,ddLamG_r1,Filler);
    ddLamRz1(:,:,i) = ddLamR1;
end

clear Theta_ft ddLamG_r1 p1 y2 ddLamR1 poly_r1 Thetadummy

for i =1:size(DLamL1,2)   
    dLamLr1 = reshape(dLamH3r_1(:,i,:),size(dLamH3r_1,1),size(dLamH3r_1,3));

    
    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr1)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr1)-buffer,max(Theta_fr1)+buffer,25*size(Theta_fr1,1));
    
    ddLamG_rr1 = zeros(length(dLamLr1),25*size(Theta_fr1,1));
    for k = 1:size(Theta_fr1,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr1(k)));
        ddLamG_rr1(:,T(1)) = dLamLr1(:,k);
    end
    for l = 1:size(ddLamG_rr1,1)
        for j = 1:size(Theta_fr1)
            poly_r1(l,j+1)= dLamLr1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr1)+2,1);
        Thetadummy(2:end-1) = Theta_fr1.';
        Thetadummy(1) = min(Theta_fr1)-buffer;
        Thetadummy(end) = max(Theta_fr1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ddLamG_rr1(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr1)+buffer,2*pi-pi/4,20);
    
    Theta_finalr1t = cat(2,Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr1,1),20);    
    ddLamRr1 = cat(2,Filler, ddLamG_rr1,Filler);
    ddLamRrr1(:,:,i) = ddLamRr1;
    
end

clear ddLamG_rr1 p1 y2 poly_r1 ddLamRr1 Thetadummy 


 for i =1:size(DLamL2,2)

    dLamL2 = reshape(dLamH3_2(:,i,:),size(dLamH3_2,1),size(dLamH3_2,3));

    Theta_f_0 = linspace(-pi/4,min(Theta_f2)-buffer,20);
    
    Theta_ft2 = linspace(min(Theta_f2)-buffer,max(Theta_f2)+buffer,20*size(Theta_f2,1));
    
    ddLamG_r2 = zeros(length(dLamL2),20*size(Theta_f2,1));
    for k = 1:size(Theta_f2,1)
        [O,T] = sort(abs(Theta_ft2-Theta_f2(k)));
        ddLamG_r2(:,T(1)) = dLamL2(:,k);
    end
    
    for l = 1:size(ddLamG_r2,1)
        for j = 1:size(Theta_f2)
            poly_r1(l,j+1)= dLamL2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f2)+2,1);
        Thetadummy(2:end-1) = Theta_f2.';
        Thetadummy(1) = min(Theta_f2)-buffer;
        Thetadummy(end) = max(Theta_f2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft2);
        ddLamG_r2(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f2)+buffer,pi-pi/4,20);
 

    Theta_final2t = cat(2,Theta_f_0, Theta_ft2,Theta_f_180);
    Filler = zeros(size(dLamL2,1),20);
    ddLamR2 = cat(2,Filler,ddLamG_r2,Filler);
    ddLamRz2(:,:,i) = ddLamR2;
 end

 clear ddLamG_r2 y2 p1 ddLamR2 Thetadummy poly_r1

for i =1:size(DLamL2,2)
    dLamLr2 = reshape(dLamH3r_2(:,i,:),size(dLamH3r_2,1),size(dLamH3r_2,3));

    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr2)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr2)-buffer,max(Theta_fr2)+buffer,20*size(Theta_fr2,1));
    
    ddLamG_rr2 = zeros(length(dLamLr2),20*size(Theta_fr2,1));
    for k = 1:size(Theta_fr2,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr2(k)));
        ddLamG_rr2(:,T(1)) = dLamLr2(:,k);
    end
     for l = 1:size(ddLamG_rr2,1)
        for j = 1:size(Theta_fr2)
            poly_r1(l,j+1)= dLamLr2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_fr2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr2)+2,1);
        Thetadummy(2:end-1) = Theta_fr2.';
        Thetadummy(1) = min(Theta_fr2)-buffer;
        Thetadummy(end) = max(Theta_fr2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ddLamG_rr2(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr2),2*pi-pi/4,20);
    
    Theta_finalr2t = cat(2, Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr2,1),20);    
    ddLamRr2 = cat(2,Filler, ddLamG_rr2,Filler);
    ddLamRrr2(:,:,i) = ddLamRr2;    
end

clear p1 y2 ddLamG_r2 ddLamRr2
clear Theta_fr_2pi Theta_f_180 Theta_f_0 Theta_180_flr
clear poly_r1 O T Thetadummy

%% ^ Im

disp('Im poly_int');

%plevel = polyfit level
plevel = 2;
buffer = 10/180*pi;

for i =1:size(ImDLamF1,2)
    ImdLamL1 = reshape(ImdLamH3_1(:,i,:),size(ImdLamH3_1,1),size(ImdLamH3_1,3));
    
    Theta_f_0 = linspace(-pi/4,min(Theta_f1)-buffer,20);
    
    Theta_ft = linspace(min(Theta_f1)-buffer,max(Theta_f1)+buffer,25*size(Theta_f1,1));
    
    ImddLamG_r1 = zeros(length(ImdLamL1),25*size(Theta_f1,1));
    for k = 1:size(Theta_f1,1)
        [O,T] = sort(abs(Theta_ft-Theta_f1(k)));
        ImddLamG_r1(:,T(1)) = ImdLamL1(:,k);
    end
    for l = 1:size(ImddLamG_r1,1)
        for j = 1:size(Theta_f1)
            poly_r1(l,j+1)= ImdLamL1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f1)+2,1);
        Thetadummy(2:end-1) = Theta_f1.';
        Thetadummy(1) = min(Theta_f1)-buffer;
        Thetadummy(end) = max(Theta_f1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft);
        ImddLamG_r1(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f1)+buffer,pi-pi/4,20);
    
    Theta_final1t = cat(2,Theta_f_0, Theta_ft,Theta_f_180);
    Filler = zeros(size(dLamL1,1),20);
    ImddLamR1 = cat(2,Filler,ImddLamG_r1,Filler);
    ImddLamRz1(:,:,i) = ImddLamR1;
end

clear Theta_ft ImddLamG_r1 p1 y2 ImddLamR1 poly_r1 Thetadummy

for i =1:size(ImDLamL1,2)   
    ImdLamLr1 = reshape(ImdLamH3r_1(:,i,:),size(ImdLamH3r_1,1),size(ImdLamH3r_1,3));

    
    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr1)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr1)-buffer,max(Theta_fr1)+buffer,25*size(Theta_fr1,1));
    
    ImddLamG_rr1 = zeros(length(ImdLamLr1),25*size(Theta_fr1,1));
    for k = 1:size(Theta_fr1,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr1(k)));
        ImddLamG_rr1(:,T(1)) = dLamLr1(:,k);
    end
    for l = 1:size(ImddLamG_rr1,1)
        for j = 1:size(Theta_fr1)
            poly_r1(l,j+1)= ImdLamLr1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr1)+2,1);
        Thetadummy(2:end-1) = Theta_fr1.';
        Thetadummy(1) = min(Theta_fr1)-buffer;
        Thetadummy(end) = max(Theta_fr1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ImddLamG_rr1(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr1)+buffer,2*pi-pi/4,20);
    
    Theta_finalr1t = cat(2,Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr1,1),20);    
    ImddLamRr1 = cat(2,Filler, ImddLamG_rr1,Filler);
    ImddLamRrr1(:,:,i) = ImddLamRr1;
    
end

clear ImddLamG_rr1 p1 y2 poly_r1 ImddLamRr1 Thetadummy 


 for i =1:size(ImDLamL2,2)

    ImdLamL2 = reshape(ImdLamH3_2(:,i,:),size(ImdLamH3_2,1),size(ImdLamH3_2,3));

    Theta_f_0 = linspace(-pi/4,min(Theta_f2)-buffer,20);
    
    Theta_ft2 = linspace(min(Theta_f2)-buffer,max(Theta_f2)+buffer,20*size(Theta_f2,1));
    
    ImddLamG_r2 = zeros(length(ImdLamL2),20*size(Theta_f2,1));
    for k = 1:size(Theta_f2,1)
        [O,T] = sort(abs(Theta_ft2-Theta_f2(k)));
        ImddLamG_r2(:,T(1)) = ImdLamL2(:,k);
    end
    
    for l = 1:size(ImddLamG_r2,1)
        for j = 1:size(Theta_f2)
            poly_r1(l,j+1)= ImdLamL2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f2)+2,1);
        Thetadummy(2:end-1) = Theta_f2.';
        Thetadummy(1) = min(Theta_f2)-buffer;
        Thetadummy(end) = max(Theta_f2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft2);
        ImddLamG_r2(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f2)+buffer,pi-pi/4,20);
 

    Theta_final2t = cat(2,Theta_f_0, Theta_ft2,Theta_f_180);
    Filler = zeros(size(dLamL2,1),20);
    ImddLamR2 = cat(2,Filler,ImddLamG_r2,Filler);
    ImddLamRz2(:,:,i) = ImddLamR2;
 end

 clear ImddLamG_r2 y2 p1 ImddLamR2 Thetadummy poly_r1

for i =1:size(ImDLamL2,2)
    ImdLamLr2 = reshape(ImdLamH3r_2(:,i,:),size(ImdLamH3r_2,1),size(ImdLamH3r_2,3));

    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr2)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr2)-buffer,max(Theta_fr2)+buffer,20*size(Theta_fr2,1));
    
    ImddLamG_rr2 = zeros(length(ImdLamLr2),20*size(Theta_fr2,1));
    for k = 1:size(Theta_fr2,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr2(k)));
        ImddLamG_rr2(:,T(1)) = ImdLamLr2(:,k);
    end
     for l = 1:size(ImddLamG_rr2,1)
        for j = 1:size(Theta_fr2)
            poly_r1(l,j+1)= ImdLamLr2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_fr2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr2)+2,1);
        Thetadummy(2:end-1) = Theta_fr2.';
        Thetadummy(1) = min(Theta_fr2)-buffer;
        Thetadummy(end) = max(Theta_fr2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ImddLamG_rr2(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr2),2*pi-pi/4,20);
    
    Theta_finalr2t = cat(2, Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr2,1),20);    
    ImddLamRr2 = cat(2,Filler, ImddLamG_rr2,Filler);
    ImddLamRrr2(:,:,i) = ImddLamRr2;    
end

clear p1 y2 ddLamG_r2 ddLamRr2
clear Theta_fr_2pi Theta_f_180 Theta_f_0 Theta_180_flr
clear poly_r1 O T Thetadummy

%{
%% polar slice disp
disp('filter&nofilter');

dmin_1 = min(min(min(cat(1,ddLamRz1,ddLamRrr1))));
dmax_1 = max(max(max(cat(1,ddLamRz1,ddLamRrr1))));
dmin_2 = min(min(min(cat(1,ddLamRz2,ddLamRrr2))));
dmax_2 = max(max(max(cat(1,ddLamRz2,ddLamRrr2))));

for i =100:10:180%:size(DLamF1,2)
%     %   not displaying figures
    level = 2;

    fig = figure('Visible','off');
    set(gcf,'pos',[10,10,2000,500]);
    holder = [-19.65;19.65];
    Theta_final_1t = cat(2,Theta_final1t,Theta_finalr1t);
    Theta_final_2t = cat(2,Theta_final2t,Theta_finalr2t);
    location = max(max(xCon/200))/size(DLamF1,2)*i;
%     ax1=subplot(2,4,1);
%     contourf(xCon/200,yCon/200,DLamD_f(:,:,9), level);
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-b');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['a)'],'interpreter','latex','FontSize',15);
%     %1.4955
%     
%     ax2=subplot(2,4,2);
%     dLamL1 = reshape(dLamH3_1(:,i,:),size(dLamH3_1,1),size(dLamH3_1,3));
%     polarcont(R,Theta_f1,dLamL1,8);
% %     caxis([min(min(min(dLamH3_1))),max(max(max(dLamH3_1)))]);
%     hold on
%     dLamLr1 = reshape(dLamH3r_1(:,i,:),size(dLamH3r_1,1),size(dLamH3r_1,3));
%     polarcont(R,Theta_fr1,dLamLr1,8);
% %     caxis([min(min(min(dLamH3r_1))),max(max(max(dLamH3r_1)))]);
%     polar(repmat(Theta_f1.',2,1),repmat(holder,1,length(Theta_f1)),'--r');
%     polar(repmat(Theta_f1(2),2,1),holder,'-b')
%     title(['b)'],'interpreter','latex','FontSize',15);
%     
%     ax3 = subplot(2,4,3);
%     fdLamRzr1 = cat(2, squeeze(fdLamRz1(:,:,i)),squeeze(fdLamRrr1(:,:,i)));
%     polarcont(R,Theta_final_1,fdLamRzr1,level);
% %     caxis([fmin_1,fmax_1]);
%     hold on
%     polar(repmat(Theta_f1.',2,1),repmat(holder,1,length(Theta_f1)),'--r');
%     polar(repmat(Theta_f1(2),2,1),holder,'-b');
%     title(['c)'],'interpreter','latex','FontSize',15);
% 
%     ax4 = subplot(2,4,4);
%     ddLamRzr1 = cat(2, squeeze(ddLamRz1(:,:,i)),squeeze(ddLamRrr1(:,:,i)));
%     polarcont(R,Theta_final_1t,ddLamRzr1,level);
% %     caxis([dmin_1,dmax_1]);
%     hold on 
%     polar(repmat(Theta_f1.',2,1),repmat(holder,1,length(Theta_f1)),'--r');
%     polar(repmat(Theta_f1(2),2,1),holder,'-b');
%     title(['d)'],'interpreter','latex','FontSize',15);

%     ax5=subplot(1,5,1);
%     contourf(xCon/200,-yCon/200,DLamC_f(:,:,4), level);
%     location = max(max(xCon/200))/size(DLamF1,2)*i;
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-b');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['a)'],'interpreter','latex','FontSize',15);
%     %5.041
%     

%     ax4=subplot(2,8,[1,2]);
% %     DLamC_fmax4 = max(max(DLamC_f(:,:,4)));
% %     DLamC_fmin4 = min(min(DLamC_f(:,:,4)));
%     DLamC4 = DLamC_f(:,:,17);
% %     DLamC4(per*DLamC_fmax4>DLamC4 & DLamC4>per*DLamC_fmin4) = 0;
%     contourf(xCon/200,-yCon/200,DLamC4, level);
%     location = xCon(1,i)/200;
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-b');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['a) 3.478'],'interpreter','latex','FontSize',15);    
%     caxis([-20*10^-3 20*10^-3]);
    
    ax5=subplot(2,8,[1,2,9,10]);
%     DLamC_fmax5 = max(max(DLamC_f(:,:,1)));
%     DLamC_fmin5 = min(min(DLamC_f(:,:,1)));
    DLamC5 = DLamC_f(:,:,3);
%     DLamC5(per*DLamC_fmax5>DLamC5 & DLamC5>per*DLamC_fmin5) = 0;
    contourf(xCon/200, -yCon/200, DLamC5, level);
    hold on
    plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-g');
    axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
      title(['a) 4.3378'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    caxis([-20*10^-3 20*10^-3]);

%     ax6=subplot(2,8,[9,10]);
% %     DLamC_fmax10 = max(max(DLamC_f(:,:,10)));
% %     DLamC_fmin10 = min(min(DLamC_f(:,:,10)));
%     DLamC10 = DLamC_f(:,:,5);
% %     DLamC10(per*DLamC_fmax10>DLamC10 & DLamC10>per*DLamC_fmin10) = 0;
%     contourf(xCon/200, -yCon/200, DLamC10, level);
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-k');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['c) 4.5781'],'interpreter','latex','FontSize',15);  
%     caxis([-20*10^-3 20*10^-3]);

    ax7=subplot(2,8,[3,4,11,12]);
%     DLamC_fmax19 = max(max(DLamC_f(:,:,6)));
%     DLamC_fmin19 = min(min(DLamC_f(:,:,6)));
    DLamC19 = DLamC_f(:,:,13);
%     DLamC19(per*DLamC_fmax19>DLamC19 & DLamC19>per*DLamC_fmin19) = 0;
    contourf(xCon/200, -yCon/200, DLamC19, level);
    hold on
    plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-m');
    axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
    title(['b) 4.7745'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    caxis([-20*10^-3 20*10^-3]);

    ax8=subplot(2,8,[5,6,13,14]);
    dLamL2 = reshape(dLamH3_2(:,i,:),size(dLamH3_2,1),size(dLamH3_2,3));
    if mean(mean(DLamL2)) == 0 
        level = 0;
    else
        level = 2;
    end
    polarcont(R,Theta_f2,dLamL2,level);
%     caxis([min(min(min(dLamH3_2))),max(max(max(dLamH3_2)))]);
    hold on
    dLamLr2 = reshape(dLamH3r_2(:,i,:),size(dLamH3r_2,1),size(dLamH3r_2,3));
    if mean(mean(DLamLr2)) == 0
        level = 0;
    else
        level = 2;
    end
    polarcont(R,Theta_fr2,dLamLr2,level);
%     caxis([min(min(min(dLamH3r_2))),max(max(max(dLamH3r_2)))]);
    polar(repmat(Theta_f2.',2,1),repmat(holder,1,length(Theta_f2)),'--r');
    polar(repmat(Theta_fr2(1),2,1),holder,'-b');
    polar(repmat(Theta_fr2(3),2,1),holder,'-g');
    polar(repmat(Theta_fr2(4),2,1),holder,'-k');
    polar(repmat(Theta_fr2(6),2,1),holder,'-m');
    caxis([-20*10^-3 20*10^-3]);
    title(['c)'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    
%     ax9 = subplot(2,8,[5,6,13,14]);
%     fdLamRzr2 = cat(2, squeeze(fdLamRz2(:,:,i)),squeeze(fdLamRrr2(:,:,i)));
%     polarcont(R,Theta_final_2,fdLamRzr2,level);
% %     caxis([fmin_2,fmax_2]);
%     hold on
%     polar(repmat(Theta_f2.',2,1),repmat(holder,1,length(Theta_f2)),'--r');
%     polar(repmat(Theta_fr2(13),2,1),holder,'-b')
%     polar(repmat(Theta_fr2(9),2,1),holder,'-g');
%     polar(repmat(Theta_fr2(5),2,1),holder,'-k');
%     polar(repmat(Theta_fr2(2),2,1),holder,'-m');
%     caxis([-20*10^-3 20*10^-3]);
%     title(['f)'],'interpreter','latex','FontSize',15);
    
    ax10 = subplot(2,8,[7,8,15,16]);
    ddLamRzr2 = cat(2, squeeze(ddLamRz2(:,:,i)),squeeze(ddLamRrr2(:,:,i)));
    if mean(mean(ddLamRzr2)) == 0
        level = 0;
    else
        level = 2;
    end    
    polarcont(R,Theta_final_2t,ddLamRzr2,level);
    hold on
    polar(repmat(Theta_f2.',2,1),repmat(holder,1,length(Theta_f2)),'--r');
    polar(repmat(Theta_fr2(1),2,1),holder,'-b');
    polar(repmat(Theta_fr2(3),2,1),holder,'-g');
    polar(repmat(Theta_fr2(4),2,1),holder,'-k');
    polar(repmat(Theta_fr2(6),2,1),holder,'-m');
    caxis([-20*10^-3 20*10^-3]);
    title(['d)'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    
    colormap(gray);
    saveas(fig, fullfile('/home/kevin/Desktop/polar_0.2_0.5_z//', ['polarcont_comb',num2str(i),'.jpg']));
end

clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax 11
clear DLamC5 DLamC_fmax5 DLamC_fmin5
clear DLamC4 DLamC_fmax4 DLamC_fmin4
clear DLamC19 DLamC_fmax19 DLamC_fmin19
clear DLamC10 DLamC_fmax10 DLamC_fmin10

%% 3D
disp('3D');
fig = figure('visible','off')
set(gcf,'pos',[10, 10, 800*2, 2*300]);


fdLamH3_2 = zeros(size(dLamH3_2,1),size(dLamH3_2,3),size(dLamH3_2,2));
fdLamH3r_2 = zeros(size(dLamH3r_2,1),size(dLamH3r_2,3),size(dLamH3r_2,2));

for i = 1:size(dLamH3_2,2)
    fdLamH3_2(:,:,i) = squeeze(dLamH3_2(:,i,:));
end
for i = 1:size(dLamH3r_2,2)
    fdLamH3r_2(:,:,i) = squeeze(dLamH3r_2(:,i,:));
end

ImfdLamH3_2 = zeros(size(ImdLamH3_2,1),size(ImdLamH3_2,3),size(ImdLamH3_2,2));
ImfdLamH3r_2 = zeros(size(ImdLamH3r_2,1),size(ImdLamH3r_2,3),size(ImdLamH3r_2,2));

for i = 1:size(ImdLamH3_2,2)
    ImfdLamH3_2(:,:,i) = squeeze(ImdLamH3_2(:,i,:));
    ImfdLam3_2 = ImfdLamH3_2(:,:,i);
    ImfdLam3_2(ImfdLam3_2>0) = 1;
    ImfdLam3_2(ImfdLam3_2<0) = -1;
    ImfdLamH3_2(:,:,i) = ImfdLam3_2;
end
for i = 1:size(ImdLamH3r_2,2)
    ImfdLamH3r_2(:,:,i) = squeeze(ImdLamH3r_2(:,i,:));
    ImfdLam3r_2 = ImfdLamH3r_2(:,:,i);
    ImfdLam3r_2(ImfdLam3r_2>0) = 1;
    ImfdLam3r_2(ImfdLam3r_2<0) = -1;
    ImfdLamH3r_2(:,:,i) = ImfdLam3r_2;
end

a = size(fdLamH3_2,1);
b = size(fdLamH3_2,2);

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = R(j)*cos(Theta_f2(k));
        y(j,k) = R(j)*sin(Theta_f2(k));
    end
end

x = repmat(x,1,1,size(fdLamH3_2,3));
y = repmat(y,1,1,size(fdLamH3_2,3));
z = zeros(size(fdLamH3_2));
for i =1:size(fdLamH3_2,3)
    z(:,:,i) = i;
end
    
fdLamH3_2 = squeeze(fdLamH3_2);
fdLamH3_2 = smooth3(fdLamH3_2);
pfdLamH3_2 =(fdLamH3_2>0);
pfdLamH3_2 = pfdLamH3_2.*fdLamH3_2; 
nfdLamH3_2 =(fdLamH3_2<0);
nfdLamH3_2 = nfdLamH3_2.*fdLamH3_2;

ax1 = subplot(2,8,[1,2,3]);
[fo, vo] = isosurface(x,y,z,pfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',fo,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');



%circle  
a = size(fdLamH3r_2,1);
b = size(fdLamH3r_2,2);

xr = zeros(a,b);
yr = zeros(a,b);

for j = 1:a
    for k = 1:b
        xr(j,k) = R(j)*cos(Theta_fr2(k));
        yr(j,k) = R(j)*sin(Theta_fr2(k));
    end
end

xr = repmat(xr,1,1,size(fdLamH3r_2,3));
yr = repmat(yr,1,1,size(fdLamH3r_2,3));
z = zeros(size(fdLamH3r_2));
for i =1:size(fdLamH3r_2,3)
    z(:,:,i) = i;
end
    

fdLamH3r_2 = squeeze(fdLamH3r_2);
fdLamH3r_2 = smooth3(fdLamH3r_2);

%recall Remainder matrix is that of flipped, thus <0 instead of >0
pfdLamH3r_2 =(fdLamH3r_2<0);
pfdLamH3r_2 = pfdLamH3r_2.*fdLamH3r_2; 
nfdLamH3r_2 =(fdLamH3r_2>0);
nfdLamH3r_2 = nfdLamH3r_2.*fdLamH3r_2;

[f1, v1] = isosurface(xr,yr,z,pfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nfdLamH3r_2);
% [f6, v6] = reducepatch(f1,v1,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');

C = [0,0, -1];
Theta = 0:0.01:2*pi;
X = C(1)+ 19.65*cos(Theta);
Y = C(2)+ 19.65*sin(Theta);
Z = C(3)+zeros(size(X));
plot3(X,Y,Z);

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
    title(['a) Linear int: Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
ax10 = subplot(2,8,4);
[fo, vo] = isosurface(x,y,z,pfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',fo,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');
[f1, v1] = isosurface(xr,yr,z,pfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nfdLamH3r_2);
% [f6, v6] = reducepatch(f1,v1,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
view(45,90);


clear f1 v1 pr2 fo vo p2 pr6 f6 v6 p5 f5 v5

%
%
ax2= subplot(2,8,[9,10,11]);
ImfdLamH3_2 = squeeze(ImfdLamH3_2);
ImfdLamH3_2 = smooth3(ImfdLamH3_2);
pImfdLamH3_2 =(ImfdLamH3_2>0);
pImfdLamH3_2 = pImfdLamH3_2.*ImfdLamH3_2; 
nImfdLamH3_2 =(ImfdLamH3_2<0);
nImfdLamH3_2 = nImfdLamH3_2.*ImfdLamH3_2; 

[fo, vo] = isosurface(x,y,z,pImfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nImfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

ImfdLamH3r_2 = squeeze(ImfdLamH3r_2);
ImfdLamH3r_2 = smooth3(ImfdLamH3r_2);
pImfdLamH3r_2 =(ImfdLamH3r_2<0);
pImfdLamH3r_2 = pImfdLamH3r_2.*ImfdLamH3r_2; 
nImfdLamH3r_2 =(ImfdLamH3r_2>0);
nImfdLamH3r_2 = nImfdLamH3r_2.*ImfdLamH3r_2; 


[f1, v1] = isosurface(xr,yr,z,pImfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImfdLamH3r_2);
% [f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
    title(['b) Linear int: Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
%
ax11 = subplot(2,8,12);
[fo, vo] = isosurface(x,y,z,pImfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nImfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

[f1, v1] = isosurface(xr,yr,z,pImfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImfdLamH3r_2);
% [f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
view(45,90);

clear f1 v1 pr2 fo vo p2 p3 pr3 f2 v2 f3 v3 pr6 f6 v6 p5 f5 v5

% %
% ax3 = subplot(1,6,3);
% [fo, vo] = isosurface(x,y,z,pfdLamH3_2);
% p2 = patch('Faces',fo,'Vertices',vo);
% set(p2, 'FaceColor','red','EdgeColor','none');
% hold all
% [f5, v5] = isosurface(x,y,z,nfdLamH3_2);
% p5 = patch('Faces',f5,'Vertices',v5);
% set(p5, 'FaceColor','blue','EdgeColor','none');
% 
% [f1, v1] = isosurface(xr,yr,z,pfdLamH3r_2);
% pr2 = patch('Faces',f1,'Vertices',v1);
% set(pr2, 'FaceColor','red','EdgeColor','none');
% 
% [f6, v6] = isosurface(xr,yr,z,nfdLamH3r_2);
% pr6 = patch('Faces',f6,'Vertices',v6);
% set(pr6, 'FaceColor','blue','EdgeColor','none');
% 
% [f2, v2] = isosurface(x,y,z,pImfdLamH3_2);
% p3 = patch('Faces',f2,'Vertices',v2);
% set(p3, 'FaceColor','red','EdgeColor','none');
% 
% [f7, v7] = isosurface(x,y,z,nImfdLamH3_2);
% p7 = patch('Faces',f7,'Vertices',v7);
% set(p7, 'FaceColor','blue','EdgeColor','none');
% 
% [f3, v3] = isosurface(xr,yr,z,pImfdLamH3r_2);
% pr3 = patch('Faces',f3,'Vertices',v3);
% set(pr3, 'FaceColor','red','EdgeColor','none');
% 
% [f8, v8] = isosurface(xr,yr,z,nImfdLamH3r_2);
% pr8 = patch('Faces',f8,'Vertices',v8);
% set(pr8, 'FaceColor','blue','EdgeColor','none');
% plot3(X,Y,Z);
% 
% view(3);
% axis tight
% daspect([1,1,1])
% camlight 
% lighting gouraud
% title(['c) No int'],'interpreter','latex','FontSize',15);
% 
% clear f1 v1 pr2 fo vo p2 p3 pr3 f2 v2 f3 v3 f5 f6 f7 f8 p5 p6 p7 p8
% clear pr8 p7 pr6 p5 
%

%
%
a = size(ddLamRz2,1);
b = size(ddLamRz2,2);

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = R(j)*cos(Theta_final2t(k));
        y(j,k) = R(j)*sin(Theta_final2t(k));
    end
end

x = repmat(x,1,1,size(ddLamRz2,3));
y = repmat(y,1,1,size(ddLamRz2,3));
z = zeros(size(ddLamRz2));
for i =1:size(ddLamRz2,3)
    z(:,:,i) = i;
end
    
ddLamRz2 = squeeze(ddLamRz2);
ddLamRz2 = smooth3(ddLamRz2);
pddLamRz2 =(ddLamRz2>0);
pddLamRz2 = pddLamRz2.*ddLamRz2; 
nddLamRz2 =(ddLamRz2<0);
nddLamRz2 = nddLamRz2.*ddLamRz2; 


ax4 = subplot(2,8,[5,6,7]);
[fo, vo] = isosurface(x,y,z,pddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

a = size(ddLamRrr2,1);
b = size(ddLamRrr2,2);

xr = zeros(a,b);
yr = zeros(a,b);

for j = 1:a
    for k = 1:b
        xr(j,k) = R(j)*cos(Theta_finalr2t(k));
        yr(j,k) = R(j)*sin(Theta_finalr2t(k));
    end
end

xr = repmat(xr,1,1,size(ddLamRrr2,3));
yr = repmat(yr,1,1,size(ddLamRrr2,3));
z = zeros(size(ddLamRrr2));
for i =1:size(ddLamRrr2,3)
    z(:,:,i) = i;
end
    

ddLamRrr2 = squeeze(ddLamRrr2);
ddLamRrr2 = smooth3(ddLamRrr2);
pddLamRrr2 =(ddLamRrr2<0);
pddLamRrr2 = pddLamRrr2.*ddLamRrr2; 
nddLamRrr2 =(ddLamRrr2>0);
nddLamRrr2 = nddLamRrr2.*ddLamRrr2; 


[f1, v1] = isosurface(xr,yr,z,pddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');

C = [0,0, -1];
Theta = 0:0.01:2*pi;
X = C(1)+ 19.65*cos(Theta);
Y = C(2)+ 19.65*sin(Theta);
Z = C(3)+zeros(size(X));
plot3(X,Y,Z);

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
title(['c)',num2str(plevel),'-poly int: Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
ax12 = subplot(2,8,8);
[fo, vo] = isosurface(x,y,z,pddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

[f1, v1] = isosurface(xr,yr,z,pddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

axis tight
daspect([1,1,1])
camlight 
camlight right
% camlight left
lighting gouraud

view(45,90);
clear f1 v1 pr2 fo vo p2 pr6 p5 f5 v5 f6 v6 

%
%
ax5= subplot(2,8,[13,14,15]);
ImddLamRz2 = squeeze(ImddLamRz2);
ImddLamRz2 = smooth3(ImddLamRz2);
pImddLamRz2 =(ImddLamRz2>0);
pImddLamRz2 = pImddLamRz2.*ImddLamRz2; 
nImddLamRz2 =(ImddLamRz2<0);
nImddLamRz2 = nImddLamRz2.*ImddLamRz2; 

[fo, vo] = isosurface(x,y,z,pImddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nImddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

ImddLamRrr2 = squeeze(ImddLamRrr2);
ImddLamRrr2 = smooth3(ImddLamRrr2);
pImddLamRrr2 =(ImddLamRrr2<0);
pImddLamRrr2 = pImddLamRrr2.*ImddLamRrr2; 
nImddLamRrr2 =(ImddLamRrr2>0);
nImddLamRrr2 = nImddLamRrr2.*ImddLamRrr2; 

[f1, v1] = isosurface(xr,yr,z,pImddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');

plot3(X,Y,Z);
view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
title(['d)',num2str(plevel),'-poly int: Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
ax14 = subplot(2,8,16);
[fo, vo] = isosurface(x,y,z,pImddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nImddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

[f1, v1] = isosurface(xr,yr,z,pImddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud

view(45,90);

%
%
% ax6 = subplot(1,6,6);
% [fo, vo] = isosurface(x,y,z,ddLamRz2);
% [fo, vo] = reducepatch(fo,vo,.5);
% p2 = patch('Faces',fo,'Vertices',vo);
% set(p2, 'FaceColor','blue','EdgeColor','none');
% hold all
% [f1, v1] = isosurface(xr,yr,z,ddLamRrr2);
% [f1, v1] = reducepatch(f1,v1,.5);
% pr2 = patch('Faces',f1,'Vertices',v1);
% set(pr2, 'FaceColor','red','EdgeColor','none');
% [f2, v2] = isosurface(x,y,z,ImddLamRz2);
% [f2, v2] = reducepatch(f2,v2,.5);
% p3 = patch('Faces',f2,'Vertices',v2);
% set(p3, 'FaceColor','green','EdgeColor','none');
% [f3, v3] = isosurface(xr,yr,z,ImddLamRrr2);
% [f3, v3] = reducepatch(f3,v3,.5);
% pr3 = patch('Faces',f3,'Vertices',v3);
% set(pr3, 'FaceColor','green','EdgeColor','none');
% plot3(X,Y,Z);
% view(3);
% axis tight
% daspect([1,1,1])
% camlight 
% lighting gouraud
% title(['f)',num2str(plevel),'-poly int'],'interpreter','latex','FontSize',15);

campos(ax1,[-35, -45, -0]);
campos(ax2,[-35, -45, -0]);
campos(ax4,[-35, -45, -0]);
campos(ax5,[-35, -45, -0]);


camroll(ax1,90);
camroll(ax2,90);
camroll(ax4,90);
camroll(ax5,90);
% %%
% 
% saveas(fig, fullfile('/home/kevin/Desktop//', ['isosurface_f.fig']));
saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_3.jpg']));
% 
% view(ax1,[-90,0]);
% view(ax2,[-90,0]);
% % view(ax3,[-90,0]);
% view(ax4,[-90,0]);
% view(ax5,[-90,0]);
% % view(ax6,[-90,0]);
% 
% 
% saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_5.jpg']));
% view(ax1,[-60,25]);
% view(ax2,[-60,25]);
% % view(ax3,[-60,25]);
% view(ax4,[-60,25]);
% view(ax5,[-60,25]);
% % view(ax6,[-60,25]);
% 
% 
% saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_6.jpg']));
% 
% view(ax1,[-60,90]);
% view(ax2,[-60,90]);
% % view(ax3,[-60,90]);
% view(ax4,[-60,90]);
% view(ax5,[-60,90]);
% % view(ax6,[-60,90]);
% 
% saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_7.jpg']));
% 
% %
%
%
%}
%
%%     contour(xCon/200,yCon/200,DLamC(:,:,k),3)%% OMD -mode calculation function
function [xCon, yCon, St, DLamSt] = slicedata(omd_rank,St_chosen,iptx_start,iptx_end,snapshotsV_temp,visz,visy)
snapshotsV_f = snapshotsV_temp(:,iptx_start:iptx_end);
setlength = size(snapshotsV_f,2);
if setlength > 500
    setlength = 500;
end

%limiting the number of setlenght (snapshot) lengths in order to 
%decrease the computational time, REMOVE this later for higher quality 
%flow data (might increases the accuracy)

visv = snapshotsV_f(:,1:setlength);
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
    interpolated_visv(i,:)= X(i,:);
end


%prepare data into form for contour plot
p = length(xdummy);
n = length(ydummy);
xCon = reshape(xnew,[n,p]);
yCon = reshape(ynew,[n,p]);

B = interpolated_visv(:,1:end-1);
A = interpolated_visv(:,2:end);

[L,M] = omd(A,B,omd_rank,[],'alternating'); 
[Y,D_nu]=eig(M);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd

Vand=zeros(size(M));
for loop1=1:size(B,2)
    Vand(:,loop1)=diag(D_nu).^(loop1-1); %creating vandermonde with increasing powers
end
U=L;
Psi = U*Y;

Z(:,1) =(Y^-1)*U'*B(:,1);

for i = 1:size(Vand,1)
    for j = 1:size(Vand,1)
        alpha_omd(j,i) = Z(j,1)*Vand(j,i);
    end
end

% 
alpha_omd_abs_I = abs(real(alpha_omd));
alpha_omd_abs_I_t = alpha_omd_abs_I.';

I = sum(alpha_omd_abs_I_t);
[C,T] = sort(I,'descend');

dt = 1/720;
freq = 0.1965/15/(2*pi)*imag(log(diag(D_nu))/dt);


%%%% closes to 0.2
[c St] = sort(abs(St_chosen-freq));
k = St(1);
St = freq(St(1));
% disp(St);
DLamSt = reshape(Psi(:,k),[n,p]);

% 
%%%% I method
% dummy2 = [c, St02,alpha_omd];
% dummy2 = dummy2(dummy2(:,1)<0.19,:);
% 
% sum_alpha = sum(dummy2(:,3:end),2);
% [cdummy St02dummy] = sort(sum_alpha);
% disp(St02dummy(1));
% St02 = dummy2(St02dummy(1),2);

% [d St006] = sort(abs(0.06-freq));

% dummy6 = [d, St006,alpha_omd];
% dummy6 = dummy6(dummy6(:,1)<0.059,:);
% 
% sum_alpha = sum(dummy6(:,3:end),2);
% [cdummy St006dummy] = sort(sum_alpha);
% St006 = dummy6(St006dummy(1),2);


%calculating DLam (DMD modes for St02 and St006)
% DLamSt02 = reshape(Psi(:,k),[n,p]);
% 
% i = St006;
% DLamSt006 = reshape(Psi(:,i),[n,p]);
%     


end

%% calc_cop_t function
function [ R, Theta] = calc_cop_t( p )

% calc_cop_t Function to give the coordinates of the COP for pressure in mm
% and radians
% timeseries

[np,nt]=size(p);

% set up the vector of radius values
rr=linspace(11,88,8);
r=repmat(rr',[8 1]);
th=[6*pi/4 7*pi/4 0 pi/4 pi/2 3*pi/4 pi 5*pi/4];
thet=reshape(repmat(th,[8 1]),64,1);

% Set up the vector of areas corresponding to each tapping
area = zeros(64,1);
for it=1:length(r)
    if r(it) == 11
        ro = mean([r(it+1),r(it)]);
        area(it) = pi*ro^2/8;
    elseif r(it) == 88
        ro = 98.25; ri = mean([r(it-1),r(it)]);
        area(it) = pi*(ro^2 - ri^2)/8;
    else
        ro = mean([r(it+1),r(it)]); ri = mean([r(it-1),r(it)]);
        area(it) = pi*(ro^2 - ri^2)/8;
    end
end
sum(area);
area=repmat(area,[1 nt]);
size(area);


% Evaluate the vector of forces, the total force and the terms FRcos(theta)
% and FRsin(theta)
f = p.*area;
F = sum(f);
size(F);
FRcos = (f'*(r.*cos(thet)))';
FRsin = (f'*(r.*sin(thet)))';

size(FRcos);
% Evaluate Theta and R
% Evaluate Theta in the range 0-pi for a continuous time series

% xcop=FRcos./F;
% ycop=FRsin./F;
% R = sqrt((FRcos.^2+FRsin.^2)./F.^2);
% Theta = acos(FRcos./(F.*R));
% Non-dimensionalise R
Theta = atan2(FRsin,FRcos);
R = abs(FRcos./cos(Theta)./F);
end


%% function for polar contour
function [C,h] = polarcont(r,theta,z,N,s)

[a,b] = size(z);

if a ~= length(r)
    error('r is not the same length as the first dimension of z')
end

if b ~= length(theta)
    error('theta is not the same length as the second dimension of z')
end

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = r(j)*cos(theta(k));
        y(j,k) = r(j)*sin(theta(k));
    end
end

if nargin == 3
    [C,h] = contourf(x,y,z);
elseif nargin == 4
    [C,h] = contourf(x,y,z,N);
elseif nargin == 5
    [C,h] = contourf(x,y,z,N,s);
else
    error('Incorrect number of inputs')
end 
end

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
%
%
%}
kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end


for i =3%1:length(iptx_f)/2
    cx1 = subplot(1,4,3)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['Theta:',num2str(Theta(i)),' c) Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(1,4,4)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['d) Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end



hold off
clear cx1 cx2 cx3 cx4 cx5 cx6 cx7 cx8




%%
for i =14%1:length(iptx_f)/2
    fig = figure('Visible','on');
    set(gcf,'pos',[10,10,3600,800]);
   
    hold all 
    cx1 = subplot(2,4,1)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['a) Real, level 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(2,4,5)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['e) Imag, level = 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % lower contour level
    cx3 = subplot(2,4,2)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level);
    title(['b) Real, level 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx4 = subplot(2,4,6)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level);
    title(['f) Imag, level = 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % filter data 
    cx5 = subplot(2,4,3)
    contourf(xCon/200,yCon/200,DLamC(:,:,i), level);
    title(['c) Real, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%       ylim([0.1,inf]);    
    caxis([-20*10^-3 20*10^-3]);


    cx6 = subplot(2,4,7)
    contourf(xCon/200,yCon/200,DLamD(:,:,i),level);
    title(['g) Imag, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    colormap gray
    caxis([-20*10^-3 20*10^-3]);

    
    % corrected filter data
    cx7 = subplot(2,4,4)
    contourf(xCon/200,yCon/200,DLamC_f(:,:,i), level);
    title(['d) Real, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);

    cx8 = subplot(2,4,8)
    contourf(xCon/200,yCon/200,DLamD_f(:,:,i),level);
    title(['h) Imag, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);


%       if Theta(i)>pi
%           ylim([cx5,cx6,cx7,cx8],[-inf, 0]);
%       else
%           ylim([cx5,cx6,cx7,cx8],[0, inf]);
%       end


    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end
hold off
clear cx1 cx2 cx3 cx4 cx5 cx6 cx7 cx8

%%
k =1;
fig = figure('Visible','on');
set(gcf,'pos',[10,10,1000,150]);

for i =13;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    colormap gray

end

k=2; 

for i =4;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    

    colormap gray

end


%{
%% slice images
disp('slices');
for i =14%1:length(iptx_f)/2
    fig = figure('Visible','on');
    set(gcf,'pos',[10,10,3600,800]);
   
    hold all 
    cx1 = subplot(2,4,1)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['a) Real, level 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx2 = subplot(2,4,5)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level+6);
    title(['e) Imag, level = 8'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % lower contour level
    cx3 = subplot(2,4,2)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level);
    title(['b) Real, level 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    
    cx4 = subplot(2,4,6)
    contourf(xCon/200,yCon/200,imag(DLamSt(:,:,i)),level);
    title(['f) Imag, level = 2'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
    
    % filter data 
    cx5 = subplot(2,4,3)
    contourf(xCon/200,yCon/200,DLamC(:,:,i), level);
    title(['c) Real, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%       ylim([0.1,inf]);    
    caxis([-20*10^-3 20*10^-3]);


    cx6 = subplot(2,4,7)
    contourf(xCon/200,yCon/200,DLamD(:,:,i),level);
    title(['g) Imag, filter'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    colormap gray
    caxis([-20*10^-3 20*10^-3]);

    
    % corrected filter data
    cx7 = subplot(2,4,4)
    contourf(xCon/200,yCon/200,DLamC_f(:,:,i), level);
    title(['d) Real, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);

    cx8 = subplot(2,4,8)
    contourf(xCon/200,yCon/200,DLamD_f(:,:,i),level);
    title(['h) Imag, correction'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);
%     ylim([0.1,inf]);
    caxis([-20*10^-3 20*10^-3]);


%       if Theta(i)>pi
%           ylim([cx5,cx6,cx7,cx8],[-inf, 0]);
%       else
%           ylim([cx5,cx6,cx7,cx8],[0, inf]);
%       end


    colormap gray
%     saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd_z/', ['Theta:',num2str(Theta(i)),' St=',num2str(St(i)),'.jpeg']));
   
end
hold off
clear cx1 cx2 cx3 cx4 cx5 cx6 cx7 cx8

%%
k =1;
fig = figure('Visible','on');
set(gcf,'pos',[10,10,1000,300]);

for i =13;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    
    colormap gray

end

k=2; 

for i =4;

    hold all 
    cx1 = subplot(1,2,k)
    contourf(xCon/200,yCon/200,real(DLamSt(:,:,i)), level+6);
    title(['\theta:' num2str(Theta(i))],'fontweight','normal','FontSize',22,'fontname','arial');
    axis equal
    axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]);    

    colormap gray

end

%}
%}
%% image matching between real and imag
disp('image/color matching');
dummy1 = [St.',Theta.',[1:size(DLamSt,3)].'];
dummy1 = dummy1(dummy1(:,2)<2.3,:); %either >3.3 or <1.8 modes (switch ylim accordingly)
dummy1 = dummy1(dummy1(:,1)>.15 & dummy1(:,1)<.25,:);
dummy1 = sortrows(dummy1,2);
%0.03 0.08 for 0.06
%0.15 0.25 for 0.2
%0.3 0.5 for 0.4
dummy2 = [St.',Theta.',[1:size(DLamSt,3)].'];
dummy2 = dummy2(dummy2(:,2)>3.3,:); %either >3.3 or <1.8 modes (switch ylim accordingly)
dummy2 = dummy2(dummy2(:,1)>.15 & dummy2(:,1)<.25,:);
dummy2 = sortrows(dummy2,2);

DLamF1(:,:,1) = DLamC_f((end+1)/2:end,:,dummy1(1,3)); %starting mode
DLamG1(:,:,1) = DLamC_f(1:(end+1)/2,:,dummy1(1,3)); %starting mode


ImDLamF1(:,:,1) = DLamD_f((end+1)/2:end,:,dummy1(1,3)); %starting mode
ImDLamG1(:,:,1) = DLamD_f(1:(end+1)/2,:,dummy1(1,3)); %starting mode

% switch the focus to DLamF1 if need be and vice versa
for i = 2:size(dummy1,1)
    k = dummy1(i,3);    
    
        DLamCdummy = DLamC_f(1:((end+1)/2),:,k);
        DLamDdummy = DLamD_f(1:((end+1)/2),:,k);
    
    DLamCdummy2 = DLamCdummy/norm(DLamCdummy);
    DLamDdummy2 = DLamDdummy/norm(DLamDdummy);
    
    DLamFdummy = DLamG1(:,:,i-1)/norm(DLamG1(:,:,i-1));
    
    if sum(max(abs(DLamFdummy-DLamCdummy2)))<= sum(max(abs(DLamFdummy-DLamDdummy2)))
        DLamG1(:,:,i) = DLamCdummy;
        DLamF1(:,:,i) = DLamC_f((end+1)/2:end,:,k);
        ImDLamG1(:,:,i) = DLamDdummy;
        ImDLamG1(:,:,i) = DLamD_f((end+1)/2:end,:,k);
    else
        DLamG1(:,:,i) = DLamDdummy;
        DLamF1(:,:,i) = DLamD_f((end+1)/2:end,:,k);
        ImDLamG1(:,:,i) = DLamCdummy;
        ImDLamF1(:,:,i) = DLamC_f((end+1)/2:end,:,k);
    end
end
clear DLamCdummy DLamCdummy2 DLamDdummy DLamDdummy2 DLamFdummy


%
DLamF2(:,:,1) = DLamC_f((end+1)/2:end,:,dummy2(1,3)); %starting mode
DLamG2(:,:,1) = DLamC_f(1:(end+1)/2,:,dummy2(1,3)); %starting mode

ImDLamF2(:,:,1) = DLamD_f((end+1)/2:end,:,dummy2(1,3)); %starting mode
ImDLamG2(:,:,1) = DLamD_f(1:(end+1)/2,:,dummy2(1,3)); %starting mode

for i = 2:size(dummy2,1)
    k = dummy2(i,3);    
    
        DLamCdummy = DLamC_f(1:(end+1)/2,:,k);
        DLamDdummy = DLamD_f(1:(end+1)/2,:,k);
    
    DLamCdummy2 = DLamCdummy/norm(DLamCdummy);
    DLamDdummy2 = DLamDdummy/norm(DLamDdummy);
    
    DLamFdummy = DLamG2(:,:,i-1)/norm(DLamG2(:,:,i-1));
    
    if sum(max(abs(DLamFdummy-DLamCdummy2)))<= sum(max(abs(DLamFdummy-DLamDdummy2)))
        DLamG2(:,:,i) = DLamCdummy;
        DLamF2(:,:,i) = DLamC_f((end+1)/2:end,:,k);
        ImDLamG2(:,:,i) = DLamDdummy;
        ImDLamF2(:,:,i) = DLamD_f((end+1)/2:end,:,k);
    else
        DLamG2(:,:,i) = DLamDdummy;
        DLamF2(:,:,i) = DLamD_f((end+1)/2:end,:,k);
        ImDLamG2(:,:,i) = DLamCdummy;
        ImDLamF2(:,:,i) = DLamC_f((end+1)/2:end,:,k);
    end
end
clear DLamCdummy DLamCdummy2 DLamDdummy DLamDdummy2 DLamFdummy

clear DLamCdummy DLamCdummy2 DLamDdummy DLamDdummy2 DLamFdummy



%%
fig = figure('visible','on');
set(gcf,'pos',[10,10,900*size(DLamF1,3),800]);

DLamMix = zeros(size(DLamC_f,1),size(DLamC_f,2),size(DLamG2,3));
for i= 1:size(DLamG2,3)
    DLamMix(1:(end+1)/2,:,i) = DLamG2(:,:,i);
    DLamMix((end+1)/2:end,:,i) = DLamF2(:,:,i);
end

for i = 1:size(DLamG2,3)
    subplot(2,((size(DLamF2,3))+1)/2,i)
    contourf(xCon/200,yCon/200,DLamMix(:,:,i),level);
    title(['\theta:' num2str(dummy2(i,2))],'fontweight','normal','FontSize',22,'fontname','arial');
    caxis([-20*10^-3 20*10^-3]);

    colormap(gray)
end
% 
% %  saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd//', ['filtered.jpeg']));



%% contourslice stack
% xNew = xCon((end+1)/2:end,:);
% yNew = yCon((end+1)/2:end,:);
% 
% for i=1:size(dummy2,1)-1
%     fig(i) = figure('visible','on');
%     hold all
%     contour(xNew/200,yNew/200,DLamG2(:,:,i),1)
%     saveas(fig(i), fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['Stack: Theta:',num2str(dummy2(i,2)),' St=',num2str(dummy2(i),1),'.jpeg']));
%     
%     contour(xNew/200,yNew/200,DLamG2(:,:,i+1),1)
%     saveas(fig(i), fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['Stack: Theta:',num2str(dummy2(i,2)),' St=',num2str(dummy2(i),1),'.jpeg']));
% 
%     close all
% end
% 
% 
% fig = figure('visible','on');
% for i=1:size(dummy2,1)
%     hold all
%     contour(xNew/200,yNew/200,DLamG2(:,:,i),1)
% end
% saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['ALL Stack: Theta:',num2str(dummy2(1,2)),'~',num2str(dummy2(end,2)),'.jpeg']));
% 
% fig = figure('visible','on');
% 
% for i=1:size(dummy2,1)
%     hold all
%     contour(xNew/200,yNew/200,DLamG2(:,:,i),1)
% end
% saveas(fig, fullfile('/home/kevin/Desktop/006_saved_5omd/contour_stack/', ['180 ALL Stack: Theta:',num2str(dummy2(1,2)),'~',num2str(dummy2(end,2)),'.jpeg']));

%% setting up polar 
disp('polar');
N = size(DLamSt,3);

k = size(DLamSt,1);
R = zeros(k,1);
for i = 1:k
    R(i) = -19.65 +19.65*2/(k-1)*(i-1);
end 
R=R.';

%%
disp('no filter');

Theta_f1 = dummy1(:,2);
dummyG1 = [Theta_f1,[1:size(DLamF1,3)].'];
dummyG1 = sortrows(dummyG1,1);
R = R(R>=0);
Theta_f1 = dummyG1(:,1);

Theta_f2 = dummy2(:,2)-pi;
dummyG2 = [Theta_f2,[1:size(DLamF2,3)].'];
dummyG2 = sortrows(dummyG2,1);
R = R(R>=0);
Theta_f2 = dummyG2(:,1);

for i = 1:size(DLamF1,3)
    DLamL1(:,:,i) = DLamF1(:,:,dummyG1(i,2));
    ImDLamL1(:,:,i) = ImDLamF1(:,:,dummyG1(i,2));
end

for i = 1:size(DLamF2,3)
    DLamL2(:,:,i) = DLamF2(:,:,dummyG2(i,2));
    ImDLamL2(:,:,i) = ImDLamF2(:,:,dummyG2(i,2));

end

% for the other 180 degrees
Theta_fr1 = dummy1(:,2)+pi;
dummyGr1 = [Theta_fr1,[1:size(DLamG1,3)].'];
dummyGr1 = sortrows(dummyGr1,1);
Theta_fr1 = dummyGr1(:,1);

for i = 1:size(DLamG1,3)
    DLamLr1(:,:,i)=DLamG1(:,:,dummyGr1(i,2));
    ImDLamLr1(:,:,i)=ImDLamG1(:,:,dummyGr1(i,2));

end

Theta_fr2 = dummy2(:,2);
dummyGr2 = [Theta_fr2,[1:size(DLamG2,3)].'];
dummyGr2 = sortrows(dummyGr2,1);
Theta_fr2 = dummyGr2(:,1);

for i = 1:size(DLamG2,3)
    DLamLr2(:,:,i)=DLamG2(:,:,dummyGr2(i,2));
    ImDLamLr2(:,:,i)= ImDLamG2(:,:,dummyGr2(i,2));
end

% for i =1:size(DLamF1,2)
% %       not displaying figures
%     fig = figure('Visible','off');
%     set(gcf,'pos',[10,10,900,900]);
%     
%     hold all 
%     dLamL = reshape(DLamL1(:,i,:),size(DLamL1,1),size(DLamL1,3));
%     polarcont(R,Theta_f1,dLamL,level);    
%     colormap(gray)
%     saveas(fig, fullfile('/home/kevin/Desktop/meh/polar<3/<pi/', ['polarcont',num2str(i),'.jpeg']));
% 
% end
clear dLamF dLamG 
%% data shuffle (cart->polar
disp('cart->polar');
level = 4;
for i = 1:size(DLamF1,2)
    dLamL1 = reshape(DLamL1(:,i,:),size(DLamL1,1),size(DLamL1,3));
    dLamH3_1(:,i,:) = dLamL1;
    ImdLamL1 = reshape(ImDLamL1(:,i,:),size(ImDLamL1,1),size(ImDLamL1,3));
    ImdLamH3_1(:,i,:) = ImdLamL1;
    
end
for i = 1:size(DLamF2,2)
    dLamL2 = reshape(DLamL2(:,i,:),size(DLamL2,1),size(DLamL2,3));
    dLamH3_2(:,i,:) = dLamL2;
    ImdLamL2 = reshape(ImDLamL2(:,i,:),size(ImDLamL2,1),size(ImDLamL2,3));
    ImdLamH3_2(:,i,:) = ImdLamL2;
end
% for the other 180 degrees
for i = 1:size(DLamG1,2)
    dLamLr1 = reshape(DLamLr1(:,i,:),size(DLamLr1,1),size(DLamLr1,3));
    dLamH3r_1(:,i,:)=dLamLr1;
    ImdLamLr1 = reshape(ImDLamLr1(:,i,:),size(ImDLamLr1,1),size(ImDLamLr1,3));
    ImdLamH3r_1(:,i,:)=ImdLamLr1;
end
for i = 1:size(DLamG2,2)
    dLamLr2 = reshape(DLamLr2(:,i,:),size(DLamLr2,1),size(DLamLr2,3));
    dLamH3r_2(:,i,:)=dLamLr2;
    ImdLamLr2 = reshape(ImDLamLr2(:,i,:),size(ImDLamLr2,1),size(ImDLamLr2,3));
    ImdLamH3r_2(:,i,:)=ImdLamLr2;
end
close all
clear  dLamH dLamL1 dLamL2 dLamLr1 dLamLr2
clear  dLamH ImdLamL1 ImdLamL2 ImdLamLr1 ImdLamLr2

%% 3D - poly_int
disp('poly_int');

%plevel = polyfit level
plevel = 2;
buffer = 10/180*pi;

for i =1:size(DLamF1,2)
    dLamL1 = reshape(dLamH3_1(:,i,:),size(dLamH3_1,1),size(dLamH3_1,3));
    
    Theta_f_0 = linspace(-pi/4,min(Theta_f1)-buffer,20);
    
    Theta_ft = linspace(min(Theta_f1)-buffer,max(Theta_f1)+buffer,25*size(Theta_f1,1));
    
    ddLamG_r1 = zeros(length(dLamL1),25*size(Theta_f1,1));
    for k = 1:size(Theta_f1,1)
        [O,T] = sort(abs(Theta_ft-Theta_f1(k)));
        ddLamG_r1(:,T(1)) = dLamL1(:,k);
    end
    for l = 1:size(ddLamG_r1,1)
        for j = 1:size(Theta_f1)
            poly_r1(l,j+1)= dLamL1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f1)+2,1);
        Thetadummy(2:end-1) = Theta_f1.';
        Thetadummy(1) = min(Theta_f1)-buffer;
        Thetadummy(end) = max(Theta_f1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft);
        ddLamG_r1(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f1)+buffer,pi-pi/4,20);
    
    Theta_final1t = cat(2,Theta_f_0, Theta_ft,Theta_f_180);
    Filler = zeros(size(dLamL1,1),20);
    ddLamR1 = cat(2,Filler,ddLamG_r1,Filler);
    ddLamRz1(:,:,i) = ddLamR1;
end

clear Theta_ft ddLamG_r1 p1 y2 ddLamR1 poly_r1 Thetadummy

for i =1:size(DLamL1,2)   
    dLamLr1 = reshape(dLamH3r_1(:,i,:),size(dLamH3r_1,1),size(dLamH3r_1,3));

    
    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr1)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr1)-buffer,max(Theta_fr1)+buffer,25*size(Theta_fr1,1));
    
    ddLamG_rr1 = zeros(length(dLamLr1),25*size(Theta_fr1,1));
    for k = 1:size(Theta_fr1,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr1(k)));
        ddLamG_rr1(:,T(1)) = dLamLr1(:,k);
    end
    for l = 1:size(ddLamG_rr1,1)
        for j = 1:size(Theta_fr1)
            poly_r1(l,j+1)= dLamLr1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr1)+2,1);
        Thetadummy(2:end-1) = Theta_fr1.';
        Thetadummy(1) = min(Theta_fr1)-buffer;
        Thetadummy(end) = max(Theta_fr1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ddLamG_rr1(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr1)+buffer,2*pi-pi/4,20);
    
    Theta_finalr1t = cat(2,Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr1,1),20);    
    ddLamRr1 = cat(2,Filler, ddLamG_rr1,Filler);
    ddLamRrr1(:,:,i) = ddLamRr1;
    
end

clear ddLamG_rr1 p1 y2 poly_r1 ddLamRr1 Thetadummy 


 for i =1:size(DLamL2,2)

    dLamL2 = reshape(dLamH3_2(:,i,:),size(dLamH3_2,1),size(dLamH3_2,3));

    Theta_f_0 = linspace(-pi/4,min(Theta_f2)-buffer,20);
    
    Theta_ft2 = linspace(min(Theta_f2)-buffer,max(Theta_f2)+buffer,20*size(Theta_f2,1));
    
    ddLamG_r2 = zeros(length(dLamL2),20*size(Theta_f2,1));
    for k = 1:size(Theta_f2,1)
        [O,T] = sort(abs(Theta_ft2-Theta_f2(k)));
        ddLamG_r2(:,T(1)) = dLamL2(:,k);
    end
    
    for l = 1:size(ddLamG_r2,1)
        for j = 1:size(Theta_f2)
            poly_r1(l,j+1)= dLamL2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f2)+2,1);
        Thetadummy(2:end-1) = Theta_f2.';
        Thetadummy(1) = min(Theta_f2)-buffer;
        Thetadummy(end) = max(Theta_f2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft2);
        ddLamG_r2(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f2)+buffer,pi-pi/4,20);
 

    Theta_final2t = cat(2,Theta_f_0, Theta_ft2,Theta_f_180);
    Filler = zeros(size(dLamL2,1),20);
    ddLamR2 = cat(2,Filler,ddLamG_r2,Filler);
    ddLamRz2(:,:,i) = ddLamR2;
 end

 clear ddLamG_r2 y2 p1 ddLamR2 Thetadummy poly_r1

for i =1:size(DLamL2,2)
    dLamLr2 = reshape(dLamH3r_2(:,i,:),size(dLamH3r_2,1),size(dLamH3r_2,3));

    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr2)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr2)-buffer,max(Theta_fr2)+buffer,20*size(Theta_fr2,1));
    
    ddLamG_rr2 = zeros(length(dLamLr2),20*size(Theta_fr2,1));
    for k = 1:size(Theta_fr2,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr2(k)));
        ddLamG_rr2(:,T(1)) = dLamLr2(:,k);
    end
     for l = 1:size(ddLamG_rr2,1)
        for j = 1:size(Theta_fr2)
            poly_r1(l,j+1)= dLamLr2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_fr2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr2)+2,1);
        Thetadummy(2:end-1) = Theta_fr2.';
        Thetadummy(1) = min(Theta_fr2)-buffer;
        Thetadummy(end) = max(Theta_fr2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ddLamG_rr2(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr2),2*pi-pi/4,20);
    
    Theta_finalr2t = cat(2, Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr2,1),20);    
    ddLamRr2 = cat(2,Filler, ddLamG_rr2,Filler);
    ddLamRrr2(:,:,i) = ddLamRr2;    
end

clear p1 y2 ddLamG_r2 ddLamRr2
clear Theta_fr_2pi Theta_f_180 Theta_f_0 Theta_180_flr
clear poly_r1 O T Thetadummy

%% ^ Im

disp('Im poly_int');

%plevel = polyfit level
plevel = 2;
buffer = 10/180*pi;

for i =1:size(ImDLamF1,2)
    ImdLamL1 = reshape(ImdLamH3_1(:,i,:),size(ImdLamH3_1,1),size(ImdLamH3_1,3));
    
    Theta_f_0 = linspace(-pi/4,min(Theta_f1)-buffer,20);
    
    Theta_ft = linspace(min(Theta_f1)-buffer,max(Theta_f1)+buffer,25*size(Theta_f1,1));
    
    ImddLamG_r1 = zeros(length(ImdLamL1),25*size(Theta_f1,1));
    for k = 1:size(Theta_f1,1)
        [O,T] = sort(abs(Theta_ft-Theta_f1(k)));
        ImddLamG_r1(:,T(1)) = ImdLamL1(:,k);
    end
    for l = 1:size(ImddLamG_r1,1)
        for j = 1:size(Theta_f1)
            poly_r1(l,j+1)= ImdLamL1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f1)+2,1);
        Thetadummy(2:end-1) = Theta_f1.';
        Thetadummy(1) = min(Theta_f1)-buffer;
        Thetadummy(end) = max(Theta_f1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft);
        ImddLamG_r1(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f1)+buffer,pi-pi/4,20);
    
    Theta_final1t = cat(2,Theta_f_0, Theta_ft,Theta_f_180);
    Filler = zeros(size(dLamL1,1),20);
    ImddLamR1 = cat(2,Filler,ImddLamG_r1,Filler);
    ImddLamRz1(:,:,i) = ImddLamR1;
end

clear Theta_ft ImddLamG_r1 p1 y2 ImddLamR1 poly_r1 Thetadummy

for i =1:size(ImDLamL1,2)   
    ImdLamLr1 = reshape(ImdLamH3r_1(:,i,:),size(ImdLamH3r_1,1),size(ImdLamH3r_1,3));

    
    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr1)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr1)-buffer,max(Theta_fr1)+buffer,25*size(Theta_fr1,1));
    
    ImddLamG_rr1 = zeros(length(ImdLamLr1),25*size(Theta_fr1,1));
    for k = 1:size(Theta_fr1,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr1(k)));
        ImddLamG_rr1(:,T(1)) = dLamLr1(:,k);
    end
    for l = 1:size(ImddLamG_rr1,1)
        for j = 1:size(Theta_fr1)
            poly_r1(l,j+1)= ImdLamLr1(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f1)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr1)+2,1);
        Thetadummy(2:end-1) = Theta_fr1.';
        Thetadummy(1) = min(Theta_fr1)-buffer;
        Thetadummy(end) = max(Theta_fr1)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ImddLamG_rr1(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr1)+buffer,2*pi-pi/4,20);
    
    Theta_finalr1t = cat(2,Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr1,1),20);    
    ImddLamRr1 = cat(2,Filler, ImddLamG_rr1,Filler);
    ImddLamRrr1(:,:,i) = ImddLamRr1;
    
end

clear ImddLamG_rr1 p1 y2 poly_r1 ImddLamRr1 Thetadummy 


 for i =1:size(ImDLamL2,2)

    ImdLamL2 = reshape(ImdLamH3_2(:,i,:),size(ImdLamH3_2,1),size(ImdLamH3_2,3));

    Theta_f_0 = linspace(-pi/4,min(Theta_f2)-buffer,20);
    
    Theta_ft2 = linspace(min(Theta_f2)-buffer,max(Theta_f2)+buffer,20*size(Theta_f2,1));
    
    ImddLamG_r2 = zeros(length(ImdLamL2),20*size(Theta_f2,1));
    for k = 1:size(Theta_f2,1)
        [O,T] = sort(abs(Theta_ft2-Theta_f2(k)));
        ImddLamG_r2(:,T(1)) = ImdLamL2(:,k);
    end
    
    for l = 1:size(ImddLamG_r2,1)
        for j = 1:size(Theta_f2)
            poly_r1(l,j+1)= ImdLamL2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_f2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_f2)+2,1);
        Thetadummy(2:end-1) = Theta_f2.';
        Thetadummy(1) = min(Theta_f2)-buffer;
        Thetadummy(end) = max(Theta_f2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ft2);
        ImddLamG_r2(l,:)=y2;
    end
    
    Theta_f_180 = linspace(max(Theta_f2)+buffer,pi-pi/4,20);
 

    Theta_final2t = cat(2,Theta_f_0, Theta_ft2,Theta_f_180);
    Filler = zeros(size(dLamL2,1),20);
    ImddLamR2 = cat(2,Filler,ImddLamG_r2,Filler);
    ImddLamRz2(:,:,i) = ImddLamR2;
 end

 clear ImddLamG_r2 y2 p1 ImddLamR2 Thetadummy poly_r1

for i =1:size(ImDLamL2,2)
    ImdLamLr2 = reshape(ImdLamH3r_2(:,i,:),size(ImdLamH3r_2,1),size(ImdLamH3r_2,3));

    Theta_180_flr = linspace(pi-pi/4,min(Theta_fr2)-buffer,20);
    
    Theta_ftr = linspace(min(Theta_fr2)-buffer,max(Theta_fr2)+buffer,20*size(Theta_fr2,1));
    
    ImddLamG_rr2 = zeros(length(ImdLamLr2),20*size(Theta_fr2,1));
    for k = 1:size(Theta_fr2,1)
        [O,T] = sort(abs(Theta_ftr-Theta_fr2(k)));
        ImddLamG_rr2(:,T(1)) = ImdLamLr2(:,k);
    end
     for l = 1:size(ImddLamG_rr2,1)
        for j = 1:size(Theta_fr2)
            poly_r1(l,j+1)= ImdLamLr2(l,j);
        end
        poly_r1(:,1) = 0;
        poly_r1(:,length(Theta_fr2)+2) = 0;
        
        Thetadummy = zeros(length(Theta_fr2)+2,1);
        Thetadummy(2:end-1) = Theta_fr2.';
        Thetadummy(1) = min(Theta_fr2)-buffer;
        Thetadummy(end) = max(Theta_fr2)+buffer;
        
        p1 = polyfit(Thetadummy.',poly_r1(l,:),plevel);
        y2 = polyval(p1,Theta_ftr);
        ImddLamG_rr2(l,:)=y2;
    end
    
    Theta_fr_2pi = linspace(max(Theta_fr2),2*pi-pi/4,20);
    
    Theta_finalr2t = cat(2, Theta_180_flr, Theta_ftr,Theta_fr_2pi);
    Filler = zeros(size(dLamLr2,1),20);    
    ImddLamRr2 = cat(2,Filler, ImddLamG_rr2,Filler);
    ImddLamRrr2(:,:,i) = ImddLamRr2;    
end

clear p1 y2 ddLamG_r2 ddLamRr2
clear Theta_fr_2pi Theta_f_180 Theta_f_0 Theta_180_flr
clear poly_r1 O T Thetadummy

%{
%% polar slice disp
disp('filter&nofilter');

dmin_1 = min(min(min(cat(1,ddLamRz1,ddLamRrr1))));
dmax_1 = max(max(max(cat(1,ddLamRz1,ddLamRrr1))));
dmin_2 = min(min(min(cat(1,ddLamRz2,ddLamRrr2))));
dmax_2 = max(max(max(cat(1,ddLamRz2,ddLamRrr2))));

for i =100:10:180%:size(DLamF1,2)
%     %   not displaying figures
    level = 2;

    fig = figure('Visible','off');
    set(gcf,'pos',[10,10,2000,500]);
    holder = [-19.65;19.65];
    Theta_final_1t = cat(2,Theta_final1t,Theta_finalr1t);
    Theta_final_2t = cat(2,Theta_final2t,Theta_finalr2t);
    location = max(max(xCon/200))/size(DLamF1,2)*i;
%     ax1=subplot(2,4,1);
%     contourf(xCon/200,yCon/200,DLamD_f(:,:,9), level);
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-b');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['a)'],'interpreter','latex','FontSize',15);
%     %1.4955
%     
%     ax2=subplot(2,4,2);
%     dLamL1 = reshape(dLamH3_1(:,i,:),size(dLamH3_1,1),size(dLamH3_1,3));
%     polarcont(R,Theta_f1,dLamL1,8);
% %     caxis([min(min(min(dLamH3_1))),max(max(max(dLamH3_1)))]);
%     hold on
%     dLamLr1 = reshape(dLamH3r_1(:,i,:),size(dLamH3r_1,1),size(dLamH3r_1,3));
%     polarcont(R,Theta_fr1,dLamLr1,8);
% %     caxis([min(min(min(dLamH3r_1))),max(max(max(dLamH3r_1)))]);
%     polar(repmat(Theta_f1.',2,1),repmat(holder,1,length(Theta_f1)),'--r');
%     polar(repmat(Theta_f1(2),2,1),holder,'-b')
%     title(['b)'],'interpreter','latex','FontSize',15);
%     
%     ax3 = subplot(2,4,3);
%     fdLamRzr1 = cat(2, squeeze(fdLamRz1(:,:,i)),squeeze(fdLamRrr1(:,:,i)));
%     polarcont(R,Theta_final_1,fdLamRzr1,level);
% %     caxis([fmin_1,fmax_1]);
%     hold on
%     polar(repmat(Theta_f1.',2,1),repmat(holder,1,length(Theta_f1)),'--r');
%     polar(repmat(Theta_f1(2),2,1),holder,'-b');
%     title(['c)'],'interpreter','latex','FontSize',15);
% 
%     ax4 = subplot(2,4,4);
%     ddLamRzr1 = cat(2, squeeze(ddLamRz1(:,:,i)),squeeze(ddLamRrr1(:,:,i)));
%     polarcont(R,Theta_final_1t,ddLamRzr1,level);
% %     caxis([dmin_1,dmax_1]);
%     hold on 
%     polar(repmat(Theta_f1.',2,1),repmat(holder,1,length(Theta_f1)),'--r');
%     polar(repmat(Theta_f1(2),2,1),holder,'-b');
%     title(['d)'],'interpreter','latex','FontSize',15);

%     ax5=subplot(1,5,1);
%     contourf(xCon/200,-yCon/200,DLamC_f(:,:,4), level);
%     location = max(max(xCon/200))/size(DLamF1,2)*i;
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-b');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['a)'],'interpreter','latex','FontSize',15);
%     %5.041
%     

%     ax4=subplot(2,8,[1,2]);
% %     DLamC_fmax4 = max(max(DLamC_f(:,:,4)));
% %     DLamC_fmin4 = min(min(DLamC_f(:,:,4)));
%     DLamC4 = DLamC_f(:,:,17);
% %     DLamC4(per*DLamC_fmax4>DLamC4 & DLamC4>per*DLamC_fmin4) = 0;
%     contourf(xCon/200,-yCon/200,DLamC4, level);
%     location = xCon(1,i)/200;
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-b');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['a) 3.478'],'interpreter','latex','FontSize',15);    
%     caxis([-20*10^-3 20*10^-3]);
    
    ax5=subplot(2,8,[1,2,9,10]);
%     DLamC_fmax5 = max(max(DLamC_f(:,:,1)));
%     DLamC_fmin5 = min(min(DLamC_f(:,:,1)));
    DLamC5 = DLamC_f(:,:,3);
%     DLamC5(per*DLamC_fmax5>DLamC5 & DLamC5>per*DLamC_fmin5) = 0;
    contourf(xCon/200, -yCon/200, DLamC5, level);
    hold on
    plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-g');
    axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
      title(['a) 4.3378'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    caxis([-20*10^-3 20*10^-3]);

%     ax6=subplot(2,8,[9,10]);
% %     DLamC_fmax10 = max(max(DLamC_f(:,:,10)));
% %     DLamC_fmin10 = min(min(DLamC_f(:,:,10)));
%     DLamC10 = DLamC_f(:,:,5);
% %     DLamC10(per*DLamC_fmax10>DLamC10 & DLamC10>per*DLamC_fmin10) = 0;
%     contourf(xCon/200, -yCon/200, DLamC10, level);
%     hold on
%     plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-k');
%     axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
%     title(['c) 4.5781'],'interpreter','latex','FontSize',15);  
%     caxis([-20*10^-3 20*10^-3]);

    ax7=subplot(2,8,[3,4,11,12]);
%     DLamC_fmax19 = max(max(DLamC_f(:,:,6)));
%     DLamC_fmin19 = min(min(DLamC_f(:,:,6)));
    DLamC19 = DLamC_f(:,:,13);
%     DLamC19(per*DLamC_fmax19>DLamC19 & DLamC19>per*DLamC_fmin19) = 0;
    contourf(xCon/200, -yCon/200, DLamC19, level);
    hold on
    plot([location,location],[min(min(yCon/200)),max(max(yCon/200))],'-m');
    axis equal
%     axis([6.3352/200 346.6375/200 -147.8552/200 150.1393/200]); 
    title(['b) 4.7745'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    caxis([-20*10^-3 20*10^-3]);

    ax8=subplot(2,8,[5,6,13,14]);
    dLamL2 = reshape(dLamH3_2(:,i,:),size(dLamH3_2,1),size(dLamH3_2,3));
    if mean(mean(DLamL2)) == 0 
        level = 0;
    else
        level = 2;
    end
    polarcont(R,Theta_f2,dLamL2,level);
%     caxis([min(min(min(dLamH3_2))),max(max(max(dLamH3_2)))]);
    hold on
    dLamLr2 = reshape(dLamH3r_2(:,i,:),size(dLamH3r_2,1),size(dLamH3r_2,3));
    if mean(mean(DLamLr2)) == 0
        level = 0;
    else
        level = 2;
    end
    polarcont(R,Theta_fr2,dLamLr2,level);
%     caxis([min(min(min(dLamH3r_2))),max(max(max(dLamH3r_2)))]);
    polar(repmat(Theta_f2.',2,1),repmat(holder,1,length(Theta_f2)),'--r');
    polar(repmat(Theta_fr2(1),2,1),holder,'-b');
    polar(repmat(Theta_fr2(3),2,1),holder,'-g');
    polar(repmat(Theta_fr2(4),2,1),holder,'-k');
    polar(repmat(Theta_fr2(6),2,1),holder,'-m');
    caxis([-20*10^-3 20*10^-3]);
    title(['c)'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    
%     ax9 = subplot(2,8,[5,6,13,14]);
%     fdLamRzr2 = cat(2, squeeze(fdLamRz2(:,:,i)),squeeze(fdLamRrr2(:,:,i)));
%     polarcont(R,Theta_final_2,fdLamRzr2,level);
% %     caxis([fmin_2,fmax_2]);
%     hold on
%     polar(repmat(Theta_f2.',2,1),repmat(holder,1,length(Theta_f2)),'--r');
%     polar(repmat(Theta_fr2(13),2,1),holder,'-b')
%     polar(repmat(Theta_fr2(9),2,1),holder,'-g');
%     polar(repmat(Theta_fr2(5),2,1),holder,'-k');
%     polar(repmat(Theta_fr2(2),2,1),holder,'-m');
%     caxis([-20*10^-3 20*10^-3]);
%     title(['f)'],'interpreter','latex','FontSize',15);
    
    ax10 = subplot(2,8,[7,8,15,16]);
    ddLamRzr2 = cat(2, squeeze(ddLamRz2(:,:,i)),squeeze(ddLamRrr2(:,:,i)));
    if mean(mean(ddLamRzr2)) == 0
        level = 0;
    else
        level = 2;
    end    
    polarcont(R,Theta_final_2t,ddLamRzr2,level);
    hold on
    polar(repmat(Theta_f2.',2,1),repmat(holder,1,length(Theta_f2)),'--r');
    polar(repmat(Theta_fr2(1),2,1),holder,'-b');
    polar(repmat(Theta_fr2(3),2,1),holder,'-g');
    polar(repmat(Theta_fr2(4),2,1),holder,'-k');
    polar(repmat(Theta_fr2(6),2,1),holder,'-m');
    caxis([-20*10^-3 20*10^-3]);
    title(['d)'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');
    
    colormap(gray);
    saveas(fig, fullfile('/home/kevin/Desktop/polar_0.2_0.5_z//', ['polarcont_comb',num2str(i),'.jpg']));
end

clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax 11
clear DLamC5 DLamC_fmax5 DLamC_fmin5
clear DLamC4 DLamC_fmax4 DLamC_fmin4
clear DLamC19 DLamC_fmax19 DLamC_fmin19
clear DLamC10 DLamC_fmax10 DLamC_fmin10

%% 3D
disp('3D');
fig = figure('visible','off')
set(gcf,'pos',[10, 10, 800*2, 2*300]);


fdLamH3_2 = zeros(size(dLamH3_2,1),size(dLamH3_2,3),size(dLamH3_2,2));
fdLamH3r_2 = zeros(size(dLamH3r_2,1),size(dLamH3r_2,3),size(dLamH3r_2,2));

for i = 1:size(dLamH3_2,2)
    fdLamH3_2(:,:,i) = squeeze(dLamH3_2(:,i,:));
end
for i = 1:size(dLamH3r_2,2)
    fdLamH3r_2(:,:,i) = squeeze(dLamH3r_2(:,i,:));
end

ImfdLamH3_2 = zeros(size(ImdLamH3_2,1),size(ImdLamH3_2,3),size(ImdLamH3_2,2));
ImfdLamH3r_2 = zeros(size(ImdLamH3r_2,1),size(ImdLamH3r_2,3),size(ImdLamH3r_2,2));

for i = 1:size(ImdLamH3_2,2)
    ImfdLamH3_2(:,:,i) = squeeze(ImdLamH3_2(:,i,:));
    ImfdLam3_2 = ImfdLamH3_2(:,:,i);
    ImfdLam3_2(ImfdLam3_2>0) = 1;
    ImfdLam3_2(ImfdLam3_2<0) = -1;
    ImfdLamH3_2(:,:,i) = ImfdLam3_2;
end
for i = 1:size(ImdLamH3r_2,2)
    ImfdLamH3r_2(:,:,i) = squeeze(ImdLamH3r_2(:,i,:));
    ImfdLam3r_2 = ImfdLamH3r_2(:,:,i);
    ImfdLam3r_2(ImfdLam3r_2>0) = 1;
    ImfdLam3r_2(ImfdLam3r_2<0) = -1;
    ImfdLamH3r_2(:,:,i) = ImfdLam3r_2;
end

a = size(fdLamH3_2,1);
b = size(fdLamH3_2,2);

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = R(j)*cos(Theta_f2(k));
        y(j,k) = R(j)*sin(Theta_f2(k));
    end
end

x = repmat(x,1,1,size(fdLamH3_2,3));
y = repmat(y,1,1,size(fdLamH3_2,3));
z = zeros(size(fdLamH3_2));
for i =1:size(fdLamH3_2,3)
    z(:,:,i) = i;
end
    
fdLamH3_2 = squeeze(fdLamH3_2);
fdLamH3_2 = smooth3(fdLamH3_2);
pfdLamH3_2 =(fdLamH3_2>0);
pfdLamH3_2 = pfdLamH3_2.*fdLamH3_2; 
nfdLamH3_2 =(fdLamH3_2<0);
nfdLamH3_2 = nfdLamH3_2.*fdLamH3_2;

ax1 = subplot(2,8,[1,2,3]);
[fo, vo] = isosurface(x,y,z,pfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',fo,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');



%circle  
a = size(fdLamH3r_2,1);
b = size(fdLamH3r_2,2);

xr = zeros(a,b);
yr = zeros(a,b);

for j = 1:a
    for k = 1:b
        xr(j,k) = R(j)*cos(Theta_fr2(k));
        yr(j,k) = R(j)*sin(Theta_fr2(k));
    end
end

xr = repmat(xr,1,1,size(fdLamH3r_2,3));
yr = repmat(yr,1,1,size(fdLamH3r_2,3));
z = zeros(size(fdLamH3r_2));
for i =1:size(fdLamH3r_2,3)
    z(:,:,i) = i;
end
    

fdLamH3r_2 = squeeze(fdLamH3r_2);
fdLamH3r_2 = smooth3(fdLamH3r_2);

%recall Remainder matrix is that of flipped, thus <0 instead of >0
pfdLamH3r_2 =(fdLamH3r_2<0);
pfdLamH3r_2 = pfdLamH3r_2.*fdLamH3r_2; 
nfdLamH3r_2 =(fdLamH3r_2>0);
nfdLamH3r_2 = nfdLamH3r_2.*fdLamH3r_2;

[f1, v1] = isosurface(xr,yr,z,pfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nfdLamH3r_2);
% [f6, v6] = reducepatch(f1,v1,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');

C = [0,0, -1];
Theta = 0:0.01:2*pi;
X = C(1)+ 19.65*cos(Theta);
Y = C(2)+ 19.65*sin(Theta);
Z = C(3)+zeros(size(X));
plot3(X,Y,Z);

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
    title(['a) Linear int: Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
ax10 = subplot(2,8,4);
[fo, vo] = isosurface(x,y,z,pfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',fo,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');
[f1, v1] = isosurface(xr,yr,z,pfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nfdLamH3r_2);
% [f6, v6] = reducepatch(f1,v1,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
view(45,90);


clear f1 v1 pr2 fo vo p2 pr6 f6 v6 p5 f5 v5

%
%
ax2= subplot(2,8,[9,10,11]);
ImfdLamH3_2 = squeeze(ImfdLamH3_2);
ImfdLamH3_2 = smooth3(ImfdLamH3_2);
pImfdLamH3_2 =(ImfdLamH3_2>0);
pImfdLamH3_2 = pImfdLamH3_2.*ImfdLamH3_2; 
nImfdLamH3_2 =(ImfdLamH3_2<0);
nImfdLamH3_2 = nImfdLamH3_2.*ImfdLamH3_2; 

[fo, vo] = isosurface(x,y,z,pImfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nImfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

ImfdLamH3r_2 = squeeze(ImfdLamH3r_2);
ImfdLamH3r_2 = smooth3(ImfdLamH3r_2);
pImfdLamH3r_2 =(ImfdLamH3r_2<0);
pImfdLamH3r_2 = pImfdLamH3r_2.*ImfdLamH3r_2; 
nImfdLamH3r_2 =(ImfdLamH3r_2>0);
nImfdLamH3r_2 = nImfdLamH3r_2.*ImfdLamH3r_2; 


[f1, v1] = isosurface(xr,yr,z,pImfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImfdLamH3r_2);
% [f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
    title(['b) Linear int: Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
%
ax11 = subplot(2,8,12);
[fo, vo] = isosurface(x,y,z,pImfdLamH3_2);
% [fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all
[f5, v5] = isosurface(x,y,z,nImfdLamH3_2);
% [f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

[f1, v1] = isosurface(xr,yr,z,pImfdLamH3r_2);
% [f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImfdLamH3r_2);
% [f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
view(45,90);

clear f1 v1 pr2 fo vo p2 p3 pr3 f2 v2 f3 v3 pr6 f6 v6 p5 f5 v5

% %
% ax3 = subplot(1,6,3);
% [fo, vo] = isosurface(x,y,z,pfdLamH3_2);
% p2 = patch('Faces',fo,'Vertices',vo);
% set(p2, 'FaceColor','red','EdgeColor','none');
% hold all
% [f5, v5] = isosurface(x,y,z,nfdLamH3_2);
% p5 = patch('Faces',f5,'Vertices',v5);
% set(p5, 'FaceColor','blue','EdgeColor','none');
% 
% [f1, v1] = isosurface(xr,yr,z,pfdLamH3r_2);
% pr2 = patch('Faces',f1,'Vertices',v1);
% set(pr2, 'FaceColor','red','EdgeColor','none');
% 
% [f6, v6] = isosurface(xr,yr,z,nfdLamH3r_2);
% pr6 = patch('Faces',f6,'Vertices',v6);
% set(pr6, 'FaceColor','blue','EdgeColor','none');
% 
% [f2, v2] = isosurface(x,y,z,pImfdLamH3_2);
% p3 = patch('Faces',f2,'Vertices',v2);
% set(p3, 'FaceColor','red','EdgeColor','none');
% 
% [f7, v7] = isosurface(x,y,z,nImfdLamH3_2);
% p7 = patch('Faces',f7,'Vertices',v7);
% set(p7, 'FaceColor','blue','EdgeColor','none');
% 
% [f3, v3] = isosurface(xr,yr,z,pImfdLamH3r_2);
% pr3 = patch('Faces',f3,'Vertices',v3);
% set(pr3, 'FaceColor','red','EdgeColor','none');
% 
% [f8, v8] = isosurface(xr,yr,z,nImfdLamH3r_2);
% pr8 = patch('Faces',f8,'Vertices',v8);
% set(pr8, 'FaceColor','blue','EdgeColor','none');
% plot3(X,Y,Z);
% 
% view(3);
% axis tight
% daspect([1,1,1])
% camlight 
% lighting gouraud
% title(['c) No int'],'interpreter','latex','FontSize',15);
% 
% clear f1 v1 pr2 fo vo p2 p3 pr3 f2 v2 f3 v3 f5 f6 f7 f8 p5 p6 p7 p8
% clear pr8 p7 pr6 p5 
%

%
%
a = size(ddLamRz2,1);
b = size(ddLamRz2,2);

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = R(j)*cos(Theta_final2t(k));
        y(j,k) = R(j)*sin(Theta_final2t(k));
    end
end

x = repmat(x,1,1,size(ddLamRz2,3));
y = repmat(y,1,1,size(ddLamRz2,3));
z = zeros(size(ddLamRz2));
for i =1:size(ddLamRz2,3)
    z(:,:,i) = i;
end
    
ddLamRz2 = squeeze(ddLamRz2);
ddLamRz2 = smooth3(ddLamRz2);
pddLamRz2 =(ddLamRz2>0);
pddLamRz2 = pddLamRz2.*ddLamRz2; 
nddLamRz2 =(ddLamRz2<0);
nddLamRz2 = nddLamRz2.*ddLamRz2; 


ax4 = subplot(2,8,[5,6,7]);
[fo, vo] = isosurface(x,y,z,pddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

a = size(ddLamRrr2,1);
b = size(ddLamRrr2,2);

xr = zeros(a,b);
yr = zeros(a,b);

for j = 1:a
    for k = 1:b
        xr(j,k) = R(j)*cos(Theta_finalr2t(k));
        yr(j,k) = R(j)*sin(Theta_finalr2t(k));
    end
end

xr = repmat(xr,1,1,size(ddLamRrr2,3));
yr = repmat(yr,1,1,size(ddLamRrr2,3));
z = zeros(size(ddLamRrr2));
for i =1:size(ddLamRrr2,3)
    z(:,:,i) = i;
end
    

ddLamRrr2 = squeeze(ddLamRrr2);
ddLamRrr2 = smooth3(ddLamRrr2);
pddLamRrr2 =(ddLamRrr2<0);
pddLamRrr2 = pddLamRrr2.*ddLamRrr2; 
nddLamRrr2 =(ddLamRrr2>0);
nddLamRrr2 = nddLamRrr2.*ddLamRrr2; 


[f1, v1] = isosurface(xr,yr,z,pddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');

C = [0,0, -1];
Theta = 0:0.01:2*pi;
X = C(1)+ 19.65*cos(Theta);
Y = C(2)+ 19.65*sin(Theta);
Z = C(3)+zeros(size(X));
plot3(X,Y,Z);

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
title(['c)',num2str(plevel),'-poly int: Real'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
ax12 = subplot(2,8,8);
[fo, vo] = isosurface(x,y,z,pddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

[f1, v1] = isosurface(xr,yr,z,pddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);

axis tight
daspect([1,1,1])
camlight 
camlight right
% camlight left
lighting gouraud

view(45,90);
clear f1 v1 pr2 fo vo p2 pr6 p5 f5 v5 f6 v6 

%
%
ax5= subplot(2,8,[13,14,15]);
ImddLamRz2 = squeeze(ImddLamRz2);
ImddLamRz2 = smooth3(ImddLamRz2);
pImddLamRz2 =(ImddLamRz2>0);
pImddLamRz2 = pImddLamRz2.*ImddLamRz2; 
nImddLamRz2 =(ImddLamRz2<0);
nImddLamRz2 = nImddLamRz2.*ImddLamRz2; 

[fo, vo] = isosurface(x,y,z,pImddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nImddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

ImddLamRrr2 = squeeze(ImddLamRrr2);
ImddLamRrr2 = smooth3(ImddLamRrr2);
pImddLamRrr2 =(ImddLamRrr2<0);
pImddLamRrr2 = pImddLamRrr2.*ImddLamRrr2; 
nImddLamRrr2 =(ImddLamRrr2>0);
nImddLamRrr2 = nImddLamRrr2.*ImddLamRrr2; 

[f1, v1] = isosurface(xr,yr,z,pImddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');

plot3(X,Y,Z);
view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud
title(['d)',num2str(plevel),'-poly int: Imag'],'fontweight','normal','interpreter','latex','FontSize',22,'fontname','arial');


%
%
ax14 = subplot(2,8,16);
[fo, vo] = isosurface(x,y,z,pImddLamRz2);
[fo, vo] = reducepatch(fo,vo,.5);
p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');
hold all 
[f5, v5] = isosurface(x,y,z,nImddLamRz2);
[f5, v5] = reducepatch(f5,v5,.5);
p5 = patch('Faces',f5,'Vertices',v5);
set(p5, 'FaceColor','blue','EdgeColor','none');

[f1, v1] = isosurface(xr,yr,z,pImddLamRrr2);
[f1, v1] = reducepatch(f1,v1,.5);
pr2 = patch('Faces',f1,'Vertices',v1);
set(pr2, 'FaceColor','red','EdgeColor','none');

[f6, v6] = isosurface(xr,yr,z,nImddLamRrr2);
[f6, v6] = reducepatch(f6,v6,.5);
pr6 = patch('Faces',f6,'Vertices',v6);
set(pr6, 'FaceColor','blue','EdgeColor','none');
plot3(X,Y,Z);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud

view(45,90);

%
%
% ax6 = subplot(1,6,6);
% [fo, vo] = isosurface(x,y,z,ddLamRz2);
% [fo, vo] = reducepatch(fo,vo,.5);
% p2 = patch('Faces',fo,'Vertices',vo);
% set(p2, 'FaceColor','blue','EdgeColor','none');
% hold all
% [f1, v1] = isosurface(xr,yr,z,ddLamRrr2);
% [f1, v1] = reducepatch(f1,v1,.5);
% pr2 = patch('Faces',f1,'Vertices',v1);
% set(pr2, 'FaceColor','red','EdgeColor','none');
% [f2, v2] = isosurface(x,y,z,ImddLamRz2);
% [f2, v2] = reducepatch(f2,v2,.5);
% p3 = patch('Faces',f2,'Vertices',v2);
% set(p3, 'FaceColor','green','EdgeColor','none');
% [f3, v3] = isosurface(xr,yr,z,ImddLamRrr2);
% [f3, v3] = reducepatch(f3,v3,.5);
% pr3 = patch('Faces',f3,'Vertices',v3);
% set(pr3, 'FaceColor','green','EdgeColor','none');
% plot3(X,Y,Z);
% view(3);
% axis tight
% daspect([1,1,1])
% camlight 
% lighting gouraud
% title(['f)',num2str(plevel),'-poly int'],'interpreter','latex','FontSize',15);

campos(ax1,[-35, -45, -0]);
campos(ax2,[-35, -45, -0]);
campos(ax4,[-35, -45, -0]);
campos(ax5,[-35, -45, -0]);


camroll(ax1,90);
camroll(ax2,90);
camroll(ax4,90);
camroll(ax5,90);
% %%
% 
% saveas(fig, fullfile('/home/kevin/Desktop//', ['isosurface_f.fig']));
saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_3.jpg']));
% 
% view(ax1,[-90,0]);
% view(ax2,[-90,0]);
% % view(ax3,[-90,0]);
% view(ax4,[-90,0]);
% view(ax5,[-90,0]);
% % view(ax6,[-90,0]);
% 
% 
% saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_5.jpg']));
% view(ax1,[-60,25]);
% view(ax2,[-60,25]);
% % view(ax3,[-60,25]);
% view(ax4,[-60,25]);
% view(ax5,[-60,25]);
% % view(ax6,[-60,25]);
% 
% 
% saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_6.jpg']));
% 
% view(ax1,[-60,90]);
% view(ax2,[-60,90]);
% % view(ax3,[-60,90]);
% view(ax4,[-60,90]);
% view(ax5,[-60,90]);
% % view(ax6,[-60,90]);
% 
% saveas(fig,fullfile('/home/kevin/Desktop//',['isosurface_7.jpg']));
% 
% %
%
%
%}
%
%%     contour(xCon/200,yCon/200,DLamC(:,:,k),3)%% OMD -mode calculation function
function [xCon, yCon, St, DLamSt] = slicedata(omd_rank,St_chosen,iptx_start,iptx_end,snapshotsV_temp,visz,visy)
snapshotsV_f = snapshotsV_temp(:,iptx_start:iptx_end);
setlength = size(snapshotsV_f,2);
if setlength > 500
    setlength = 500;
end

%limiting the number of setlenght (snapshot) lengths in order to 
%decrease the computational time, REMOVE this later for higher quality 
%flow data (might increases the accuracy)

visv = snapshotsV_f(:,1:setlength);
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
    interpolated_visv(i,:)= X(i,:);
end


%prepare data into form for contour plot
p = length(xdummy);
n = length(ydummy);
xCon = reshape(xnew,[n,p]);
yCon = reshape(ynew,[n,p]);

B = interpolated_visv(:,1:end-1);
A = interpolated_visv(:,2:end);

[L,M] = omd(A,B,omd_rank,[],'alternating'); 
[Y,D_nu]=eig(M);%eigenvector Y and correspondong eigenvalues D_nu of F_dmd

Vand=zeros(size(M));
for loop1=1:size(B,2)
    Vand(:,loop1)=diag(D_nu).^(loop1-1); %creating vandermonde with increasing powers
end
U=L;
Psi = U*Y;

Z(:,1) =(Y^-1)*U'*B(:,1);

for i = 1:size(Vand,1)
    for j = 1:size(Vand,1)
        alpha_omd(j,i) = Z(j,1)*Vand(j,i);
    end
end

% 
alpha_omd_abs_I = abs(real(alpha_omd));
alpha_omd_abs_I_t = alpha_omd_abs_I.';

I = sum(alpha_omd_abs_I_t);
[C,T] = sort(I,'descend');

dt = 1/720;
freq = 0.1965/15/(2*pi)*imag(log(diag(D_nu))/dt);


%%%% closes to 0.2
[c St] = sort(abs(St_chosen-freq));
k = St(1);
St = freq(St(1));
% disp(St);
DLamSt = reshape(Psi(:,k),[n,p]);

% 
%%%% I method
% dummy2 = [c, St02,alpha_omd];
% dummy2 = dummy2(dummy2(:,1)<0.19,:);
% 
% sum_alpha = sum(dummy2(:,3:end),2);
% [cdummy St02dummy] = sort(sum_alpha);
% disp(St02dummy(1));
% St02 = dummy2(St02dummy(1),2);

% [d St006] = sort(abs(0.06-freq));

% dummy6 = [d, St006,alpha_omd];
% dummy6 = dummy6(dummy6(:,1)<0.059,:);
% 
% sum_alpha = sum(dummy6(:,3:end),2);
% [cdummy St006dummy] = sort(sum_alpha);
% St006 = dummy6(St006dummy(1),2);


%calculating DLam (DMD modes for St02 and St006)
% DLamSt02 = reshape(Psi(:,k),[n,p]);
% 
% i = St006;
% DLamSt006 = reshape(Psi(:,i),[n,p]);
%     


end

%% calc_cop_t function
function [ R, Theta] = calc_cop_t( p )

% calc_cop_t Function to give the coordinates of the COP for pressure in mm
% and radians
% timeseries

[np,nt]=size(p);

% set up the vector of radius values
rr=linspace(11,88,8);
r=repmat(rr',[8 1]);
th=[6*pi/4 7*pi/4 0 pi/4 pi/2 3*pi/4 pi 5*pi/4];
thet=reshape(repmat(th,[8 1]),64,1);

% Set up the vector of areas corresponding to each tapping
area = zeros(64,1);
for it=1:length(r)
    if r(it) == 11
        ro = mean([r(it+1),r(it)]);
        area(it) = pi*ro^2/8;
    elseif r(it) == 88
        ro = 98.25; ri = mean([r(it-1),r(it)]);
        area(it) = pi*(ro^2 - ri^2)/8;
    else
        ro = mean([r(it+1),r(it)]); ri = mean([r(it-1),r(it)]);
        area(it) = pi*(ro^2 - ri^2)/8;
    end
end
sum(area);
area=repmat(area,[1 nt]);
size(area);


% Evaluate the vector of forces, the total force and the terms FRcos(theta)
% and FRsin(theta)
f = p.*area;
F = sum(f);
size(F);
FRcos = (f'*(r.*cos(thet)))';
FRsin = (f'*(r.*sin(thet)))';

size(FRcos);
% Evaluate Theta and R
% Evaluate Theta in the range 0-pi for a continuous time series

% xcop=FRcos./F;
% ycop=FRsin./F;
% R = sqrt((FRcos.^2+FRsin.^2)./F.^2);
% Theta = acos(FRcos./(F.*R));
% Non-dimensionalise R
Theta = atan2(FRsin,FRcos);
R = abs(FRcos./cos(Theta)./F);
end


%% function for polar contour
function [C,h] = polarcont(r,theta,z,N,s)

[a,b] = size(z);

if a ~= length(r)
    error('r is not the same length as the first dimension of z')
end

if b ~= length(theta)
    error('theta is not the same length as the second dimension of z')
end

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = r(j)*cos(theta(k));
        y(j,k) = r(j)*sin(theta(k));
    end
end

if nargin == 3
    [C,h] = contourf(x,y,z);
elseif nargin == 4
    [C,h] = contourf(x,y,z,N);
elseif nargin == 5
    [C,h] = contourf(x,y,z,N,s);
else
    error('Incorrect number of inputs')
end 
end

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
%
%
%}

