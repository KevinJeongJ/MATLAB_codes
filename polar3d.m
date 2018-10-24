%% 3D
function [p2] = polar3d(dLamH, R, Theta)
disp('3D - visualization (Polar)');
fig = figure('visible','on') % turn the visibility of the plot on or off
set(gcf,'pos',[10, 10, 800, 300]);
%{
This code allows for the 3D polar isocontour to be created. 
    dLamH = 3D matrix of modes 
    mode = 2D matrix at each theta locations. its dimension is "diameter" x "height"
    for this code to work, sort the modes in increasing "theta"
    dLamH = values at "diameter" x " height" x "theta" 
    
    R = radius of the base 
    
    Theta = ascending theta values 
%}

fdLamH3 = zeros(size(dLamH,1),size(dLamH,3),size(dLamH,2));
%% creating an empty 3D matrix 

for i = 1:size(dLamH,2)
    fdLamH3(:,:,i) = squeeze(dLamH(:,i,:));
end

% creating the polar grid
a = size(fdLamH3,1);
b = size(fdLamH3,2);

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = R(j)*cos(Theta(k));
        y(j,k) = R(j)*sin(Theta(k));
    end
end

x = repmat(x,1,1,size(fdLamH3,3));
y = repmat(y,1,1,size(fdLamH3,3));
z = zeros(size(fdLamH3));

for i =1:size(fdLamH3,3)
    z(:,:,i) = i;
end
    
fdLamH3 = squeeze(fdLamH3);
fdLamH3 = smooth3(fdLamH3);


[fo, vo] = isosurface(x,y,z,fdLamH3);
% [fo, vo] = reducepatch(fo,vo,.5);
% only apply the above line if the data is too big for MATLAB to handle
% within given time frame 


p2 = patch('Faces',fo,'Vertices',vo);
set(p2, 'FaceColor','red','EdgeColor','none');

%% visual aspects

view(3);
axis tight
daspect([1,1,1])
camlight
camlight right
% camlight left
lighting gouraud

camroll(90);

end
