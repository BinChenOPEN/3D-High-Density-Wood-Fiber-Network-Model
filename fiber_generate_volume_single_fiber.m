function [volume,volume_solid,middle_line] = fiber_generate_volume_single_fiber(radius,cell_wall_thick,Length,angle)
%% this script can generate a the 3D volume that includes a single fiber
fiber_axis = [sin(angle),cos(angle),0]; % axis of the fiber, where angle is the angle of the fiber

size_volume = ceil([Length,Length,radius(2)*2*(1)+1]); % initial volume of one fiber
Initial_point = [Length/2+2*Length*(rand(1)-0.5),Length/2+2*Length*(rand(1)-0.5),round(size_volume(3)/2)];
[x,y,z] = ndgrid(1:size_volume(1),1:size_volume(2),1:size_volume(3));

%% two type of volume will be generated
volume = zeros(size_volume);
volume_solid = volume;

vector_point = [x(:),y(:),z(:)]-Initial_point;
projection_axis = vector_point(:,[1,2])*fiber_axis(1,[1,2])';
projection_plane = vector_point - projection_axis*fiber_axis;

rz2 = projection_plane(:,3).^2; %% The projection of a point in z axis
rxy2 = sum(projection_plane.^2,2)-rz2; %% The projection of a point in xy plane

%% to check is that point in or outside of the fiber
dist_indicator_out = rxy2/radius(1)^2+rz2/radius(2)^2;
dist_indicator_in = rxy2/(radius(1)-cell_wall_thick)^2+rz2/(radius(2)-cell_wall_thick)^2;

volume(find(dist_indicator_out<=1 & dist_indicator_in>=1 & abs(projection_axis)<=Length*1.5)) = 1;
volume_solid(find(dist_indicator_out<=1 & abs(projection_axis)<=1.5*Length)) = 1;

fiber_axis = [sin(angle),cos(angle)];
[x,y] = ndgrid(1:size_volume(1),1:size_volume(2));
vector_point = [x(:),y(:)]-Initial_point([1,2]);
projection_axis = vector_point*fiber_axis';
projection_plane = vector_point - projection_axis*fiber_axis;
rz2 = sqrt(sum(projection_plane.^2,2));

% figure, surf(reshape(r2,size_volume([1,2]))), shading interp, axis equal,

indx_fiber = find(rz2<2 & abs(projection_axis)<2*Length*1.1);

middle_line.fiber_axis      = fiber_axis;
middle_line.Initial_point   = Initial_point-[0,0,1];
middle_line.fiber_length    = Length*1.5;
middle_line.indx_fiber      = indx_fiber;
% middle_line.projection_axis = projection_axis;

volume = volume(:,:,2:end-1);
volume_solid = volume_solid(:,:,2:end-1);


