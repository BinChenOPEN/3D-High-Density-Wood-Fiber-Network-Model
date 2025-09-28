% This project is used to generate realistic 3D wood fiber network
% microstructure using a geometrical method

% Author: Bin Chen
% E-mail: binchen@kth.se
% Update: 2025/09/28

clear all, close all,

%% set the basic parameters
fiber_num = 50;
cell_wall_thick = 4;
length_small = 500; % size of the volume in x-y plane
length_large = length_small+200; % generate a large one and then keep only the middle
fiber_align_mode = 1;% choose between 1 and 2 for isotropic and anisotropic fiber distribution
fiber_width_mean = 60;
fiber_width_variation = 20;
mold_type = 1; % The value is selected between 1 and 2. 1 is for spherical shape,
% while 2 is the wedge shape
is_molded_WFN = 1; % The value is selected between 0 and 1. 1 will genreate
% molded WFN, while 0 will skip that.
angle_std = pi/8; % angle standard deviation for gaussian distribution.


%% The folders to save the results
folder_save_all = ['Resuts/random_',num2str(fiber_num)];

folder_save_results = fullfile(folder_save_all,'results');
mkdir(folder_save_results);
folder_single_solid_fibers = fullfile(folder_save_all,'fiber_volume_solid');

% mkdir(folder_single_fibers);
mkdir(folder_single_solid_fibers)
ft = 18;

%% generate single fibers
middle_line_file = fullfile(folder_save_all,'middle_line_fiber.mat');

half_offset = (length_large-length_small)/2;

for i = 1:fiber_num
    ratio_volume = 0;
    
    radius(1) = [round(fiber_width_mean+2*fiber_width_variation*(rand(1)-0.5))]; %% parameters to control the fiber size
    radius(2) = round(radius(1)/fiber_width_mean)+cell_wall_thick;
    
    filename_solid = fullfile(folder_single_solid_fibers,['fiber_',num2str(i,'%04d'),'.mat']);
    
    
    file_fiber_solid_all{i} = filename_solid;
    if ~exist(file_fiber_solid_all{i})
        
        switch fiber_align_mode
            case 1
                angle = pi*(rand(1)-0.5);
            case 2
                angle = angle_std*(randn(1));
        end
        
        %  if in the beginning the fibers occupy only a small region. We
        %  have to skii that fiber, until in the end the fiber is roughy
        %  inside the volume
        while ratio_volume<0.02
            [volume,volume_solid,middle_line] = fiber_generate_volume_single_fiber(radius,cell_wall_thick,length_large,angle);
            
            x = half_offset+1:half_offset+length_small;
            y = half_offset+1:half_offset+length_small;
            volume_center_solid = volume_solid(x,y,:);
            ratio_volume = sum(sum(sum(volume_center_solid)))/prod(size(volume_center_solid));
        end
        
        middle_line_fiber.angle(i) = angle;
        middle_line_fiber.radius(i,:) = radius;
        middle_line_fiber.middle_line{i} = middle_line;
        
        save(filename_solid,'volume_solid');
        if i == fiber_num
            save(middle_line_file,'middle_line_fiber');
        end
    end
end

%% calculate the summation of the thickness at each layer
sum_thick_all_file = fullfile(folder_save_all,'sum_thick_all.mat');
if ~exist(sum_thick_all_file)
    
    for i = 1:fiber_num
        
        % load(filemale_all{i});
        load(file_fiber_solid_all{i})
        sum_thick_all(:,:,i) = sum(volume_solid,3);
    end
    save(sum_thick_all_file,'sum_thick_all', '-v7.3');
else
    load(sum_thick_all_file)
end
param_file = fullfile(folder_save_all,'Params.mat');
Params.log_file = fullfile(folder_save_all,'Porosity.txt');
Params.param_file = param_file;

if ~exist(param_file)
    tic
    %% Optimize the fiber location to reduce the porosity
    Params.radius = radius;
    Params.cell_wall_thick = cell_wall_thick;
    if exist(param_file)
        
        load(param_file);
    end
    x_offset = zeros(2*size(sum_thick_all,3),1);
    Params.ratioEst_init = thick_optimize(x_offset,sum_thick_all,length_large,length_small,Params);
    if 1
        
        x_offset_l = -half_offset*ones(2*size(sum_thick_all,3),1);
        x_offset_u = half_offset*ones(2*size(sum_thick_all,3),1);
        options = optimoptions('particleswarm','SwarmSize',5,'Display','iter',...
            'MaxIterations',50,'UseParallel',true, 'FunctionTolerance',0.001);
        [x_offset_opt,fval,exitflag] = particleswarm(@(x_offset)thick_optimize(x_offset,sum_thick_all,length_large,length_small,Params),...
            length(x_offset),x_offset_l,x_offset_u);
        Params.x_offset_opt = x_offset_opt;
    end
    toc
else
    load(param_file);
end

results_file = fullfile(folder_save_results,'Results.mat');
indx_fiber  = 1:fiber_num;
x_offset_opt = Params.x_offset_opt;
load(middle_line_file);

%% obtain the initial sparse structure
Thickness = sum(middle_line_fiber.radius(indx_fiber,2))*2-length(middle_line_fiber.radius);
volume_all_solid = zeros(length_small,length_small,Thickness,'single');
x_offset_round = reshape(round(x_offset_opt),[],2);
sum_thick_all_new = zeros(length_small,length_small,size(sum_thick_all,3));
half_offset = (length_large-length_small)/2;


k = 0;
thick_accum = 0;
for i = indx_fiber
    k = k+1;
    %         load(filemale_all{i});
    load(file_fiber_solid_all{i})
    
    x = half_offset+1+x_offset_round(i,1):x_offset_round(i,1)+half_offset+length_small;
    y = half_offset+1+x_offset_round(i,2):x_offset_round(i,2)+half_offset+length_small;
    
    sum_thick_all_new(:,:,k) = sum_thick_all(x,y,i);
    
    volume_solid_move = single(volume_solid(x,y,:));
    %             volume_all_solid = cat(3,volume_all_solid,volume_solid_move);
    volume_all_solid(:,:,thick_accum+(1:size(volume_solid_move,3))) ...
        = volume_solid_move;
    thick_accum = thick_accum + size(volume_solid_move,3);
end
sum_all = sum(sum_thick_all_new,3);
[max_thick,indx] = max(sum_all(:));

Params.sum_thick_all_new = sum_thick_all_new;
Params.sum_all = sum_all;
% save (fullfile(folder_save_results{t},'volume_all_solid.mat'),'volume_all_solid', '-v7.3');

% saveas(volume_all,'volume_all.tif');
L_compress_thick = max_thick;
Lc = L_compress_thick;
size_vol_all = size(volume_all_solid);
L0 = size_vol_all(3);
Lf = sum(volume_all_solid,3);
Params.Lf = Lf;
Params.porosity_expect = 1-mean(Lf(:))/Lc;
%     figure,surf(Lf),shading interp, axis equal, view([0,90]),axis off,colorbar,set(gcf,'color','w'),set(gca,'fontsize',ft)
%     exportgraphics(gcf,fullfile(folder_save_results{t},'Lf_optimize.pdf'))


%% compute the disp fields
w_all = zeros(size_vol_all);
sum_sub_volume = zeros(size_vol_all(1:2));

comp_strain = (L0-Lc)./(L0-Lf);
for i = 1:size(volume_all_solid,3)
    sum_sub_volume = sum_sub_volume+double(volume_all_solid(:,:,i));
    w_all(:,:,i) = -(i-sum_sub_volume).*comp_strain;
end
% smooth the displacement fields
h = fspecial3('average',[11,11,1]);
w_all_filter = convn(w_all,h,'same');
%% Generate compressed volume structure
clear volum_compress volum_compress_solid;

% compress layer by layer using interpolation
[y,z] = ndgrid(1:size_vol_all(2),1:size_vol_all(3));
parfor i = 1:size_vol_all(1)
    
    z_new = z+squeeze(w_all_filter(i,:,:));
    slice_volume = double(squeeze(volume_all_solid(i,:,:)));
    
    F = scatteredInterpolant(y(:),z_new(:),slice_volume(:),'linear');
    
    [y_comp,z_comp] = ndgrid(1:size_vol_all(2),1:Lc);
    Vq = F(y_comp(:),z_comp(:));
    volum_compress_solid(i,:,:) = reshape(Vq,[size_vol_all(2),Lc]);
end

volum_compress_solid_center = volum_compress_solid(10:end-9,10:end-9,:);
volum_compress_solid_center_new = volum_compress_solid_center;
volum_compress_solid_center_new(find(volum_compress_solid_center>0.5)) = 1;
volum_compress_solid_center_new(find(volum_compress_solid_center<=0.5)) = 0;
porosity_solid = 1-(sum(volum_compress_solid_center_new(:))/prod(size(volum_compress_solid_center_new)));

Params.porosity_solid = porosity_solid; % The porosity of the materil
save(fullfile(folder_save_results,'Params.mat'),'Params');
save(fullfile(folder_save_results,'volum_compress_solid_center_new.mat'),'volum_compress_solid_center_new', '-v7.3');

fprintf('%n structures have been generated...\n')


%% Generate the molded WFN with a given molded shape
if is_molded_WFN
    switch mold_type
        % generate the 3D displacement map
        case 1
            %% curved mold
            [x,y,z] = ndgrid(1:size(w_all_filter,1),1:size(w_all_filter,2),1:size(w_all_filter,3));
            x = x-(1+size(w_all_filter,1))/2;
            y = y-(1+size(w_all_filter,1))/2;
            
            dist = sqrt(x.^2+y.^2);
            r = 2000;
            Mold_Disp = r-r.*cos(asin(dist/r));
        case 2
            %% weidge mold
            [x,y,z] = ndgrid(1:size(w_all_filter,1),1:size(w_all_filter,2),1:size(w_all_filter,3));
            x = x-(1+size(w_all_filter,1))/2;
            y = y-(1+size(w_all_filter,1))/2;
            Mold_Disp = abs(y)*0.2;
    end
    
    %% Combine the compressed displacement maps
    w_all_curve = w_all_filter+Mold_Disp;
    save(fullfile(folder_save_results,'w_all_curve.mat'),'w_all_curve', '-v7.3');
    
    Lc_curve = ceil(Lc+max(Mold_Disp(:)));
    
    %% image interpolation to get the compressed WFNs
    parfor i = 1:size_vol_all(1)
        [y,z] = ndgrid(1:size_vol_all(2),1:size_vol_all(3));
        z_new = z+squeeze(w_all_curve(i,:,:));
        slice_volume = double(squeeze(volume_all_solid(i,:,:)));
        
        F = scatteredInterpolant(y(:),z_new(:),slice_volume(:),'linear','none');
        
        [y_comp,z_comp] = ndgrid(1:size_vol_all(2),1:Lc_curve);
        Vq = F(y_comp(:),z_comp(:));
        volum_compress_curve(i,:,:) = reshape(Vq,[size_vol_all(2),Lc_curve]);
    end
    volum_compress_curve_clean = volum_compress_curve;
    [x,y,z] = ndgrid(1:size(volum_compress_curve,1),...
        1:size(volum_compress_curve,2),...
        1:size(volum_compress_curve,3));
    % The compressed WFN network is not perfect due to interpolation. We need
    % to remove all the bad points outside the compressed WFNs
    volum_compress_curve_clean(z<Mold_Disp(:,:,1:size(volum_compress_curve,3))) = 0;
    volum_compress_curve_clean(z>Lc+Mold_Disp(:,:,1:size(volum_compress_curve,3))) = 0;
    % Save the 3D model
    volum_compress_curve_clean_center = volum_compress_curve_clean(10:end-9,10:end-9,:);
    save (fullfile(folder_save_results,'volum_compress_curve_clean_center.mat'),'volum_compress_curve_clean_center', '-v7.3');
end

