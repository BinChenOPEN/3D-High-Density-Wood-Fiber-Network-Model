function porosity = thick_optimize(x_offset,sum_thick_all,Length0,Length,Params)
% This script can optimize the fiber location in xy plane to reduce the porosity of the
% materials
half_offset = (Length0-Length)/2;
x_offset_round = reshape(round(x_offset),[],2);

sum_all = zeros(Length,Length);
for i = 1:size(sum_thick_all,3)
    x = half_offset+1+x_offset_round(i,1):x_offset_round(i,1)+half_offset+Length;
    y = half_offset+1+x_offset_round(i,2):x_offset_round(i,2)+half_offset+Length;
    sum_all = sum_all + sum_thick_all(x,y,i);
end
max_thick = max(sum_all(:));   
porosity = 1-mean(sum_all(:))/max_thick;
if 1
    if ~exist(Params.log_file)
        porosity_old = 1;
    else
        fid = fopen(Params.log_file, 'r');  % 'w' = write mode
        porosity_old = fscanf(fid, '%f');
        fclose(fid);
    end
        
    % if the porosity decrease, save the parameters, otherwise skip
    if (porosity_old-porosity)/porosity_old>0.01
        fid = fopen(Params.log_file, 'w');  % 'w' = write mode
        fprintf(fid, '%5f',porosity);
        fprintf('Estimated porosity: %f\n',porosity);
        Params.x_offset_opt = round(x_offset);
        Params.ratioEst_optim = porosity;
        save(Params.param_file,'Params');
        fclose(fid);
    end
end