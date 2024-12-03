% function [cfar_output] = CACFAR(scaled_fft2d)
% % scaled_fft2d = log10(scaled_fft2d + 1);
%  % scaled_fft2d = scaled_fft2d(:, 1:32); 
% num_guard_cells = 2;    
% num_train_cells = 8;    
% pfa = 10^(-5);          
% g = num_train_cells*(pfa^(-1/num_train_cells) - 1); 
% 
% cfar_output = zeros(size(scaled_fft2d));
% cfar_threshold = zeros(size(scaled_fft2d));
% for row = 1:size(scaled_fft2d, 1)
%     for i = 1:size(scaled_fft2d, 2)
%         if i <= num_train_cells + num_guard_cells
%             train_region = scaled_fft2d(row, i + num_guard_cells + 1 : min(i + num_guard_cells + num_train_cells, size(scaled_fft2d, 2)));
%             noise_level = mean(train_region);
%             threshold = noise_level * g;
%         elseif i > size(scaled_fft2d, 2) - (num_train_cells + num_guard_cells)
%             train_region = scaled_fft2d(row, max(1, i - num_guard_cells - num_train_cells) : i - num_guard_cells - 1);
%             noise_level = mean(train_region);
%             threshold = noise_level * g;
%         else
%             train_region_left = scaled_fft2d(row, i - num_train_cells - num_guard_cells : i - num_guard_cells - 1);
%             train_region_right = scaled_fft2d(row, i + num_guard_cells + 1 : i + num_guard_cells + num_train_cells);
%             train_region = [train_region_left train_region_right];
%             noise_level = mean(train_region);
%             threshold = noise_level * g*2;
%         end
% 
%         cfar_threshold(row, i) = threshold;  
%         if scaled_fft2d(row, i) > threshold
%             cfar_output(row, i) = scaled_fft2d(row, i); 
%         else
%             cfar_output(row, i) = 0;  
%         end
%     end
% end
% 
% % fft2d_data_32=fft2d_data(:,1:32);
% % cfar_output_cmplx=zeros(size(fft2d_data_32,1),size(fft2d_data_32,2));
% % [row1,col1]=find(cfar_output);
% % cfar_output_cmplx(row1,col1)=fft2d_data_32(row1,col1);
% end

function [cfar_output,cfar_output_cmplx] = CACFAR(scaled_fft2d,fft2d_data)
scaled_fft2d = log10(scaled_fft2d + 1);
% scaled_fft2d = scaled_fft2d(:, 1:32); 
num_guard_cells = 2;    
num_train_cells = 20;    
pfa = 10^(-3);          
g = num_train_cells*(pfa^(-1/num_train_cells) - 1); 

cfar_output = zeros(size(scaled_fft2d));
cfar_threshold = zeros(size(scaled_fft2d));
for row = 1:size(scaled_fft2d, 1)
    for i = 1:size(scaled_fft2d, 2)
        if i <= num_train_cells + num_guard_cells
            train_region = scaled_fft2d(row, i + num_guard_cells + 1 : min(i + num_guard_cells + num_train_cells, size(scaled_fft2d, 2)));
            noise_level = mean(train_region);
            threshold = noise_level * g;
        elseif i > size(scaled_fft2d, 2) - (num_train_cells + num_guard_cells)
            train_region = scaled_fft2d(row, max(1, i - num_guard_cells - num_train_cells) : i - num_guard_cells - 1);
            noise_level = mean(train_region);
            threshold = noise_level * g;
        else
            train_region_left = scaled_fft2d(row, i - num_train_cells - num_guard_cells : i - num_guard_cells - 1);
            train_region_right = scaled_fft2d(row, i + num_guard_cells + 1 : i + num_guard_cells + num_train_cells);
            train_region = [train_region_left train_region_right];
            noise_level = mean(train_region);
            threshold = noise_level * g*2;
        end

        cfar_threshold(row, i) = threshold;  
        if scaled_fft2d(row, i) > threshold
            cfar_output(row, i) = scaled_fft2d(row, i); 
        else
            cfar_output(row, i) = 0;  
        end
    end
end

% fft2d_data_32=fft2d_data(:,1:32);
fft2d_data_32=fft2d_data;
cfar_output_cmplx=zeros(size(fft2d_data_32,1),size(fft2d_data_32,2));
[row1,col1]=find(cfar_output);
cfar_output_cmplx(row1,col1) = fft2d_data_32(row1,col1);
end