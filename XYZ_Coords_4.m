clc;
clearvars;
close all;
load(".....path to radar data....");
%% 
% c1 = [];c2 = [];c3 = [];

for receiver = 1:3
    ov = overall_data(receiver:3:end, :, :);
    % b = [];
    c = [];
    for i = 1:100
        first_data = reshape(ov(i, :, :), 128, 64);
        fft_result = fft(first_data, [], 2);
        modified_fft_result = zeros(128, 64);
        % 1D FFT
        for j = 1:128
            skibidi = fft_result(j, :); 
            skibidi(1) = 0;
            modified_fft_result(j, :) = skibidi; 
        end
       

        %  2D FFT
        for k=1:64
                    modified_fft_result_2d(:,k) = fft(modified_fft_result(:,k),128);
        end

        scaled_fft2d=abs(modified_fft_result_2d);
        [cfar_output,cfar_output_cmplx] = CACFAR(scaled_fft2d,modified_fft_result_2d);
        modified_fft_result_2d_final = cfar_output_cmplx;
        % 
        for x = 1:128
            s = modified_fft_result_2d_final(x, :);
            t1 = max(abs(s(1:32)));
            t2 = std(abs(s(1:32)));
            mask = abs(s) > (t1 - t2) & abs(s) < (t1 + t2);
            mask(33:end) = false;  
            s(~mask) = 0;
            modified_fft_result_2d_final(x, :) = s;
        end

        c = [c; modified_fft_result_2d_final];
    end
    
    if receiver == 1
        c11 = c;
    elseif receiver == 2
        c22 = c;
    elseif receiver == 3
        c33 = c;
    end
end
%% 
figure;
imagesc(abs(c11)');
clim([0 2]);
colorbar;
title('Receiver 1');
 figure;
imagesc(abs(c22)');
clim([0 2]);
colorbar;
title('Receiver 2');
figure;
imagesc(abs(c33)');
clim([0 2]);
colorbar;
title('Receiver 3');
%% 
azimuth_data = zeros(12800, 64);
elevation_data = zeros(12800, 64);

for i = 1:100
    for x = 1:128
        
        phase_diff_azimuth = angle(c33((i-1)*128 + x,:)) - angle(c11((i-1)*128 + x,:));
        
        phase_diff_elevation = angle(c33((i-1)*128 + x,:)) - angle(c22((i-1)*128 + x,:));
        
        d_h = 2.5*10^-3; 
        d_v = 2.5*10^-3;   
        lambda = 5*10^-3;
        
         azimuth = (phase_diff_azimuth * lambda) / (2 * pi * d_h);
         azimuth = unwrap(azimuth);
         azimuth = max(min(azimuth, 1), -1);
        
        elevation = (phase_diff_elevation * lambda) / (2 * pi * d_v);
        elevation = unwrap(elevation);
        elevation = max(min(elevation, 1), -1);
        
        
        azimuth = asin(azimuth);
        elevation = asin(elevation);
        
        row_idx = (i-1)*128 + x;
        azimuth_data(row_idx, :) = azimuth;
        elevation_data(row_idx, :) = elevation;
    end
end
%% 
azimuth_data = azimuth_data(:,1:32);
elevation_data = elevation_data(:,1:32);

azimuth_deg = rad2deg(azimuth_data);
elevation_deg = rad2deg(elevation_data);

c22 = c22(:,1:32);
range_data = zeros(12800,32);

for i = 1:12800
    x = c22(i,:);
    for j = 1:32
        if(x(j) ~= 0)
            x(j) = j*7.5;
        end
    end
    range_data(i,:) = x;
end

% x_coords = zeros(size(range_data));
% y_coords = zeros(size(range_data));
% z_coords = zeros(size(range_data));
%% 
% Single Figure for Animation
mean_r1=[];
xx=[];yy=[];zz=[];


figure;
hold on;
grid on;
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
axis([-100 100 0 150 -40 40]);
title('3D Point Cloud Animation - Frame-wise (128 chirps per frame)');
% view(11,37);
% view(0,90);
trail_length =4;  

chirps_per_frame = 128;
num_frames = floor(size(range_data, 1) / chirps_per_frame);

scatter_handles = gobjects(trail_length, 1);

for i = 1:num_frames
    fx=[];fy=[];fz=[];
    for h = scatter_handles(:)'
        if isvalid(h)
            delete(h);
        end
    end

    if i == 1
        plot3(0, 0, 0, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
            'DisplayName', 'Radar Position');
    
        line([0 110], [0 0], [0 0], 'Color', 'k', 'LineStyle', ':');
        line([0 0], [0 110], [0 0], 'Color', 'k', 'LineStyle', ':');
        line([0 0], [0 0], [0 110], 'Color', 'k', 'LineStyle', ':');
    end
        frame_start = ((i-1) * chirps_per_frame) + 1;
        frame_end = i * chirps_per_frame;

        range_data1 = range_data(frame_start:frame_end, :);
        azimuth_data1 = azimuth_data(frame_start:frame_end, :);
        elevation_data1 = elevation_data(frame_start:frame_end, :);

      non_zero_r = range_data1(range_data1 ~= 0);
      non_zero_az=azimuth_data1(azimuth_data1~=0);
      non_zero_el=elevation_data1(elevation_data1~=0);

     mean_r = mean(non_zero_r);
     mean_az = mean(non_zero_az);
     mean_el= mean(non_zero_el);
mean_r1=[mean_r1;mean_r];

     cent_x = mean_r * cos(mean_el) * sin(mean_az);
     cent_y = mean_r * cos(mean_el) * cos(mean_az);
     cent_z = mean_r * sin(mean_el);
     xx=[xx;cent_x];yy=[yy;cent_y];zz=[zz;cent_z];

% for i = 1:size(range_data, 1)
for x = 1:128
    for y = 1:32
        if range_data1(x, y) ~= 0
            r = range_data1(x, y);
            az = azimuth_data1(x, y);
            el = elevation_data1(x, y);
            
            if ~isnan(az) && ~isnan(el)
                x_coords(x, y) = r * cos(el) * sin(az);
                y_coords(x, y) = r * cos(el) * cos(az);
                z_coords(x, y) = r * sin(el);
            end
            std_x=ceil(std(x_coords(x_coords~=0)));
             std_y=ceil(std(y_coords(y_coords~=0)));
              std_z=ceil(std(z_coords(z_coords~=0)));
            if ((cent_x-std_x)<x_coords(x, y) && x_coords(x, y)<(cent_x+std_x)) || ((cent_y-std_y)<y_coords(x, y) && y_coords(x, y)<(cent_y+std_y)) || ((cent_z-std_z)<z_coords(x, y) && z_coords(x, y)<(cent_z+std_z))
                fx=[fx,x_coords(x, y)];
                fy=[fy,y_coords(x, y)];
                fz=[fz,z_coords(x, y)];
            end
        end
    end
end

    for j = max(1, i-trail_length):i
      
        frame_x = fx(:);
        frame_y = fy(:);
        frame_z = fz(:);
        
        valid_points = frame_x ~= 0 & frame_y ~= 0 & frame_z ~= 0;
        alpha = (j - (i-trail_length)) / trail_length;
        
        scatter_idx = mod(j-1, trail_length) + 1;
        scatter_handles(scatter_idx) = scatter3(frame_x(valid_points), ...
            frame_y(valid_points), frame_z(valid_points), ...
            'filled', 'MarkerFaceAlpha', alpha);
        % clear scatter_handles;
    end
    title(['Frame: ' num2str(i)]);
    grid on;

    drawnow;
    pause(0.3);  
end

figure;plot(xx);title("X");
figure;plot(yy);title("Y");
figure;plot(zz);title("Z");
