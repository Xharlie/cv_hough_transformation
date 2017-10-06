function line_detected_img = lineFinder(orig_img, hough_img, hough_threshold, local_max_diameter)
    rho_diameter = local_max_diameter(1);
    theta_diameter = local_max_diameter(2);
    fh1 = figure();
    imshow(orig_img);
    [y_max,x_max] = size(orig_img);
    diagonal = (x_max.^2+y_max.^2).^(0.5);
    [rho_num_bins, theta_num_bins] = size(hough_img);
    candidate_array=[];
    for rho_bin = 1:rho_num_bins
        for theta_bin = 1:theta_num_bins
            if hough_img(rho_bin, theta_bin) >= hough_threshold
                candidate_array = [candidate_array; [rho_bin, theta_bin]];
            end
        end
    end
    local_max=local_filter(candidate_array, rho_diameter, theta_diameter, hough_img);
    idx = find(local_max>0);
    rho_bin = candidate_array(idx,1);
    theta_bin = candidate_array(idx,2);
    rho = 2 * diagonal ./ rho_num_bins .* (rho_bin - 0.5) - diagonal;
    theta = pi ./ theta_num_bins .* (theta_bin - 0.5);  
    for i = 1:(size(rho,1))
        draw_line(rho(i), theta(i), y_max, x_max)
    end 
    line_detected_img = saveAnnotatedImg(fh1);
    delete(fh1);
end

function valid_array=local_filter(candidate_array, rho_diameter, theta_diameter, hough_img)
    valid_array= ones(size(candidate_array,1),1);
    for i=1:(size(candidate_array,1)-1)
        if valid_array(i) == 0
            continue;
        end
        for j=i+1:size(candidate_array,1)
            if abs(candidate_array(i,1) - candidate_array(j,1)) < rho_diameter ...
                && abs(candidate_array(i,2) - candidate_array(j,2)) < theta_diameter
                if hough_img(candidate_array(i,1),candidate_array(i,2)) >= ...
                    hough_img(candidate_array(j,1),candidate_array(j,2))
                    valid_array(j) = 0;
                else 
                    valid_array(i) = 0;
                    break;
                end
            end
        end
    end
end

function draw_line(rho, theta, y_max, x_max)
    hold on;
    line_array = [];
    point_x = 0;
    point_y = rho ./ cos(theta);
%     if point_y >= -0.01 && point_y <= y_max
        line_array = [point_x, point_y];
%     end
    point_x = x_max;
    point_y = (rho  +  x_max .* sin(theta)) ./ cos(theta);
%     if point_y >= -0.01 && point_y <= y_max
        line_array = [line_array; [point_x,point_y]];
%     end
%     point_x = -1 .* rho ./ sin(theta)
%     point_y = 0
%     if point_x >= -0.01 && point_x <= x_max
%         line_array = [line_array; [point_x,point_y]];
%     end
%     point_x = (y_max * cos(theta)  - rho) ./ sin(theta)
%     point_y = y_max
%     if point_x >= -0.01 && point_x <= x_max
%         line_array = [line_array; [point_x,point_y]];
%     end
    plot(line_array(:, 1), line_array(:, 2),...
        'LineWidth',1, 'Color', [0, 1, 0]);
end

function annotated_img = saveAnnotatedImg(fh)
figure(fh); % Shift the focus back to the figure fh

% The figure needs to be undocked
set(fh, 'WindowStyle', 'normal');

% The following two lines just to make the figure true size to the
% displayed image. The reason will become clear later.
img = getimage(fh);
truesize(fh, [size(img, 1), size(img, 2)]);

% getframe does a screen capture of the figure window, as a result, the
% displayed figure has to be in true size. 
frame = getframe(fh);
frame = getframe(fh);
pause(0.5); 
% Because getframe tries to perform a screen capture. it somehow 
% has some platform depend issues. we should calling
% getframe twice in a row and adding a pause afterwards make getframe work
% as expected. This is just a walkaround. 
annotated_img = frame.cdata;
end