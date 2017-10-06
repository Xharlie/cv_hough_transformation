function cropped_line_img = lineSegmentFinder(orig_img, hough_img, hough_threshold)
    fh1 = figure();
    imshow(orig_img);
    [y_max,x_max] = size(orig_img);
    diagonal = (x_max.^2+y_max.^2).^(0.5);
    [rho_num_bins, theta_num_bins] = size(hough_img);
    for rho_bin = 1:rho_num_bins
        for theta_bin = 1:theta_num_bins
            if hough_img(rho_bin, theta_bin) >= hough_threshold
                rho = 2 * diagonal ./ rho_num_bins .* (rho_bin - 0.5) - diagonal;
                theta = pi ./ theta_num_bins .* (theta_bin - 0.5);
                draw_line_segment(orig_img, rho, theta);
            end
        end
    end
    cropped_line_img = saveAnnotatedImg(fh1);
    delete(fh1); 
end

function draw_line_segment(orig_img, rho, theta)
   % Get all pixels associated with Hough transform cell.
   [ys, xs] = corresponding_pixels(orig_img, theta, rho);
   if (size(ys,1)==0)
       return
   end
   % Rotate the pixel locations about (1,1) so that they lie
   % approximately along a vertical line.
   T = [cos(-theta) sin(-theta); -sin(-theta) cos(-theta)];
   xy = [xs + (rho ./ sin(theta)), ys] * T;
   x = sort(xy(:,1));    
%    hold on;    
%    plot(xs ,ys, 'r*');
   % Find the gaps larger than the threshold.
   diff_x = [diff(x); Inf];
   idx = [0; find(diff_x > 2)];
   for p = 1:length(idx) - 1
      x1 = x(idx(p) + 1); x2 = x(idx(p + 1));
      linelength = x2 - x1;
      if linelength >= 20
         point1 = [x1 0]; point2 = [x2 0];
         % Rotate the end-point locations back to the original
         % angle.
         Tinv = inv(T);
         anchor = [rho ./ sin(theta), 0];
         point1 = point1 * Tinv - anchor; point2 = point2 * Tinv - anchor;
         line_array = [point1; point2];
         hold on;
         plot(line_array(:, 1), line_array(:, 2),...
            'LineWidth',1, 'Color', [0, 1, 0]);
      end
   end
end

function [ys, xs] = corresponding_pixels(ori_img, theta, rho)
    [y, x, val] = find(ori_img);
    rho_point = y.*cos(theta) - x.*sin(theta);
    index = find(abs(rho_point - rho)<= 1);
    y = y(index); x = x(index);
    index_2 = find(check(ori_img, y,x) > 0);
    ys = y(index_2); xs = x(index_2);
end

function flag=check(ori_img, ylist, xlist) 
    flag = zeros(size(ylist));
    [y_max, x_max] = size(ori_img);
    for i = 1:size(flag)
        y = ylist(i);
        x = xlist(i);
        if (y + 1 <= y_max) && (y - 1 >= 1) && (x + 1 <= x_max) && (x - 1 >= 1) 
            com = [ori_img(y-1,x-1), ori_img(y+1,x+1), ...
                ori_img(y-1,x+1), ori_img(y-1,x+1)];
            if max(com) - min(com) >= 4.5
                flag(i)=1;
            end
        end
    end
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


