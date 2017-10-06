function hough_img = generateHoughAccumulator(img, theta_num_bins, rho_num_bins)
    [y_max,x_max] = size(img);
    diagonal = (x_max.^2+y_max.^2).^(0.5);
    accumulator = zeros(rho_num_bins, theta_num_bins);
    for y = 1:y_max
        for x = 1:x_max
            if img(y,x)>50
                accumulator = accumulator + find_theta_rho(...
                    y,x,diagonal,theta_num_bins, rho_num_bins);
            end
        end
    end
    maximum = max(accumulator(:));
    hough_img = accumulator ./maximum  .* 255;
    
%     fh = figure; subplot(2, 2, 1); imshow(img); title('Original');
%     subplot(2, 2, 2); imshow(hough_img);
end

function binary_accumulator=find_theta_rho(y,x,diagonal,theta_num_bins, rho_num_bins)
    binary_accumulator = zeros(rho_num_bins, theta_num_bins);
    A = (x.^2+y.^2).^(0.5);
    sine = y ./ A;
    phase = asin(sine);
    if x > 0
        phase = pi - phase;
    end
    for theta_bin=1:theta_num_bins
%         for j=1:100
            theta = pi ./ theta_num_bins .* (theta_bin - 1 + 0.5);
            rho = A .* sin(phase + theta);
            rho_bin = ceil((rho + diagonal) ./ (diagonal .* 2 ./ rho_num_bins));
            binary_accumulator(rho_bin, theta_bin) = 1;
%         end
    end
end