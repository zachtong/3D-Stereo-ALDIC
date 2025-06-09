function [X_coords, Y_coords] = generateGridMatrices(DICpara)
    % check valid
    % if mod(x_end - x_start, step) ~= 0 || mod(y_end - y_start, step) ~= 0
    %     error('The provided range cannot be evenly divided by the step size.');
    % end

    x_start = DICpara.gridxyROIRange.gridx(1);
    x_end = DICpara.gridxyROIRange.gridx(2);
    y_start = DICpara.gridxyROIRange.gridy(1);
    y_end = DICpara.gridxyROIRange.gridy(2);
    step = DICpara.winstepsize;
    
    x = x_start:step:x_end;
    y = y_start:step:y_end;

    
    X_coords = zeros(length(y), length(x));
    Y_coords = zeros(length(y), length(x));

    
    for j = 1:length(y)
        for i = 1:length(x)
            X_coords(j, i) = x(i);
            Y_coords(j, i) = y(j);
        end
    end
end