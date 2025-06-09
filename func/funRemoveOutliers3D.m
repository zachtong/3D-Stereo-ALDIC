function Displacement = funRemoveOutliers3D(Displacement,M,N,x0,y0)

fprintf('Do you clear bad points by setting upper/lower bounds? (0-yes; 1-no)  \n')
prompt = 'Input here: ';
ClearBadPointsOrNot = input(prompt);

while ClearBadPointsOrNot == 0

    u = Displacement{1,1}; v = Displacement{1,2}; w = Displacement{1,3};

    prompt = 'The displacement is almost uniform? \n1: Manually setting [upper_bound, upper_bound](by default) \n2: Remove outliers automatically (only for uniform displacement field)\n Input:';
    RemoveMode = input(prompt);

    if RemoveMode == 2

        % check
        % addon_table = matlab.addons.installedAddons;
        % InstallRequiredAddonOrNot = any(strcmp(addon_table.Name, 'Point cloud tools for Matlab'));
        % if InstallRequiredAddonOrNot == 0
        %     disp('Please install Matlab Addon "Point cloud tools for Matlab".');
        % else

        % x0 = x0';
        % y0 = y0';
        % x0_reshaped = x0(:);
        % y0_reshaped = y0(:);
        % fittingFcn = @(x) [x ones(size(x, 1), 1)] \ zeros(size(x, 1), 1);
        % distFcn = @(model, x) abs(x * model(1:3) + model(4)) / norm(model(1:3));
        % [model, inliers] = ransac([x0_reshaped,y0_reshaped,u], fittingFcn, distFcn, 23, 100);

        u_mean = mean(u); u_std = std(u);
        v_mean = mean(v); v_std = std(v);
        w_mean = mean(w); w_std = std(w);
        [row1,~] = find( u > u_mean+2*u_std );
        [row2,~] = find( u < u_mean-2*u_std );
        [row3,~] = find( v > v_mean+2*v_std );
        [row4,~] = find( v < v_mean-2*v_std );
        [row5,~] = find( w > w_mean+2*w_std );
        [row6,~] = find( w < w_mean-2*w_std );
        row = [row1; row2; row3; row4; row5; row6];

        for tempi = 1:length(row)
            u(row(tempi),1)=NaN;
            v(row(tempi),1)=NaN;
            w(row(tempi),1)=NaN;
        end

        u = inpaint_nans(u,4);
        v = inpaint_nans(v,4);
        w = inpaint_nans(w,4);

    else
        prompt = 'What is the upper bound for x-disp? Input: ';
        upperbound = input(prompt);
        [row1,~] = find(u>upperbound);
        prompt =  'What is the lower bound for x-disp? Input: ';
        lowerbound = input(prompt);
        [row2,~] = find(u<lowerbound);
        prompt = 'What is the upper bound for y-disp? Input: ';
        upperbound = input(prompt);
        [row3,~] = find(v>upperbound);
        prompt = 'What is the lower bound for y-disp? Input: ';
        lowerbound = input(prompt);
        [row4,~] = find(v<lowerbound);
        prompt = 'What is the upper bound for z-disp? Input: ';
        upperbound = input(prompt);
        [row5,~] = find(w>upperbound);
        prompt = 'What is the lower bound for z-disp? Input: ';
        lowerbound = input(prompt);
        [row6,~] = find(w<lowerbound);
        row = [row1; row2; row3; row4; row5; row6];

        for tempi = 1:length(row)
            u(row(tempi),1)=NaN;
            v(row(tempi),1)=NaN;
            w(row(tempi),1)=NaN;
        end

        u = inpaint_nans(u,4);
        v = inpaint_nans(v,4);
        w = inpaint_nans(w,4);
    end
    % --------------------------------------
    close all;
    u(HolePtInd) = nan;  v(HolePtInd) = nan; w(HolePtInd) = nan; % Set points where is a hole as nans
    figure; surf(u); colorbar; title('Displacement u','fontweight','normal');
    colormap(jet);  set(gca,'fontsize',18); set(gcf,'color','w'); a = gca; a.TickLabelInterpreter = 'latex';
    figure; surf(v); colorbar; title('Displacement v','fontweight','normal');
    colormap(jet);  set(gca,'fontsize',18); set(gcf,'color','w'); a = gca; a.TickLabelInterpreter = 'latex';
    figure; surf(w); colorbar; title('Displacement w','fontweight','normal');
    colormap(jet);  set(gca,'fontsize',18); set(gcf,'color','w'); a = gca; a.TickLabelInterpreter = 'latex';
end