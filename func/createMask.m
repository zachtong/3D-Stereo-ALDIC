function maskResult = createMask(Image)
    % 创建一个全局变量来存储最终的mask
    global MASK_RESULT
    MASK_RESULT = [];
    
    Image = Image';
    % 创建figure & axes
    hFig = figure('Name','Polygon Mask Tool',...
        'NumberTitle','off',...
        'MenuBar','none',...
        'ToolBar','none');
    hAx = axes('Parent', hFig);
    
    % ------------ 关键1：获取原图的句柄 hImg ------------
    hImg = imshow(Image, [], 'Parent', hAx);
    hold(hAx, 'on');

    % 初始化数据    
    % 创建legend的模拟对象
    hROI = patch('XData', [], 'YData', [], 'FaceColor', 'r', 'EdgeColor', 'r', 'Visible', 'off');
    hBG = patch('XData', [], 'YData', [], 'FaceColor', 'g', 'EdgeColor', 'g', 'Visible', 'off');
    
    % 添加legend
    legend([hROI, hBG], 'ROI', 'Background', 'Location', 'northeast');
    
    data.addMask = false(size(Image));
    data.cutMask = false(size(Image));
    data.hMaskOverlay = [];
    data.hOriginalImage = hImg;
    guidata(hFig, data);
    
    % 按钮 - Add Polygon (green)
    uicontrol('Parent', hFig,...
        'Style','pushbutton',...
        'String','Add Polygon',...
        'Position',[60 20 100 30],...
        'Callback', @addPolygonCallback);
    
    % 按钮 - Cut Polygon (Red)
    uicontrol('Parent', hFig,...
        'Style','pushbutton',...
        'String','Cut Polygon',...
        'Position',[160 20 100 30],...
        'Callback', @cutPolygonCallback);
    
    % 按钮 - Done
    uicontrol('Parent', hFig,...
        'Style','pushbutton',...
        'String','Done',...
        'Position',[260 20 100 30],...
        'Callback', @doneCallback);
    
    % 等待直到Done按钮被点击且MASK_RESULT被设置
    while isempty(MASK_RESULT)
        pause(0.1);
        if ~isvalid(hFig)  % 如果窗口被关闭
            MASK_RESULT = false(size(Image));  % 返回空mask
            break;
        end
    end
    
    % 清理并返回结果
    maskResult = MASK_RESULT;
    clear global MASK_RESULT
    close(hFig);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 回调函数
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addPolygonCallback(~, ~)
        data = guidata(hFig);
        % 在 hAx 上画一个红色多边形
        hPoly = drawpolygon('Parent', hAx, ...
            'LineWidth',1, 'Color','r', 'FaceAlpha',0.1);
        if ~isempty(hPoly) && isvalid(hPoly)
            % ------------ 关键2：显式使用原图句柄 createMask ------------
            newMask = createMask(hPoly, data.hOriginalImage);
            data.addMask = data.addMask | newMask;
            data = updateOverlay(data);
            guidata(hFig, data);
        end
    end

    function cutPolygonCallback(~, ~)
        data = guidata(hFig);
        % 在 hAx 上画一个绿色多边形
        hPoly = drawpolygon('Parent', hAx, ...
            'LineWidth',1, 'Color','g', 'FaceAlpha',0.1);
        if ~isempty(hPoly) && isvalid(hPoly)
            % 同样必须带 hOriginalImage
            newMask = createMask(hPoly, data.hOriginalImage);
            data.cutMask = data.cutMask | newMask;
            data = updateOverlay(data);
            guidata(hFig, data);
        end
    end

    function doneCallback(~, ~)
        data = guidata(hFig);
        finalMask = data.addMask & ~data.cutMask;
        MASK_RESULT = finalMask';  % 设置全局变量
    end

    % 叠加显示
    function data = updateOverlay(data)
        finalMask = data.addMask & ~data.cutMask;
        if ~isempty(data.hMaskOverlay) && isvalid(data.hMaskOverlay)
            delete(data.hMaskOverlay);
        end
        % 叠加一个红色
        maskRGB = cat(3, ones(size(finalMask)), zeros(size(finalMask)), zeros(size(finalMask)));
        maskRGB(~finalMask) = 0;
        % 在同一个axes上叠加
        data.hMaskOverlay = imshow(maskRGB, 'Parent', hAx);
        set(data.hMaskOverlay, 'AlphaData', 0.3);
    end
end