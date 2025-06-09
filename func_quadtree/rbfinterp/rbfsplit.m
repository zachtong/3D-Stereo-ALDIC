function interpolatedValues = rbfsplit(Coordinates,Values,interpolatedCoor,maxPoints,minSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 使用quadtree把大规模数据集分割成不同的区域,每个子叶内的插值点独立进行插值，
% 避免大规模rbfcoeff的运算
% interpolatedValues, Coordinates, Values, interpolatedCoor 都是列排列


% 四叉树划分参数
% maxPoints = 20000; % 每个区域内数据点数量的上限
% minSize = 100; % 区域大小的下限 default取决于图像大小和stepsize

% 构建四叉树
tree = buildQuadtree(Coordinates, Values, maxPoints, minSize);

% 插值
xi = interpolatedCoor(:,1);
yi = interpolatedCoor(:,2);
interpolatedValues = interpolateWithQuadtree(tree, xi,yi);
end



% 四叉树构建函数
function tree = buildQuadtree(coords, values, maxPoints, minSize, bounds, isUp, isDown, isLeft, isRight)
    a = zeros(size(coords,1),1);

if nargin < 5
    bounds = [min(coords); max(coords)];
    isUp = true;
    isDown = true;
    isLeft = true;
    isRight = true;
end

tree.coords = coords;
tree.values = values;
tree.bounds = bounds;
tree.children = {};
tree.isUp = isUp;
tree.isDown = isDown;
tree.isLeft = isLeft;
tree.isRight = isRight;

if size(coords, 1) > maxPoints && all(bounds(2,:) - bounds(1,:) > minSize)
    % 划分为四个象限
    center = (bounds(1,:) + bounds(2,:)) / 2;
    for i = 1:4
        childBounds = bounds;
        childBounds(1 + mod(i-1,2), 1) = center(1);
        childBounds(2 - floor((i-1)/2), 2) = center(2);
        idx = (coords(:,1) >= childBounds(1,1) & coords(:,1) <= childBounds(2,1)) & ...
            (coords(:,2) >= childBounds(1,2) & coords(:,2) <= childBounds(2,2));
        a = a + idx;
        if any(idx)   % 对于需要外插的节点，需要通过isUp Down Left Righ来判断位于最外侧的子叶所在的相对位置
            childIsUp = isUp && (i == 1 || i == 2);
            childIsDown = isDown && (i == 3 || i == 4);
            childIsLeft = isLeft && (i == 2 || i == 4);
            childIsRight = isRight && (i == 1 || i == 3);
            tree.children{i} = buildQuadtree(coords(idx,:), values(idx), maxPoints, minSize, childBounds, ...
                childIsUp, childIsDown, childIsLeft, childIsRight);
        end
    end
end

end

% 插值函数
function zi = interpolateWithQuadtree(tree,xi,yi)
zi = zeros(size(xi));

% 对每个叶节点进行处理
if ~isempty(tree.children) % Split the data set
    for i = 1:numel(tree.children)
        if isempty(tree.children{i}.children)
            coords = tree.children{i}.coords;
            values = tree.children{i}.values;
            bounds = tree.children{i}.bounds;
            isUp = tree.children{i}.isUp;
            isDown = tree.children{i}.isDown;
            isLeft = tree.children{i}.isLeft;
            isRight = tree.children{i}.isRight;

            % 找到属于当前叶节点的待插值点
            % 对于需要外插的节点，需要通过isUp Down Left Righ来判断位于最外侧的子叶所在的相对位置
            idx = (xi > bounds(1,1)*(1-isLeft) & xi <= bounds(2,1)*(1-isRight) + isRight*10000 ) & ...
                (yi > bounds(1,2)*(1-isUp) & yi <= bounds(2,2)*(1-isDown) + isDown*10000 );

            if any(idx)
                % 邻域扩展
                expandFactor = 0.3; % 扩展因子
                expandDist = max(max(coords) - min(coords)) * expandFactor;
                coordsExpand = coords;
                valuesExpand = values;

                for j = 1:4 % 从同一级的其他叶中找到拓展的节点
                    if ~isempty(tree.children{j}) && j ~= i
                        coordsChild = tree.children{j}.coords;
                        valuesChild = tree.children{j}.values;

                        idxChild = (coordsChild(:,1) >= bounds(1,1) - expandDist & coordsChild(:,1) <= bounds(2,1) + expandDist) & ...
                            (coordsChild(:,2) >= bounds(1,2) - expandDist & coordsChild(:,2) <= bounds(2,2) + expandDist);

                        coordsExpand = [coordsExpand; coordsChild(idxChild,:)];
                        valuesExpand = [valuesExpand; valuesChild(idxChild)];
                    end
                end

                % 去掉重复的行
                [coordsExpand,valuesExpand] = removeDuplicateRows(coordsExpand,valuesExpand);

                try
                    % 计算RBF系数
                    rbfOptions = rbfcreate(coordsExpand', valuesExpand','RBFFunction', 'thinplate');
                    % 对这些待插值点进行插值
                    zi(idx) = rbfinterp([xi(idx) yi(idx)]', rbfOptions);
                catch
                    disp('miss!');
                end
            end
        else
            % 递归处理子节点
            zi = zi + interpolateWithQuadtree(tree.children{i}, xi, yi);
        end
    end
else % Non quadtree, no split data set
    coords = tree.coords;
    values = tree.values;
    try
        % 计算RBF系数
        rbfOptions = rbfcreate(coords', values');
        % 对这些待插值点进行插值
        zi = rbfinterp([xi yi]', rbfOptions);
    catch
    end
end
end


function [uniqueMatrix1,uniqueMatrix2] = removeDuplicateRows(matrix1,matrix2)
% 使用 unique 函数找出唯一的行
[uniqueRows, ~, ic] = unique(matrix1, 'rows', 'stable');

% 找出重复行的索引
duplicateIndices = accumarray(ic, 1) > 1;

% 只保留唯一的行
uniqueMatrix1 = uniqueRows(~duplicateIndices, :);
uniqueMatrix2 = matrix2(~duplicateIndices, :);

end

