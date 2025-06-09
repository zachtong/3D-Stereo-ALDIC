function Base_Points_2D = GetBasePoints(imageLeft)
if any(imageLeft(:)>255)
imshow(uint16(imageLeft'));
else
    imshow(uint8(imageLeft'));
end
title('Click 3 points as the O,X and Y axis directions.','fontweight','normal','fontsize',16);

Point_O = ginput(1);
fprintf('Coordinates of point O is (%4.3f,%4.3f)\n',Point_O(1), Point_O(2));
Point_X = ginput(1);
fprintf('Coordinates of point X is (%4.3f,%4.3f)\n',Point_X(1), Point_X(2));
Point_Y = ginput(1);
fprintf('Coordinates of point Y is (%4.3f,%4.3f)\n',Point_Y(1), Point_Y(2));

Base_Points_2D = [Point_O;Point_X;Point_Y];

end