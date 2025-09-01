function [R,T] = GetRTMatrix(All_Points2D, Base_Points2D,Coordinate3D_ref)
% Interpolate the 3D coordinates corresponding to the three base points,
% then establish the coordinate system, and find the R matrix

Coordinate3D_ref = [Coordinate3D_ref{1},Coordinate3D_ref{2},Coordinate3D_ref{3}];

% Method 1: rbf interp (Global interpolation, high memory requirement)
% interpolation OXY corr in World_corr
rbfInterpX = rbfcreate(All_Points2D',(Coordinate3D_ref(:,1))');
rbfInterpY = rbfcreate(All_Points2D',(Coordinate3D_ref(:,2))');
rbfInterpZ = rbfcreate(All_Points2D',(Coordinate3D_ref(:,3))');

% interpolating
Base_Points(1,:) = rbfinterp(Base_Points2D',rbfInterpX);
Base_Points(2,:) = rbfinterp(Base_Points2D',rbfInterpY);
Base_Points(3,:) = rbfinterp(Base_Points2D',rbfInterpZ);
Base_Points = Base_Points';

new_x_axis = [Base_Points(2,:)-Base_Points(1,:)]'; new_x_axis = new_x_axis / norm(new_x_axis);
new_y_axis = [Base_Points(3,:)-Base_Points(1,:)]'; new_y_axis = new_y_axis / norm(new_y_axis);
new_z_axis = cross(new_x_axis,new_y_axis);
% Rotation matrix between 2 coor. sys.
R = [new_x_axis,new_y_axis,new_z_axis];
T = [Base_Points(1,:)];


% Method 2: Local linear interpolation (faster)
% Coordinate3D_ref = [Coordinate3D_ref{1}, Coordinate3D_ref{2}, Coordinate3D_ref{3}];
% 
% Fx = scatteredInterpolant(All_Points2D, Coordinate3D_ref(:,1), 'linear');
% Fy = scatteredInterpolant(All_Points2D, Coordinate3D_ref(:,2), 'linear');
% Fz = scatteredInterpolant(All_Points2D, Coordinate3D_ref(:,3), 'linear');
% 
% Base_Points = zeros(size(Base_Points2D, 1), 3);
% Base_Points(:, 1) = Fx(Base_Points2D);
% Base_Points(:, 2) = Fy(Base_Points2D);
% Base_Points(:, 3) = Fz(Base_Points2D);
% 
% new_x_axis = Base_Points(2,:) - Base_Points(1,:); 
% new_x_axis = new_x_axis' / norm(new_x_axis); 
% 
% new_y_axis = Base_Points(3,:) - Base_Points(1,:); 
% new_y_axis = new_y_axis' / norm(new_y_axis); 
% 
% new_z_axis = cross(new_x_axis, new_y_axis);
% new_z_axis = new_z_axis / norm(new_z_axis); 
% 
% 
% new_y_axis = cross(new_z_axis, new_x_axis);
% 
% R = [new_x_axis, new_y_axis, new_z_axis];
% T = Base_Points(1,:); 
end