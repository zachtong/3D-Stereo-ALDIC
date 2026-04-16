function [R,T] = GetRTMatrix(All_Points2D, Base_Points2D, Coordinate3D_ref)
% Interpolate the 3D coordinates corresponding to the three base points,
% then establish the specimen coordinate system, and return R (rotation)
% and T (translation).
%
% Uses scatteredInterpolant (linear) -- fast on large point sets.
% Previously used RBF interpolation (rbfcreate) which requires O(N^2)
% memory (for N=76104 quadtree points, a 46 GB dense linear system) and
% hours of CPU. For 3 smooth interior base points, linear interpolation
% differs from RBF by << 0.001 mm, well below the experimental noise floor.

Coordinate3D_ref = [Coordinate3D_ref{1}, Coordinate3D_ref{2}, Coordinate3D_ref{3}];

% Build interpolants once (fast: ~0.1 s for 76k points vs hours for RBF)
Fx = scatteredInterpolant(All_Points2D, Coordinate3D_ref(:,1), 'linear');
Fy = scatteredInterpolant(All_Points2D, Coordinate3D_ref(:,2), 'linear');
Fz = scatteredInterpolant(All_Points2D, Coordinate3D_ref(:,3), 'linear');

% Interpolate 3D positions at the 3 base points
Base_Points = zeros(size(Base_Points2D, 1), 3);
Base_Points(:, 1) = Fx(Base_Points2D);
Base_Points(:, 2) = Fy(Base_Points2D);
Base_Points(:, 3) = Fz(Base_Points2D);

% Build orthonormal specimen basis via Gram-Schmidt
new_x_axis = [Base_Points(2,:) - Base_Points(1,:)]'; new_x_axis = new_x_axis / norm(new_x_axis);
new_y_raw  = [Base_Points(3,:) - Base_Points(1,:)]'; new_y_raw  = new_y_raw  / norm(new_y_raw);
% Z = X x Y_raw (right-hand rule), then re-orthogonalize Y = Z x X
new_z_axis = cross(new_x_axis, new_y_raw);  new_z_axis = new_z_axis / norm(new_z_axis);
new_y_axis = cross(new_z_axis, new_x_axis); % unit length by construction

% Rotation matrix: columns are orthonormal basis vectors of specimen frame
R = [new_x_axis, new_y_axis, new_z_axis];
T = [Base_Points(1,:)];

end
