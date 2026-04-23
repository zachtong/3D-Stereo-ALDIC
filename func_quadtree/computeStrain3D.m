function [exx, eyy, exy, e1, e2, maxShear, vonMises, dwdx, dwdy] = computeStrain3D(coefficients)
%COMPUTESTRAIN3D  Compute Green-Lagrange strain from displacement-gradient coefficients.
%
%   [exx, eyy, exy, e1, e2, maxShear, vonMises, dwdx, dwdy] = computeStrain3D(coefficients)
%
%   coefficients is an Nx1 cell array where each cell is a 3x3 matrix:
%       [U,x  V,x  W,x
%        U,y  V,y  W,y
%        U,z  V,z  W,z]
%   as produced by PlaneFit3_Quadtree. Rows = partial-derivative axis
%   (x, y, z), columns = displacement component (U, V, W).
%
%   Returns per-node strain components (each Nx1):
%       exx, eyy, exy  - Green-Lagrange strain tensor entries in local frame
%       e1, e2         - principal strains (major / minor)
%       maxShear       - maximum shear strain
%       vonMises       - equivalent von Mises strain
%       dwdx, dwdy     - out-of-plane slope (for plotting / diagnostics)
%
%   The deformation gradient F uses all 9 coefficient entries:
%
%       F = [1+u_x   u_y    u_z;
%             v_x  1+v_y    v_z;
%             w_x    w_y  1+w_z]
%
%   Under 'Local' strain coordinates (PlaneFit3_Quadtree default), the
%   z-row of coefficients is zero, so u_z = v_z = w_z = 0 and F matches
%   the simpler acc-mode form used by older code.

n = size(coefficients, 1);
[exx, eyy, exy, dwdx, dwdy] = deal(zeros(n, 1));

for i = 1:n
    u_x = coefficients{i,1}(1,1); u_y = coefficients{i,1}(2,1); u_z = coefficients{i,1}(3,1);
    v_x = coefficients{i,1}(1,2); v_y = coefficients{i,1}(2,2); v_z = coefficients{i,1}(3,2);
    w_x = coefficients{i,1}(1,3); w_y = coefficients{i,1}(2,3); w_z = coefficients{i,1}(3,3);
    F = [1+u_x, u_y,   u_z;
         v_x,   1+v_y, v_z;
         w_x,   w_y,   1+w_z];
    E = 0.5 * (F' * F - eye(3));
    exx(i) = E(1,1);
    eyy(i) = E(2,2);
    exy(i) = E(1,2);
    dwdx(i) = w_x;
    dwdy(i) = w_y;
end

maxShear = sqrt((0.5*(exx - eyy)).^2 + exy.^2);
e1       = 0.5*(exx + eyy) + maxShear;
e2       = 0.5*(exx + eyy) - maxShear;
vonMises = sqrt(e1.^2 + e2.^2 - e1.*e2 + 3*maxShear.^2);

end
