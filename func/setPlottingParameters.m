function DICpara = setPlottingParameters(DICpara, mode)
%SETPLOTTINGPARAMETERS Configure plotting parameters for displacement and strain
%
% INPUTS:
%   DICpara - DIC parameter structure
%   mode    - 'auto', 'custom', or 'interactive'
%
% OUTPUTS:
%   DICpara - Updated DIC parameter structure with plotting settings

if nargin < 2
    mode = 'auto';  % Default to auto mode
end



%% Set parameters based on mode
switch lower(mode)
    case 'auto'
        fprintf('Setting automatic caxis scaling...\n');
        
        % All caxis set to 'auto'
        caxis_fields = {'caxis_u', 'caxis_v', 'caxis_w', 'caxis_mag', ...
                        'caxis_exx', 'caxis_exy', 'caxis_eyy', ...
                        'caxis_principal_max', 'caxis_principal_min', ...
                        'caxis_max_shear', 'caxis_vonMises', ...
                        'caxis_dwdx', 'caxis_dwdy'};
        
        for i = 1:length(caxis_fields)
            DICpara.(caxis_fields{i}) = 'auto';
        end
        
        % Default colormaps
        DICpara.colormap_u = 'turbo';
        DICpara.colormap_v = 'turbo';
        DICpara.colormap_w = 'turbo';
        DICpara.colormap_mag = 'turbo';
        DICpara.colormap_exx = 'turbo';
        DICpara.colormap_exy = 'turbo';
        DICpara.colormap_eyy = 'turbo';
        DICpara.colormap_vonMises = 'turbo';
        
    case 'custom'
        fprintf('Setting custom caxis values...\n');
        
        % Custom displacement caxis
        DICpara.caxis_u = [-0.01, 0.01];
        DICpara.caxis_v = [-0.01, 0.01];
        DICpara.caxis_w = [-0.01, 0.01];
        DICpara.caxis_mag = [0, 0.02];
        
        % Custom strain caxis
        DICpara.caxis_exx = [-0.01, 0.01];
        DICpara.caxis_exy = [-0.01, 0.01];
        DICpara.caxis_eyy = [-0.01, 0.01];
        DICpara.caxis_principal_max = [-0.01, 0.01];
        DICpara.caxis_principal_min = [-0.01, 0.01];
        DICpara.caxis_max_shear = [0, 0.02];
        DICpara.caxis_vonMises = [0, 0.05];
        DICpara.caxis_dwdx = [-0.01, 0.01];
        DICpara.caxis_dwdy = [-0.01, 0.01];
        
        % Custom colormaps
        DICpara.colormap_u = 'turbo';
        DICpara.colormap_v = 'turbo';
        DICpara.colormap_w = 'turbo';
        DICpara.colormap_mag = 'turbo';
        DICpara.colormap_exx = 'turbo';
        DICpara.colormap_exy = 'turbo';
        DICpara.colormap_eyy = 'turbo';
        DICpara.colormap_vonMises = 'turbo';
        
    case 'interactive'
        fprintf('Interactive caxis configuration...\n');
        
        % Ask user for each parameter
        response = questdlg('Use automatic or custom caxis scaling?', ...
                           'Caxis Configuration', 'Auto', 'Custom', 'Auto');
        
        if strcmp(response, 'Auto')
            DICpara = setPlottingParameters(DICpara, 'auto');
        else
            DICpara = setPlottingParameters(DICpara, 'custom');
            
            % Allow user to modify specific values
            fprintf('Current custom values set. You can modify them in the workspace.\n');
            fprintf('Example: DICpara.caxis_u = [-0.005, 0.005];\n');
        end
        
    otherwise
        error('setPlottingParameters: Unknown mode "%s". Use "auto", "custom", or "interactive"', mode);
end

fprintf('Plotting parameters configured in "%s" mode.\n', mode);

end
