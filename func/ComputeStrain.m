% ==================================
% To compute strain on a uniform mesh
% ----------------------------------

% Plane fitting
%D = funDerivativeOp(M,N,DICpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
%FStrain = D*reshape(ULocal,length(ULocal),1);

% Compute strain method II: Use Plane Fitting method
prompt = 'What is your half window size: ';
Rad = input(prompt);
coefficients = cell(3,1);
for i = 2:size(FinalResult.Coordinates,1) % Num of images
    for j = 1:3 % U V W
        coefficients{i,j} = PlaneFit3(M, N, Rad, FinalResult.Coordinates{i,1}, FinalResult.Displacement{i,j});
    end
end
% [Ux Uy Uz]
% [Vx Vy Vz]
% [Wx Wy Wz]


%% Update infinitesimal strain to other finite strains
% FStrainFinite = FStrain;
% for tempi = 1:4:length(FStrain)
% 
%     % Obtain each component of def grad tensor
% 
%     switch DICpara.StrainType
%         case 0 % Infinitesimal stran
%             % Do nothing
%         case 1 % Eluerian strain
%             FStrainFinite(tempi) = 1/(1-dudx)-1;
%             FStrainFinite(tempi+3) = 1/(1-dvdy)-1;
%             FStrainFinite(tempi+2) = dudy/(1-dvdy);
%             FStrainFinite(tempi+1) = dvdx/(1-dudx);
%         case 2 % Green-Lagrangian strain: E=(C-I)/2
%             FStrainFinite(tempi) = 0.5*(dudx*2-dudx^2-dvdx^2);
%             FStrainFinite(tempi+3) = 0.5*(dvdy*2-dudy^2-dvdy^2);
%             FStrainFinite(tempi+2) = 0.5*(dudy+dvdx-dudx*dudy-dvdx*dvdy);
%             FStrainFinite(tempi+1) = 0.5*(dvdx+dudy-dudy*dudx-dvdy*dvdx);
%         case 3
%             disp('Press "Ctrl+C" to modify by yourself.'); pause;
%         otherwise
%             disp('Wrong strain ty                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           pe!');
%     end
% 
% end

%FStraintemp = FStrainFinite(temp3);
%FStrainWorld = FStraintemp; FStrainWorld(2:4:end) = -FStrainWorld(2:4:end); FStrainWorld(3:4:end) = -FStrainWorld(3:4:end);
