function matchedPairs = organizeMatchedPairds_quadtree_ST2(RD_L,RD_R,DICpara)
%% Acc. mode
if DICpara.DICIncOrNot == 0
    % Initializition
    imagePairNum = size(RD_L.ResultDisp,1) + 1;
    matchedPairs = cell(imagePairNum,1);
    
    matchedPairs{1,1} = [RD_L.ResultFEMesh{1, 1}.coordinatesFEM,RD_R.Coordinates_corr];

    %%%%% Since second frame %%%%%
    for i = 1:imagePairNum-1
        matchedPairs{i+1,1}(:,1:2) = RD_L.ResultFEMesh{1,1}.coordinatesFEM + reshape(RD_L.ResultDisp{i}.U',[2,size(RD_L.ResultDisp{i}.U,1)/2])' ;
        matchedPairs{i+1,1}(:,3:4) = RD_L.ResultFEMesh{1,1}.coordinatesFEM + reshape(RD_R.ResultDisp{i}.U',[2,size(RD_R.ResultDisp{i}.U,1)/2])';
    end

%% Inc. mode
elseif DICpara.DICIncOrNot == 1
    %  TBD!!

    % % Initializition
    % imagePairNum = size(RD_L.ResultDisp,1);
    % matchedPairs = cell(imagePairNum,1);
    % 
    % matchedPairs{1,1} = [RD_L.ResultFEMesh{1, 1}.coordinatesFEM,RD_R.ResultFEMesh_corr{1, 1}];
    % 
    % %%%%% Since second frame %%%%%
    % for i = 1:imagePairNum-1
    %     matchedPairs{i+1,1}(:,1:2) = RD_L.ResultFEMesh{1,1}.coordinatesFEM + [RD_L.ResultDisp{i+1}.U_acc(1:2:end),RD_L.ResultDisp{i+1}.U_acc(2:2:end)];
    %     matchedPairs{i+1,1}(:,3:4) = RD_L.ResultFEMesh{1,1}.coordinatesFEM + [RD_R.ResultDisp{i+1}.U_acc(1:2:end),RD_R.ResultDisp{i+1}.U_acc(2:2:end)];
    % end


end