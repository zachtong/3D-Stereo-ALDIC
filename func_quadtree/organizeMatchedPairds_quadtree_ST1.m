function matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L,RD_R,DICpara)

%% Acc. & Inc. mode
% Initializition
imagePairNum = size(RD_L.ResultDisp,1) + 1;
matchedPairs = cell(imagePairNum,1);

matchedPairs{1,1} = [RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, RD_R.Coordinates_corr];

%%%%% Since second frame %%%%%
for i = 1:imagePairNum-1
    matchedPairs{i+1,1}(:,1:2) = RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM  + RD_L.ResultDisp{i,1}.U ;
    matchedPairs{i+1,1}(:,3:4) = RD_R.Coordinates_corr + RD_R.ResultDisp{i,1}.U;
end


