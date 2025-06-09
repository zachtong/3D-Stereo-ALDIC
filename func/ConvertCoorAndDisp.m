function  FinalResult = ConvertCoorAndDisp(FinalResult,RotationMatrix,TranslationMatrix)
    for i = 1:size(FinalResult.Coordinates,1)
    coortemp = FinalResult.Coordinates(i,:);
    coortemp = [coortemp{1},coortemp{2},coortemp{3}]';
    coorNew = RotationMatrix*(coortemp-TranslationMatrix');
    coorNew = coorNew';
    FinalResult.CoordinatesNew{i,1} = coorNew(:,1);
    FinalResult.CoordinatesNew{i,2} = coorNew(:,2);
    FinalResult.CoordinatesNew{i,3} = coorNew(:,3);
    end
    
    try 
    FinalResult.DisplacementNew(1,:) = FinalResult.Displacement_smooth(1,:);
    catch
    FinalResult.DisplacementNew(1,:) = FinalResult.Displacement(1,:);
    end


    for i = 2:size(FinalResult.Coordinates,1)
    FinalResult.DisplacementNew{i,1} = FinalResult.CoordinatesNew{i,1} - FinalResult.CoordinatesNew{1,1};
    FinalResult.DisplacementNew{i,2} = FinalResult.CoordinatesNew{i,2} - FinalResult.CoordinatesNew{1,2};
    FinalResult.DisplacementNew{i,3} = FinalResult.CoordinatesNew{i,3} - FinalResult.CoordinatesNew{1,3};
    end
 

end
