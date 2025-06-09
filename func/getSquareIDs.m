function square_ids = getSquareIDs(M, N, i, j, R, indexMap)

    % col_ids = repmat(col_range, length(row_range), 1); 
    % row_ids = repmat(row_range', 1, length(col_range)); 
    
    flag = 1;
    square_ids = ones((2*R+1)*(2*R+1),1);
    
%try
    for row = (i-R):(i+R)
        for col = (j-R):(j+R)
            square_ids(flag) = indexMap(row,col);
            flag = flag + 1;
        end
    end
% catch
%     disp(i);
%     disp(j);
% end

end