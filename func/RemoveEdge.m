function [FinalResult,coordinatesFEM,RD_L,xList,yList,M,N] = RemoveEdge(FinalResult,coordinatesFEM,RD_L,xList,yList,M,N)
% when removing edge, FR.coor, FR.disp, should be adjusted
% RD_L.Disp and RD_L.ResultFEMeshEachFrame should be adjusted as well due to ploting strain's need

RemoveCols_L = input('How many columns are removed from the left?(0,1,2...) \nInput: ');
RemoveCols_R = input('How many columns are removed from the right?(0,1,2...) \nInput: ');
RemoveRows_L = input('How many rows are removed from the left?(0,1,2...) \nInput: ');
RemoveRows_R = input('How many rows are removed from the right?(0,1,2...) \nInput: ');
coordinatesFEM_x = reshape(coordinatesFEM(:,1),[M,N]);
coordinatesFEM_y = reshape(coordinatesFEM(:,2),[M,N]);


for i = 1:length(FinalResult.Displacement)
    u_temp = reshape(FinalResult.Displacement{i,1},[M,N]);
    v_temp = reshape(FinalResult.Displacement{i,2},[M,N]);
    w_temp = reshape(FinalResult.Displacement{i,3},[M,N]);
    x_temp = reshape(FinalResult.Coordinates{i,1},[M,N]);
    y_temp = reshape(FinalResult.Coordinates{i,2},[M,N]);
    z_temp = reshape(FinalResult.Coordinates{i,3},[M,N]);
    if i < length(FinalResult.Displacement) % Since RD_L_disp has one less cell
        L_u_temp = reshape(RD_L.ResultDisp{i}.U(1:2:end),[M,N]);
        L_v_temp = reshape(RD_L.ResultDisp{i}.U(2:2:end),[M,N]);
        L_x_temp = reshape(RD_L.ResultFEMeshEachFrame{i}.coordinatesFEM(:,1),[M,N]);
        L_y_temp = reshape(RD_L.ResultFEMeshEachFrame{i}.coordinatesFEM(:,2),[M,N]);
    end

    if mod(RemoveCols_L,1)==0 && RemoveCols_L>0
        u_temp(1:RemoveCols_L,:) = [];
        v_temp(1:RemoveCols_L,:) = [];
        w_temp(1:RemoveCols_L,:) = [];
        x_temp(1:RemoveCols_L,:) = [];
        y_temp(1:RemoveCols_L,:) = [];
        z_temp(1:RemoveCols_L,:) = [];
        if i < length(FinalResult.Displacement) % Since RD_L_disp has one less cell
            L_u_temp(1:RemoveCols_L,:) = [];
            L_v_temp(1:RemoveCols_L,:) = [];
            L_x_temp(1:RemoveCols_L,:) = [];
            L_y_temp(1:RemoveCols_L,:) = [];
        end
        if i == 1 % only do 1 time
            xList(1:RemoveCols_L) = [];
            coordinatesFEM_x(1:RemoveCols_L,:) = [];
            coordinatesFEM_y(1:RemoveCols_L,:) = [];
        end
    end

    if mod(RemoveCols_R,1)==0 && RemoveCols_R>0
        u_temp(end-RemoveCols_R+1:end,:) = [];
        v_temp(end-RemoveCols_R+1:end,:) = [];
        w_temp(end-RemoveCols_R+1:end,:) = [];
        x_temp(end-RemoveCols_R+1:end,:) = [];
        y_temp(end-RemoveCols_R+1:end,:) = [];
        z_temp(end-RemoveCols_R+1:end,:) = [];
        if i < length(FinalResult.Displacement) % Since RD_L_disp has one less cell
            L_u_temp(end-RemoveCols_R+1:end,:) = [];
            L_v_temp(end-RemoveCols_R+1:end,:) = [];
            L_x_temp(end-RemoveCols_R+1:end,:) = [];
            L_y_temp(end-RemoveCols_R+1:end,:) = [];
        end
        if i == 1 % only do 1 time
            xList(end-RemoveCols_R+1:end) = [];
            coordinatesFEM_x(end-RemoveCols_R+1:end,:) = [];
            coordinatesFEM_y(end-RemoveCols_R+1:end,:) = [];
        end
    end

    if mod(RemoveRows_L,1)==0 && RemoveRows_L>0
        u_temp(:,1:RemoveRows_L) = [];
        v_temp(:,1:RemoveRows_L) = [];
        w_temp(:,1:RemoveRows_L) = [];
        x_temp(:,1:RemoveRows_L) = [];
        y_temp(:,1:RemoveRows_L) = [];
        z_temp(:,1:RemoveRows_L) = [];
        if i < length(FinalResult.Displacement) % Since RD_L_disp has one less cell
            L_u_temp(:,1:RemoveRows_L) = [];
            L_v_temp(:,1:RemoveRows_L) = [];
            L_x_temp(:,1:RemoveRows_L) = [];
            L_y_temp(:,1:RemoveRows_L) = [];
        end
        if i == 1 % only do 1 time
            yList(1:RemoveRows_L) = [];
            coordinatesFEM_x(1:RemoveRows_L,:) = [];
            coordinatesFEM_y(1:RemoveRows_L,:) = [];
        end
    end

    if mod(RemoveRows_R,1)==0 && RemoveRows_R>0
        u_temp(:,end-RemoveRows_R+1:end) = [];
        v_temp(:,end-RemoveRows_R+1:end) = [];
        w_temp(:,end-RemoveRows_R+1:end) = [];
        x_temp(:,end-RemoveRows_R+1:end) = [];
        y_temp(:,end-RemoveRows_R+1:end) = [];
        z_temp(:,end-RemoveRows_R+1:end) = [];
        if i < length(FinalResult.Displacement) % Since RD_L_disp has one less cell
            L_u_temp(:,end-RemoveRows_R+1:end) = [];
            L_v_temp(:,end-RemoveRows_R+1:end) = [];
            L_x_temp(:,end-RemoveRows_R+1:end) = [];
            L_y_temp(:,end-RemoveRows_R+1:end) = [];
        end
        if i == 1 % only do 1 time
            yList(end-RemoveRows_R+1:end) = [];
            coordinatesFEM_x(:,end-RemoveRows_R+1:end) = [];
            coordinatesFEM_y(:,end-RemoveRows_R+1:end) = [];
        end
    end

    FinalResult.Displacement{i,1} = reshape(u_temp,[length(xList)*length(yList),1]);
    FinalResult.Displacement{i,2} = reshape(v_temp,[length(xList)*length(yList),1]);
    FinalResult.Displacement{i,3} = reshape(w_temp,[length(xList)*length(yList),1]);
    FinalResult.Coordinates{i,1} = reshape(x_temp,[length(xList)*length(yList),1]);
    FinalResult.Coordinates{i,2} = reshape(y_temp,[length(xList)*length(yList),1]);
    FinalResult.Coordinates{i,3} = reshape(z_temp,[length(xList)*length(yList),1]);
    if i < length(FinalResult.Displacement) % Since RD_L_disp has one less cell
        temp(1,:) = reshape(L_u_temp,[1,length(xList)*length(yList)]);
        temp(2,:) = reshape(L_v_temp,[1,length(xList)*length(yList)]);
        RD_L.ResultDisp{i}.U = temp{i}(:);
        temp2(:,1) = reshape(L_x_temp,[length(xList)*length(yList),1]);
        temp2(:,2) = reshape(L_y_temp,[length(xList)*length(yList),1]);
        RD_L.ResultFEMeshEachFrame{i}.coordinatesFEM = temp2;
    end
end

M = length(xList);
N = length(yList);
coordinatesFEM_removed(:,1) = reshape(coordinatesFEM_x,[M*N,1]);
coordinatesFEM_removed(:,2) = reshape(coordinatesFEM_y,[M*N,1]);
coordinatesFEM = coordinatesFEM_removed;

end