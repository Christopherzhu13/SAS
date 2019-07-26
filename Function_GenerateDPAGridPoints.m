%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_GenerateDPAGridPoints.m
%%%%
%%%%        - Function that Generates and Saves DPA Grid Point Info for the given DPA
%%%%            - DPA Grid Point Location Information
%%%%            - DPA Grid Point Inner Boundary Index Information
%%%%
%%%%        - Generates and saves DPA Grid Point Information Table
%%%%            - DPAName_DPAGridPointInfo.xlsx
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function DPA_gridInfo = Function_GenerateDPAGridPoints(DPA_Zone, DPA_ID)



format long;
addpath('C:\Users\MunawwarSohul\Dropbox (Federated Wireless)\007 Simulations\001 ITM CodeFolder')





if strcmp(DPA_Zone, 'east') && DPA_ID == 27
    DPA_IDstr = 'Norfolk';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 28
    DPA_IDstr = 'Mayport';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 29
    DPA_IDstr = 'Pascagoula';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 30
    DPA_IDstr = 'Pensacola';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 31
    DPA_IDstr = 'WebsterField';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 15
    DPA_IDstr = 'SanDiego';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 16
    DPA_IDstr = 'Alameda';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 17
    DPA_IDstr = 'Long-Beach';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 18
    DPA_IDstr = 'Bremerton-Everett';
else
    DPA_IDstr = DPA_ID;
end     % If ends after assigning numeric values to DPA with string names






DPA_Name = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
DPArelated_dataSaveFolder = 'C:\Users\MunawwarSohul\Dropbox (Federated Wireless)\007 Simulations\002 ESC Logic System\Original DPA and ESC Data';
DPAspecific_dataSaveFolder = strcat(DPArelated_dataSaveFolder, '\', DPA_Name);
if exist(DPAspecific_dataSaveFolder, 'dir') ~= 7
    mkdir(DPAspecific_dataSaveFolder)
end













xx=1;   % Completes Input initialization














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Gather and Save the DPA Boundary Information for the given DPA
%%%%    DPA_Boundary: for each of the N DPA Boundary Points:
%%%%            - DPA_Boundary: [Longitude, Latitude, Vertex Index, EdgeIndex, Edge Point Sequence] (Nx5 Double)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(DPA_Zone, 'west') && DPA_ID == 13
    [~, DPA_Boundary_cell] = Function_GenerateDPABoundary(DPA_Zone, DPA_ID);
    DPA_Boundary = DPA_Boundary_cell{1};
    DPA_Boundary_Island_1 = DPA_Boundary_cell{2};
    DPA_Boundary_Island_2 = DPA_Boundary_cell{3};
    DPA_Boundary_Island_3 = DPA_Boundary_cell{4};
    
else
    DPA_Boundary = Function_GenerateDPABoundary(DPA_Zone, DPA_ID);
    
end

% Find the Edge points for each of the 4 edges of the DPA
DPA_Edge1_temp1 = DPA_Boundary(DPA_Boundary(:,4)==1, :);
[~, DPA_Edge1_sortIndex] = sort(DPA_Edge1_temp1(:,5));
DPA_Edge1points = DPA_Edge1_temp1(DPA_Edge1_sortIndex,:);


DPA_Edge2_temp1 = DPA_Boundary(DPA_Boundary(:,4)==2, :);
[~, DPA_Edge2_sortIndex] = sort(DPA_Edge2_temp1(:,5));
DPA_Edge2points = DPA_Edge2_temp1(DPA_Edge2_sortIndex,:);


DPA_Edge3_temp1 = DPA_Boundary(DPA_Boundary(:,4)==3, :);
[~, DPA_Edge3_sortIndex] = sort(DPA_Edge3_temp1(:,5));
DPA_Edge3points = DPA_Edge3_temp1(DPA_Edge3_sortIndex,:);


DPA_Edge4_temp1 = DPA_Boundary(DPA_Boundary(:,4)==4, :);
[~, DPA_Edge4_sortIndex] = sort(DPA_Edge4_temp1(:,5));
DPA_Edge4points = DPA_Edge4_temp1(DPA_Edge4_sortIndex,:);









xx=1;














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Gather and Save the DPA Inner Boundary Information for the given DPA
%%%%    DPA_Boundary: for each of the N DPA Boundary Points:
%%%%            - DPA_Boundary: [Longitude, Latitude, Vertex Index, EdgeIndex, Edge Point Sequence] (Nx5 Double)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DPA_innerBoundary = Function_GenerateDPAInnerBoundary(DPA_Zone, DPA_ID);










xx=1;














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Generate the DPA Grid Points for the given Grid Point Resolution
%%%%    DPA_gridInfo: for each of the N DPA Grid Points:
%%%%            - DPA_gridInfo: [Longitude, Latitude, 75km Index] (Nx3 Double)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GridPointResolution = 1;
if strcmp(DPA_Zone, 'east') && (strcmp(DPA_IDstr, 'Norfolk') || strcmp(DPA_IDstr, 'Pascagoula')) 
    dGridGranularity = 500;
elseif strcmp(DPA_Zone, 'east') && strcmp(DPA_IDstr, 'Mayport')
    dGridGranularity = 100;
elseif strcmp(DPA_Zone, 'west') && strcmp(DPA_IDstr, 'SanDiego')
    dGridGranularity = 500;
elseif strcmp(DPA_Zone, 'west') && strcmp(DPA_IDstr, 'Alameda')
    dGridGranularity = 100;
elseif strcmp(DPA_Zone, 'west') && strcmp(DPA_IDstr, 'Bremerton-Everett')
    dGridGranularity = 1000;
else
    dGridGranularity = GridPointResolution * 1000;
end     % Setting the grid resolution










% Generating the grid points covering the rectangle with vertices as below
dNorthLatitude = max(DPA_Boundary(:, 2));         dWestLongitude = min(DPA_Boundary(:, 1));
dEastLongitude = max(DPA_Boundary(:, 1));         dSouthLatitude = min(DPA_Boundary(:, 2));

tCellularGrid=ComputeCellularGridPoints(dNorthLatitude,dSouthLatitude,dEastLongitude,dWestLongitude,dGridGranularity);
gridLat = tCellularGrid.p2dGridPointLat;
gridLon = tCellularGrid.p2dGridPointLong;





% Selecting the grid points that fall within the DPA Boundary
if strcmp(DPA_Zone, 'west') && DPA_ID == 13
    
    % Taking out any grid points that fall on any of the Islands
    insideLocationIndex = inpolygon(gridLon,gridLat,DPA_Boundary(:,1),DPA_Boundary(:,2));
    insideIsland2LocationIndex = inpolygon(gridLon,gridLat,DPA_Boundary_Island_1(:,1),DPA_Boundary_Island_1(:,2));
    insideIsland3LocationIndex = inpolygon(gridLon,gridLat,DPA_Boundary_Island_2(:,1),DPA_Boundary_Island_2(:,2));
    insideIsland4LocationIndex = inpolygon(gridLon,gridLat,DPA_Boundary_Island_3(:,1),DPA_Boundary_Island_3(:,2));
    
    insidegridIndex = insideLocationIndex & ~insideIsland2LocationIndex & ~insideIsland3LocationIndex ...
        & ~insideIsland4LocationIndex;
    
    DPA_gridLat = gridLat(insidegridIndex);
    DPA_gridLon = gridLon(insidegridIndex);
    
    
else
    
    % Find the incides of the Grid Points that are located inside the DPA
    insideLocationIndex = inpolygon(gridLon,gridLat,DPA_Boundary(:, 1),DPA_Boundary(:, 2));
    
    DPA_gridLat = gridLat(insideLocationIndex);
    DPA_gridLon = gridLon(insideLocationIndex);
    
end


DPA_gridInfo = [DPA_gridLon, DPA_gridLat];          % Coordinates (Longitude, Latitude)










% Generate the DPA Part-1 Polygon and Index the Grid Points inside the Part-1 Polygon
if (strcmp(DPA_Zone, 'east') && (DPA_ID>26)) || (strcmp(DPA_Zone, 'west') && (DPA_ID>14))
    
    DPA_gridInfo(:, 3) = ones(size(DPA_gridInfo, 1), 1);
    
    
    
else
    
    % Selecting appropriate Part-1 Polygon Generation Function for the selected DPA
    group1E = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 25, 26];
    group1W = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    group2E = 17;	group3E = 18;
    
    
    
    
    if (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, group1E)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, group1W))
        
        
        % Find the DPA Edge-1 Points for the DPA Part-1 Polygon        
        DPA_Part1polygon_Edge1 = DPA_Edge1points(:,1:2);
        
        
        % Find the DPA Edge-2 Points for the DPA Part-1 Polygon
        innerBoundaryEndPoint_to_Edge2PP_Distance = sqrt( (DPA_innerBoundary(end,1) - DPA_Edge2points(:, 1)).^2 + ...
            (DPA_innerBoundary(end,2) - DPA_Edge2points(:, 2)).^2 );
        [~, Edge2_PPstart] = min(innerBoundaryEndPoint_to_Edge2PP_Distance);
        DPA_Part1polygon_Edge2 = DPA_Edge2points(1:Edge2_PPstart, 1:2);
        
        % Find the Inner Boundary Points for the DPA Part-1 Polygon
        DPA_Part1polygon_InnerBoundary = flipud(DPA_innerBoundary(:,1:2));
        
        
        % Find the DPA Edge-4 Points for the DPA Part-1 Polygon
        innerBoundaryStartPoint_to_Edge4PP_Distance = sqrt( (DPA_innerBoundary(1,1) - DPA_Edge4points(:, 1)).^2 + ...
            (DPA_innerBoundary(1,2) - DPA_Edge4points(:, 2)).^2 );
        [~, Edge4_PPstart] = min(innerBoundaryStartPoint_to_Edge4PP_Distance);
        DPA_Part1polygon_Edge4 = DPA_Edge4points(Edge4_PPstart:end, 1:2);
        
        
        
        % Generate the Part-1 Polygon
        DPA_Part1polygon = [DPA_Part1polygon_Edge1; DPA_Part1polygon_Edge2; DPA_Part1polygon_InnerBoundary; DPA_Part1polygon_Edge4; DPA_Part1polygon_Edge1(1, 1:2)];
        
        
        
        
        
        
        
        
        
        
    elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, group2E))
        
        
        % Find the DPA Edge-1 Points for the DPA Part-1 Polygon
        innerBoundaryStartPoint_to_Edge1PP_Distance = sqrt( (DPA_innerBoundary(end,1) - DPA_Edge1points(:, 1)).^2 + ...
            (DPA_innerBoundary(end,2) - DPA_Edge1points(:, 2)).^2 );
        [~, Edge1_PPstart] = min(innerBoundaryStartPoint_to_Edge1PP_Distance);
        DPA_Part1polygon_Edge1 = DPA_Edge1points(1:Edge1_PPstart, 1:2);
        
        
        % Find the Inner Boundary Points for the DPA Part-1 Polygon
        DPA_Part1polygon_InnerBoundary = flipud(DPA_innerBoundary(:,1:2));
        
        
        % Find the DPA Edge-4 Points for the DPA Part-1 Polygon
        innerBoundaryStartPoint_to_Edge4PP_Distance = sqrt( (DPA_innerBoundary(1,1) - DPA_Edge4points(:, 1)).^2 + ...
            (DPA_innerBoundary(1,2) - DPA_Edge4points(:, 2)).^2 );
        [~, Edge4_PPstart] = min(innerBoundaryStartPoint_to_Edge4PP_Distance);
        DPA_Part1polygon_Edge4 = DPA_Edge4points(Edge4_PPstart:end, 1:2);
        
        
        
        % Generate the Part-1 Polygon
        DPA_Part1polygon = [DPA_Part1polygon_Edge1; DPA_Part1polygon_InnerBoundary; DPA_Part1polygon_Edge4; DPA_Part1polygon_Edge1(1, 1:2)];
        
        
        
        
        
        
        
        
        
        
    elseif strcmp(DPA_Zone, 'east') && ismember(DPA_ID, group3E)        
        
        
        % Find the DPA Edge-1 Points for the DPA Part-1 Polygon
        innerBoundaryEndPoint_to_Edge1PP_Distance = sqrt( (DPA_innerBoundary(1,1) - DPA_Edge1points(:, 1)).^2 + ...
            (DPA_innerBoundary(1,2) - DPA_Edge1points(:, 2)).^2 );
        [~, Edge1_PPstart] = min(innerBoundaryEndPoint_to_Edge1PP_Distance);
        DPA_Part1polygon_Edge1 = DPA_Edge1points(Edge1_PPstart:end, 1:2);
        
        
        
        % Find the DPA Edge-2 Points for the DPA Part-1 Polygon
        innerBoundaryEndPoint_to_Edge2PP_Distance = sqrt( (DPA_innerBoundary(end,1) - DPA_Edge2points(:, 1)).^2 + ...
            (DPA_innerBoundary(end,2) - DPA_Edge2points(:, 2)).^2 );
        [~, Edge2_PPstart] = min(innerBoundaryEndPoint_to_Edge2PP_Distance);
        DPA_Part1polygon_Edge2 = DPA_Edge2points(1:Edge2_PPstart, 1:2);
        
        % Find the Inner Boundary Points for the DPA Part-1 Polygon
        DPA_Part1polygon_InnerBoundary = flipud(DPA_innerBoundary(:,1:2));        
        
        
        % Generate the Part-1 Polygon
        DPA_Part1polygon = [DPA_Part1polygon_Edge1; DPA_Part1polygon_Edge2; DPA_Part1polygon_InnerBoundary; DPA_Part1polygon_Edge1(1, 1:2)];
        
        
        
        
    end
    
    
    
    
    
    
    % Index the DPA Grid Points ibelonging to Part-1  and Part-2 Polygon
    % Part-1 Polygon Grid Points: DPA_gridInfo(DPA_gridInfo(:,3)==1, :)
    % Part-2 Polygon Grid Points: DPA_gridInfo(DPA_gridInfo(:,3)==0, :)
    GPinsideDPAPart1polygonIndex = inpolygon(DPA_gridInfo(:,1), DPA_gridInfo(:,2), DPA_Part1polygon(:,1), DPA_Part1polygon(:,2));
    GPinsideDPAPart1polygonIndex2 = find(GPinsideDPAPart1polygonIndex==1);
    
    DPA_gridInfo(:, 3) = zeros(size(DPA_gridInfo, 1), 1);
    DPA_gridInfo(GPinsideDPAPart1polygonIndex2, 3) = ones(size(GPinsideDPAPart1polygonIndex2, 1), 1);
    
    
    
end










% Save DPA_gridInfo
DPA_gridInfo_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_gridInfo.mat');
DPA_gridInfo_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_gridInfo_fileName);
if exist(DPA_gridInfo_fullFileName, 'file') == 2
    delete(DPA_gridInfo_fullFileName)
end
save(DPA_gridInfo_fullFileName, 'DPA_gridInfo')
clear DPA_gridInfo_fileName DPA_gridInfo_fullFileName














xx=1;       fprintf('\nESC Logic System Analysis: Generation of DPA Grid Point Information :: Complete\n')











end






