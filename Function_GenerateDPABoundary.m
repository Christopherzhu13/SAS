%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_GenerateDPABoundary.m
%%%%
%%%%        - Gather and Save the DPA Boundary Information for the given DPA
%%%%        - DPA_Boundary: for each of the N DPA Boundary Points:
%%%%            - DPA_Boundary: [Longitude, Latitude, Vertex Index, EdgeIndex, Edge Point Sequence] (Nx5 Double)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function [DPA_Boundary, DPA_Boundary_cell] = Function_GenerateDPABoundary(DPA_Zone, DPA_ID)



format long;




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
DPAESCLocation_folderName = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Data_20190531\';
DPArelated_dataSaveFolder = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Original DPA and ESC Data';
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
if (strcmp(DPA_Zone, 'east') && (DPA_ID>26)) || (strcmp(DPA_Zone, 'west') && (DPA_ID>14))
    DPAgeodata_fileName = strcat(num2str(DPA_IDstr), '.geojson');
else
    DPAgeodata_fileName = strcat(DPA_Zone, num2str(DPA_IDstr), '.geojson');
end

DPAgeodata_folderName = strcat(DPAESCLocation_folderName, '\003 DPA Geodata\DPAgeodata_20180515');
DPAgeodata_fullFileName = fullfile(DPAgeodata_folderName, DPAgeodata_fileName);

if strcmp(DPA_Zone, 'west') && DPA_ID == 13    
    
    data = loadjson(DPAgeodata_fullFileName);
    DPA_Boundary = data.features{1}.geometry.coordinates{1}(:,1:2);                 % Coordinates (Longitude, Latitude)    
    
    % Seperate the Island boundary coordinates
    DPA_Boundary_Island_1 = data.features{1}.geometry.coordinates{2}(:,1:2);        % Coordinates (Longitude, Latitude)
    DPA_Boundary_Island_2 = data.features{1}.geometry.coordinates{3}(:,1:2);        % Coordinates (Longitude, Latitude)
    DPA_Boundary_Island_3 = data.features{1}.geometry.coordinates{4}(:,1:2);        % Coordinates (Longitude, Latitude)
        
else
    data = loadjson(DPAgeodata_fullFileName);
    
    if strcmp(DPA_Zone, 'east') && ismember(DPA_ID, [30, 31])
        DPA_Boundary = data.features{1}.geometry.coordinates(1:2);                  % Coordinates (Longitude, Latitude)
    else
        DPA_Boundary = data.features{1}.geometry.coordinates{1}(:,1:2);             % Coordinates (Longitude, Latitude)
    end
end













% Load DPA Vertex Information
DPA_Analysis_SupportingInfo_fileName = 'DPA_Analysis_SupportingInfo_msohul20180815.xlsx';
DPA_Analysis_SupportingInfo_fullFileName = fullfile(DPAESCLocation_folderName, DPA_Analysis_SupportingInfo_fileName);
DPAVertexIndex_temp = readtable(DPA_Analysis_SupportingInfo_fullFileName, 'sheet', 'Sheet1');
DPAInfoIndex = find(strcmp(DPAVertexIndex_temp.DPAID, DPA_Name) == 1);
DPAVertexIndex = [DPAVertexIndex_temp.DPAVertexIndex1(DPAInfoIndex), DPAVertexIndex_temp.DPAVertexIndex2(DPAInfoIndex), ...
    DPAVertexIndex_temp.DPAVertexIndex3(DPAInfoIndex), DPAVertexIndex_temp.DPAVertexIndex4(DPAInfoIndex)];









% Get the Edge Point indices for each of the edges of the given DPA
EndPointEdge = DPAVertexIndex_temp.EndPointEdge(DPAInfoIndex);
for edgeIndex = 1 : 4
    
    EdgePointIndex{1} = DPAVertexIndex(1) : DPAVertexIndex(2)-1;
    EdgePointIndex{2} = DPAVertexIndex(2) : DPAVertexIndex(3)-1;
    EdgePointIndex{3} = DPAVertexIndex(3) : DPAVertexIndex(4)-1;
    EdgePointIndex{4} = DPAVertexIndex(4) : DPAVertexIndex(1)-1;
    
    if EndPointEdge == 1
        EdgePointIndex{1} = [DPAVertexIndex(1) : size(DPA_Boundary, 1), 1 : DPAVertexIndex(2)-1];
        
    elseif EndPointEdge == 2
        EdgePointIndex{2} = [DPAVertexIndex(2) : size(DPA_Boundary, 1), 1 : DPAVertexIndex(3)-1];
        
    elseif EndPointEdge == 3
        EdgePointIndex{3} = [DPAVertexIndex(3) : size(DPA_Boundary, 1), 1 : DPAVertexIndex(4)-1];
        
    elseif EndPointEdge == 4
        EdgePointIndex{4} = [DPAVertexIndex(4) : size(DPA_Boundary, 1), 1 : DPAVertexIndex(1)-1];
        
    end
    
end






% Add Vertex Index, Edge Index, and Edge Point Sequence to the corresponding DPA Boundary Points
DPA_Boundary(:, 3) = 0;
for edgeIndex = 1 : 4    
    
    % Add the Vertex Index to cooresponding DPA Boundary points
    DPA_Boundary(DPAVertexIndex(1, edgeIndex), 3) = edgeIndex;
    
    % Add Edge Index to cooresponding DPA Boundary points
    DPA_Boundary(EdgePointIndex{1, edgeIndex}, 4) = edgeIndex;    
    
    % Add the Edge Point Sequence to cooresponding DPA Boundary points
    DPA_Boundary(EdgePointIndex{1, edgeIndex}, 5) = ( 1 : size(EdgePointIndex{1, edgeIndex},2) )';
    
end









% Save DPA Boundary Information
DPA_Boundary_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_Boundary.mat');
DPA_Boundary_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_Boundary_fileName);
if exist(DPA_Boundary_fullFileName, 'file') == 2
    delete(DPA_Boundary_fullFileName)
end
clear DPA_Boundary_fileName

if strcmp(DPA_Zone, 'west') && DPA_ID == 13
    DPA_Boundary_cell{1} = DPA_Boundary;
    DPA_Boundary_cell{2} = DPA_Boundary_Island_1;
    DPA_Boundary_cell{3} = DPA_Boundary_Island_2;
    DPA_Boundary_cell{4} = DPA_Boundary_Island_3;
    
    save(DPA_Boundary_fullFileName, 'DPA_Boundary_cell')
    clear DPA_Boundary_fullFileName
    
else
    DPA_Boundary_cell = [];
    save(DPA_Boundary_fullFileName, 'DPA_Boundary')
    clear DPA_Boundary_fullFileName
    
end














xx=1;       fprintf('\nDPA Analysis: Generation of DPA Boundary Information :: Complete\n')













end













