%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_GenerateDPAInnerBoundary.m
%%%%
%%%%    - Gather and Save the 75km Inner Boundary Information for the given DPA
%%%%    - DPA_innerBoundary: for each of the N Inner Boundary Points:
%%%%            - DPA_innerBoundary: [Longitude, Latitude, Part Index] (Nx3 Double)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function DPA_innerBoundary = Function_GenerateDPAInnerBoundary(DPA_Zone, DPA_ID)



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
DPAESCLocation_folderName = 'C:\Users\MunawwarSohul\Dropbox (Federated Wireless)\007 Simulations\003 Related FW Documents\Data for Analysis';

DPArelated_dataSaveFolder = 'C:\Users\MunawwarSohul\Dropbox (Federated Wireless)\007 Simulations\002 ESC Logic System\Original DPA and ESC Data';
DPAspecific_dataSaveFolder = strcat(DPArelated_dataSaveFolder, '\', DPA_Name);
if exist(DPAspecific_dataSaveFolder, 'dir') ~= 7
    mkdir(DPAspecific_dataSaveFolder)
end





% Load DPA Boundary Information
DPA_Boundary_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_Boundary.mat');
DPA_Boundary_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_Boundary_fileName);
load(DPA_Boundary_fullFileName)
clear DPA_Boundary_fullFileName DPA_Boundary_fullFileName

if strcmp(DPA_Zone, 'west') && DPA_ID == 13
    DPA_Boundary = DPA_Boundary_cell{1};
    DPA_Boundary_Island_1 = DPA_Boundary_cell{2};
    DPA_Boundary_Island_2 = DPA_Boundary_cell{3};
    DPA_Boundary_Island_3 = DPA_Boundary_cell{4};
        
end














xx=1;   % Completes Input initialization














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Gather and Save the 75km Inner Boundary Information for the given DPA
%%%%    DPA_innerBoundary: for each of the N Inner Boundary Points:
%%%%            - DPA_InnerBoundary: [Longitude, Latitude, Part Index] (Nx3 Double)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(DPA_Zone, 'east')
    innerBoundary75kmLine_fileName = 'east_75km_line.geojson';
elseif strcmp(DPA_Zone, 'west')
    innerBoundary75kmLine_fileName = 'west_75km_line.geojson';
end

innerBoundary75kmLine_folderName = strcat(DPAESCLocation_folderName, '\003 DPA Geodata\DPAgeodata_20180515');
innerBoundary75kmLine_fullFileName = fullfile(innerBoundary75kmLine_folderName, innerBoundary75kmLine_fileName);

innerBoundary75kmLine_data = loadjson(innerBoundary75kmLine_fullFileName);
innerBoundary75kmLine = innerBoundary75kmLine_data.features{1}.geometry.coordinates;             % Coordinates (Longitude, Latitude)


% Find the Inner Boundary for the given DPA
DPA_innerBoundaryIndex = inpolygon(innerBoundary75kmLine(:,1), innerBoundary75kmLine(:,2), DPA_Boundary(:,1), DPA_Boundary(:,2));
DPA_innerBoundary = [innerBoundary75kmLine(DPA_innerBoundaryIndex, 1), innerBoundary75kmLine(DPA_innerBoundaryIndex, 2)];



% Indexing all parts of Inner Boundary (the parts are caused by 'inpolygon' operation)
DPA_innerBoundary(:, 3) = 0;

DPA_innerBoundary_PartSearchFlag = 0;
whileCount = 0;             startIndex = 1;
selectStartIndex = 1;       selectEndIndex = length(DPA_innerBoundaryIndex);

while DPA_innerBoundary_PartSearchFlag == 0
    
    whileCount = whileCount + 1;
    
    
    % Find the Start and End Index of the current Part of the Inner Boundary
    DPA_innerBoundary_PartStartIndex = find(DPA_innerBoundaryIndex(selectStartIndex:selectEndIndex) == 1, 1, 'first') + (selectStartIndex-1);
    selectStartIndex = DPA_innerBoundary_PartStartIndex + 1;
    DPA_innerBoundary_PartEndIndex = find(DPA_innerBoundaryIndex(selectStartIndex:selectEndIndex) == 0, 1, 'first') - 1 + (selectStartIndex-1);
    selectEndIndex = DPA_innerBoundary_PartEndIndex + 1;
    
    
    
    % Check for the end of the Boundary condition
    if isempty(DPA_innerBoundary_PartEndIndex)
        endIndex = length(DPA_innerBoundaryIndex) - DPA_innerBoundary_PartStartIndex + startIndex;
    else
        endIndex = DPA_innerBoundary_PartEndIndex - DPA_innerBoundary_PartStartIndex + startIndex;
    end
    
    
       
    
    
    % Index the Inner Boundary Points according to Part Index
    if ~isempty(selectStartIndex) && ~isempty(selectEndIndex)
        DPA_innerBoundary(startIndex : endIndex, 3) = whileCount;
        startIndex = endIndex + 1;
        selectStartIndex = selectEndIndex;
        selectEndIndex = length(DPA_innerBoundaryIndex);
        
    elseif ~isempty(selectStartIndex) && isempty(selectEndIndex)
        DPA_innerBoundary(startIndex : endIndex, 3) = whileCount;
        DPA_innerBoundary_PartSearchFlag = 1;
        
    else
        DPA_innerBoundary_PartSearchFlag = 1;
        
    end
    
    
    
end




% Save DPA inner Inner Boundary
DPA_innerBoundary_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_innerBoundary.mat');
DPA_innerBoundary_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_innerBoundary_fileName);
if exist(DPA_innerBoundary_fullFileName, 'file') == 2
    delete(DPA_innerBoundary_fullFileName)
end
save(DPA_innerBoundary_fullFileName, 'DPA_innerBoundary')
clear DPA_innerBoundary_fileName DPA_innerBoundary_fullFileName














xx=1;       fprintf('\nESC Logic System Analysis: Generation of DPA Inner Boundary Information :: Complete\n')













end













