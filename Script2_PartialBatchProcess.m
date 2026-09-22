%% Partial Batch Process Script
% If you have not run the "SingleImage" script to generate a parameters .m
% file, then this script will not work.
% This script presumes that you have already identified microglia using
% Script1. 

%% Setup
clear; clc; %clear memory and command window log
fullScriptPath = mfilename('fullpath');
ScriptFolder = fileparts(fullScriptPath);
addpath(genpath(ScriptFolder));
fprintf('Added %s and all its subfolders to the MATLAB path.\n', ScriptFolder);

app = BatchFileDataPartial_App;
waitfor(app, 'Output');
delete(app);
addpath(genpath(ParametersPath));
addpath(genpath(BatchFolder));
addpath(genpath(SaveBatchPath));

% collect command window outputs in log file for later reference if needed
logFile = fullfile(SaveBatchPath, 'memory_log.txt');
diary(logFile);
diary on;

% Get list of images
if SubFolders == 0
    imgFiles = dir(fullfile(BatchFolder, '*.tif'));
elseif SubFolders == 1
    imgFiles = dir(fullfile(BatchFolder, '**', '*.tif'));
end

fprintf('Found %d images to process.\n', numel(imgFiles));

%% Get list of FullMg 
FullMgFiles = dir(fullfile(FullMgFolder, '*.mat'));
fprintf('Found %d mat files for microglia.\n', numel(FullMgFiles));

fprintf('Batch process begins %s. \n', datetime('now'));
fprintf('----------------------------------------------------------------------------\n')

%% Batch Loop
tic; % timer
for k = 1:numel(imgFiles)

    % Reload parameters
    load(ParametersFile);

    filename = imgFiles(k).name;
    [~, name, ~] = fileparts(filename);
    filepath = fullfile(imgFiles(k).folder, imgFiles(k).name);

    subpath = erase(imgFiles(k).folder, BatchFolder);
    subpath = char(subpath);
    folderName = [name, '_data'];
    SaveFolderPath = fullfile(SaveBatchPath, subpath, folderName);
    if ~exist(SaveFolderPath, 'dir')
        mkdir(SaveFolderPath);
    end

    % Find and load matching FullMg .mat file
    matName = ['FullMg_', name, '.mat'];
    matIdx = find(strcmp({FullMgFiles.name}, matName));

    if isempty(matIdx)
        warning('No matching FullMg file found for %s. Skipping.', filename);
        continue;
    end

    matFilePath = fullfile(FullMgFiles(matIdx).folder, FullMgFiles(matIdx).name);
    load(matFilePath);
    numObjMg = numel(FullMg);

    fprintf('Processing image file: %s...\n', filename);


    %% IMPORT: LOAD TIFF FILE
    reader = bfGetReader(filepath);
    reader.setSeries(0);
    sizeX = reader.getSizeX(); % pixels in X
    sizeY = reader.getSizeY(); % pixels in Y
    sizeZ = reader.getSizeZ(); % number of Z-slices
    sizeC = reader.getSizeC(); % number of channels
    sizeXY = [sizeY, sizeX];
    ImageVol  = (sizeXY(1)*sizeXY(2)*sizeZ) * voxscale;
    ImageDataTable = table();
    ImageDataTable.ImageVol = ImageVol(:);
            
    img = zeros(sizeY, sizeX, sizeZ, sizeC, 'uint16');
    for c = 1:sizeC
        for z = 1:sizeZ
            index = reader.getIndex(z - 1, c - 1, 0);
            plane = bfGetPlane(reader, index + 1);
            img(:, :, z, c) = plane;
        end
    end

    sizeC = size(img, 4);
    imgChannels = cell(1, sizeC);
    for c = 1:sizeC
        imgChannels{c} = squeeze(img(:, :, :, c));
    end
    
    imgGlia = imgChannels{round(GliaChannel)};
    sz = size(imgGlia);
    clear imgGlia

    CellNum = (1:numel(FullMg))';
    CellDataTable = table(CellNum, ... 
        'VariableNames', {'CellNum'});

    %% COLOCALIZATION

    if (Q.ColocAns == 1)
        fprintf('Processing channels for colocalization \n')
        ColocImgs = struct();
        for i = 1:numel(ColocChannels)
            chNum = ColocChannels(i);
            chName = matlab.lang.makeValidName(ColocChannelNames{i});
            ColocImgs.(chName) = imgChannels{chNum};
        end
        
        % Compute mean fluorescence per channel
        for i = 1:numel(ColocChannels)
            chNum = ColocChannels(i);
            chName = matlab.lang.makeValidName(ColocChannelNames{i});
            MeanFIPerChannel.(chName) = mean(imgChannels{chNum}, 'all');
        end
        
        % Add mean values to image data table
        channelFields = fieldnames(MeanFIPerChannel);
        prefixedFields = cellfun(@(f) ['MeanFI_', f], channelFields, 'UniformOutput', false);
        meanValues = cellfun(@(f) MeanFIPerChannel.(f), channelFields)';
        MeanFITable = array2table(meanValues, 'VariableNames', prefixedFields);
        ImageDataTable = [ImageDataTable, MeanFITable];

        % normalize
        fields = fieldnames(ColocImgs);
        ColocImgsNorm = struct();
        
        for f = 1:numel(fields)
            fieldName = fields{f};
            ColocImgsNorm.(fieldName) = mat2gray(ColocImgs.(fieldName));
        end
    end
    clear ColocImgs

    % Coloc Threshold
    if (Q.ColocAns == 1)
        fprintf('Thresholding channels \n')
        ColocSizeFilteredImgs = struct();
        chans = fieldnames(ColocImgsNorm);
        for j = 1:numel(chans)
            chName = chans{j};
            chanImg = ColocImgsNorm.(chName);

            % threshold using specified algorithm
            % same workflow as glia thresh but looped over other channels
            switch ColocThreshMethod.(chName)
                case 'Otsu'
                    level = graythresh(chanImg(:));
                case 'Multi-Otsu'
                    level = multithresh(chanImg(:));
                case 'Median'
                    level = median(chanImg(:));
                case 'Mean'
                    level = mean(chanImg(:));
            end

            if ColocThreshMode.(chName) == 1
                % threshold slices separately
                BinaryImg = zeros(sizeXY(1), sizeXY(2), sizeZ);
                for i = 1:sizeZ
                    BinaryImg(:, :, i) = imbinarize(chanImg(:, :, i), level * ColocThresholdValues.(chName));
                end
            else
                % threshold whole stack together
                BinaryImg = imbinarize(chanImg, level * ColocThresholdValues.(chName));
            end

            % Fill holes if desired
            if ColocFillHoles.(chName) == 1
                BinaryImg = bwmorph3(BinaryImg, "fill");
            end
            
            clear chanImg;

            CC = bwconncomp(BinaryImg, 26); % Find connected components in 3D
            clear BinaryImg;
            stats = regionprops3(CC, 'Volume'); % Measure voxel count for each object
            volumes = stats.Volume * voxscale; % Convert to physical volume
            keepIdx = find(volumes >= ColocSizeCutoffs.(chName).Min & volumes <= ColocSizeCutoffs.(chName).Max);
            L = labelmatrix(CC);
            ColocSizeFilteredImgs.(chName) = ismember(L, keepIdx); 
            clear L keepIdx;
        end
    end

    % Coloc Save Orig & Thresh Images    
    if Q.ColocAns == 1
        if  ColocSaveOrig == 1
            fields = fieldnames(ColocImgsNorm);
            for i = 1:numel(fields)
                chName = fields{i};
                img3D = ColocImgsNorm.(chName);
                img2D = max(img3D, [], 3);
                fig = figure('Visible', 'off');
                imagesc(img2D);
                colormap gray;
                daspect([1 1 1]);
                filename = fullfile(SaveFolderPath, [name, '_step3_', chName, '_OrigMAX.png']);
                exportgraphics(fig, filename, 'Resolution', 300);
                clear img3D img2D;
            end
            clear ColocImgsNorm;
        end
        
        if  ColocShowImg == 1    
            fields = fieldnames(ColocSizeFilteredImgs);
            for i = 1:numel(fields)
                chName = fields{i};
                img3D = ColocSizeFilteredImgs.(chName);
                depthMap = zeros(size(img3D), 'like', double(img3D));
                for z = 1:size(img3D, 3)
                    depthMap(:,:,z) = z * double(img3D(:,:,z));
                end
                img2D = max(depthMap, [], 3);
                fig = figure('Visible', 'off');
                imagesc(img2D);
                colormap hot;
                daspect([1 1 1]);
                filename = fullfile(SaveFolderPath, [name, '_step4_', chName, '_Thresh.png']);
                exportgraphics(fig, filename, 'Resolution', 300);
                clear img3D img2D depthMap;
            end
        end
    end

    if Q.ColocAns == 1
        
        if exist('BatchColocExclude', 'var') && isstruct(BatchColocExclude)
            fprintf('Applying exclusion filtering \n')

            %Retrieve masks for exclusion filtering
            fields1 = fieldnames(BatchColocExclude); %get list of exclusion filtered object names (previously defined by user in app)
            
            for i = 1:numel(fields1) % loop over list of named exclusion filtered objects
                f1 = fields1{i};
                fields2 = fieldnames(BatchColocExclude.(f1)); %grab 'ExcludeChan' 'IfChanExcl' and 'ExclMode' substructure labels
                
                for j = 1:numel(fields2) % loop over substructure labels
                    f2 = fields2{j};
                    chanName = BatchColocExclude.(f1).(f2); % grab the value of each substructure (example: 'ExcludeChan' = 'C1', 'IfChanExcl' = 'C3', 'ExclMode' = 'Voxels')
                    
                    if isfield(ColocSizeFilteredImgs, chanName) % check that the designated filtering channel images exist in 'ColocSizeFilteredImgs' struct
                        BatchColocExclude.(f1).(f2) = ColocSizeFilteredImgs.(chanName); % copy the matching images over to 'BatchColocExclude'
                    else
                        warning('Channel "%s" not found in ColocSizeFilteredImgs.', chanName);
                    end
                end
            end 
            
            % Filter by exclusion criteria
            if ~exist('BatchColocExclude', 'var')
            else
                filterNames = fieldnames(BatchColocExclude);
                for i = 1:numel(filterNames) % loop over named exclusion filtered objects
                    fname = filterNames{i};
                    % Grab image from each field
                    bwKeep = BatchColocExclude.(fname).ExcludeChan;
                    bwIf = BatchColocExclude.(fname).IfChanExcl;
                    mode = BatchColocExclude.(fname).ExclMode;
                
                    % Ensure logical images
                    bwKeep = logical(bwKeep);
                    bwIf = logical(bwIf);
                
                    switch mode
                        case 'Voxels'
                            maskKept_vox = bwKeep & ~bwIf; % remove overlapping voxels only
                            ColocSizeFilteredImgs.([fname '_vox']) = maskKept_vox; % add to 'ColocSizeFilteredImgs'
                
                        case 'Objects'
                            % Remove entire objects that overlap anywhere
                            labeledKeep = bwlabeln(bwKeep); % returns a label matrix for connected components (objects) derived from 'ExcludeChan'
                            props = regionprops3(labeledKeep, 'VoxelIdxList'); % returns linear indices of voxels
                            maskKept_obj = false(size(bwKeep)); % create empty image to catch filtered objects
                            for objIdx = 1:height(props) % loop over voxel indices derived from 'ExcludeChan'
                                voxelIdx = props.VoxelIdxList{objIdx}; 
                                if ~any(bwIf(voxelIdx)) % only keep if no voxel overlap with 'IfChanExcl'
                                    maskKept_obj(voxelIdx) = true; % add filtered objects to empty image
                                end
                            end
                            ColocSizeFilteredImgs.([fname '_obj']) = maskKept_obj; % add exclusion filtered image to struct
                
                        case 'Both'
                            maskKept_vox = bwKeep & ~bwIf;
                            ColocSizeFilteredImgs.([fname '_vox']) = maskKept_vox;
                
                            labeledKeep = bwlabeln(bwKeep);
                            props = regionprops3(labeledKeep, 'VoxelIdxList');
                            maskKept_obj = false(size(bwKeep));
                            for objIdx = 1:height(props)
                                voxelIdx = props.VoxelIdxList{objIdx};
                                if ~any(bwIf(voxelIdx))
                                    maskKept_obj(voxelIdx) = true;
                                end
                            end
                            ColocSizeFilteredImgs.([fname '_obj']) = maskKept_obj;
                    end
                    clear bwIf bwKeep maskKept_obj maskKept_vox labeledKeep props;
                end
            end
                clear BatchColocExclude;
        end

        if exist('BatchColocInclude', 'var') && isstruct(BatchColocInclude)
            fprintf('Applying inclusion filtering \n')
            
            % Retrieve masks for inclusion filtering
            fields1 = fieldnames(BatchColocInclude);
            for i = 1:numel(fields1)
                f1 = fields1{i};
                fields2 = fieldnames(BatchColocInclude.(f1));
                for j = 1:numel(fields2)
                    f2 = fields2{j};
                    chanName = BatchColocInclude.(f1).(f2);
                    if isfield(ColocSizeFilteredImgs, chanName)
                        BatchColocInclude.(f1).(f2) = ColocSizeFilteredImgs.(chanName);
                    else
                        warning('Channel "%s" not found in ColocSizeFilteredImgs.', chanName);
                    end
                end
            end
                
            % Filter by inclusion criteria
            if ~exist('BatchColocInclude', 'var')
            else 
                filterNames = fieldnames(BatchColocInclude);
                for i = 1:numel(filterNames)
                    fname = filterNames{i};
                    bwKeep = BatchColocInclude.(fname).KeepChan;
                    bwIf = BatchColocInclude.(fname).IfChan;
                    labeledKeep = bwlabeln(bwKeep);
                    props = regionprops3(labeledKeep, 'VoxelIdxList');
                    maskKept = false(size(bwKeep));
                    for objIdx = 1:height(props)
                        voxelIdx = props.VoxelIdxList{objIdx};
                        if any(bwIf(voxelIdx)) % inverse of exclusion filtering from above
                            maskKept(voxelIdx) = true;
                        end
                    end
                    ColocSizeFilteredImgs.(fname) = maskKept;
                end
                clear bwIf bwKeep maskKept labeledKeep props;
            end
            clear BatchColocInclude;
        end
        
        if exist('BatchColocComp', 'var') && isstruct(BatchColocComp)
            fprintf('Creating composite structures \n')
            
            % Retrieve masks for multi-channel composite objects
            fields1 = fieldnames(BatchColocComp);
            for i = 1:numel(fields1)
                f1 = fields1{i};
                if ~isfield(BatchColocComp.(f1), 'Channels')
                    continue
                end
                fields3 = fieldnames(BatchColocComp.(f1).Channels);
                for j = 1:numel(fields3)
                    f3 = fields3{j};
                    chanName = BatchColocComp.(f1).Channels.(f3);
                    if isfield(ColocSizeFilteredImgs, chanName)
                        BatchColocComp.(f1).Channels.(f3) = ColocSizeFilteredImgs.(chanName);
                    else
                        warning('Channel "%s" not found in ColocSizeFilteredImgs.', chanName);
                    end
                end
            end
                
            % Multi-Channel Composite Objects
            if ~exist('BatchColocComp', 'var')
            else
                ColocCompImgs = struct(); % create empty struct to catch multi-chan composite images created below
                subNames = fieldnames(BatchColocComp); % extract names of user-defined multi-chan composites
                
                for i = 1:numel(subNames) % loop over user-defined multi-chan composites
                    subName = subNames{i};
                    subStruct = BatchColocComp.(subName); % extract user-defined criteria for each multi-chan composite as separate struct
                    channelsStruct = subStruct.Channels; % extract channels for composite
                    combineType = lower(strtrim(subStruct.Combine)); % extract user-defined rule ('intersect', 'union', or 'both')
                    channelNames = fieldnames(channelsStruct);
                    numChannels = numel(channelNames);
                
                    if numChannels < 2
                        warning('Skipping %s: fewer than 2 channels found.', subName);
                        continue;
                    end
                
                    colocIntersect = channelsStruct.(channelNames{1}); % grab first channel to combine for 'inersect' case
                    colocUnion = colocIntersect; % grab first channel for 'union' case same way
                    for c = 2:numChannels % loop over remaining channels
                        thisChannel = channelsStruct.(channelNames{c}); % grab next channel
                        colocIntersect = colocIntersect & thisChannel; % combine iteratively using intersect rule
                        colocUnion = colocUnion | thisChannel; % combine iteratively using union rule
                    end
                
                    % Add desired composite images to ColocCompImgs struct
                    switch combineType
                        case 'intersect'
                            ColocCompImgs.([subName '_Intersect']) = colocIntersect;
                
                        case 'union'
                            ColocCompImgs.([subName '_Union']) = colocUnion;
                
                        case 'both'
                            ColocCompImgs.([subName '_Intersect']) = colocIntersect;
                            ColocCompImgs.([subName '_Union']) = colocUnion;
                    end
                    clear colocIntersect colocUnion thisChannel channelsStruct subStruct
                end
            
                % Add composite images from 'ColocCompImgs' to appropriate subfield of main struct 'ColocSizeFilteredImgs'
                for i = 1:numel(ColocSizeFilteredImgs)
                    fields = fieldnames(ColocCompImgs(i));
                    for f = 1:numel(fields)
                        ColocSizeFilteredImgs(i).(fields{f}) = ColocCompImgs(i).(fields{f});
                    end
                end
                clear ColocCompImgs;
            end
            clear BatchColocComp;
        end
    end

    % Find objects in other channels and composites
    if (Q.ColocAns == 1)
        ColocConnectedComponents = struct();
        channels = fieldnames(ColocSizeFilteredImgs);
        for i = 1:numel(channels)
            chan = channels{i};
            fprintf('Finding connected components in: %s\n', chan);
            binImg = ColocSizeFilteredImgs.(chan);
            CC = bwconncomp(binImg, 26);
            ColocConnectedComponents.(chan) = CC;
            clear binImg
        end
    
        ColocObjects = struct();   
        for i = 1:numel(channels)
            chan = channels{i};
            ColocObjects.(chan) = ColocConnectedComponents.(chan).PixelIdxList;
        end

        % Compute per-object volumes (physical units) per channel
        ColocObjectVolumes = struct();
        for i = 1:numel(channels)
            chan = channels{i};
            pixelLists = ColocObjects.(chan);
            volumes = cellfun(@(idx) numel(idx) * voxscale, pixelLists)';  % one volume per object
            ColocObjectVolumes.(chan) = volumes;
        end
    
        % Get image level-data about objects in other channels and add to ImageDataTable
        objCountFields = cellfun(@(f) ['ObjCount_', f], channels, 'UniformOutput', false);
        objDensityFields = cellfun(@(f) ['ObjDensity_', f], channels, 'UniformOutput', false);
        totalVolFields   = cellfun(@(f) ['TotalVol_', f], channels, 'UniformOutput', false);
        avgVolFields     = cellfun(@(f) ['AvgVol_', f], channels, 'UniformOutput', false);
    
        objCounts = cellfun(@(f) numel(ColocObjects.(f)), channels)';
        objDensities = objCounts / ImageDataTable.ImageVol;
        totalVols    = cellfun(@(f) sum(ColocObjectVolumes.(f)), channels)';
        avgVols      = cellfun(@(f) mean(ColocObjectVolumes.(f)), channels)';
    
        allFields  = [objCountFields(:); objDensityFields(:); totalVolFields(:); avgVolFields(:)]';
        ObjCountTable = array2table([objCounts, objDensities, totalVols, avgVols], 'VariableNames', allFields);
        ImageDataTable = [ImageDataTable, ObjCountTable];
        
    end

    % Retrieve user-defined objects and rules for colocalizing with glia
    if (Q.ColocAns == 1)
        fields1 = fieldnames(BatchColocObjectsSelect);
        for i = 1:numel(fields1)
            f1 = fields1{i};
            chanName = BatchColocObjectsSelect.(f1);
            if isfield(ColocObjects, chanName)
                BatchColocObjectsSelect.(f1) = ColocObjects.(chanName);
            else
                warning('Channel "%s" not found in ColocObjects.', chanName);
            end
        end
    end

    % Colocalize objects with glia structure
     if (Q.ColocAns == 1)
        clear ColocObjects;
    
        ColocResults = table(CellDataTable.CellNum, 'VariableNames', {'CellNum'});
        channels = fieldnames(BatchColocObjectsSelect);
        numChannels = numel(channels);
        ChannelMasks = struct();
        ColocOverlapIndices = struct();
    
        for i = 1:numChannels
            chan = channels{i};
            mask = false(size(ColocSizeFilteredImgs.(chan)));
            for j = 1:numel(BatchColocObjectsSelect.(chan))
                mask(BatchColocObjectsSelect.(chan){j}) = true;
            end
            ChannelMasks.(chan) = mask;
            clear mask;
        end
    
        % Single Channel Ovelaps
        fprintf('Finding single channel overlaps \n');
        for i = 1:numChannels
            chan = channels{i};
            mask = ChannelMasks.(chan);
            channelObjects = BatchColocObjectsSelect.(chan);
        
            OverlapVolumes = zeros(numObjMg, 1);
            OverlapCounts  = zeros(numObjMg, 1);
            OverlapIndices = [];
        
            for objIdx = 1:numObjMg
                objVoxels = FullMg{objIdx};
                objSize = numel(objVoxels);
        
                % Overlap mask
                overlapping = mask(objVoxels);
                overlapCount = nnz(overlapping);
                OverlapVolumes(objIdx) = overlapCount * voxscale;
                
                overlappingObjects = cellfun(@(x) any(ismember(x, objVoxels)), channelObjects);
                if any(overlappingObjects)
                    connectedVoxels = vertcat(channelObjects{overlappingObjects});
                    OverlapIndices = [OverlapIndices; connectedVoxels(:)];
                end
                OverlapCounts(objIdx) = sum(overlappingObjects);
                clear overlappingObjects overlapping objVoxels connectedVoxels;
            end
        
            % Store overlap summary stats in table
            ColocResults.([chan '_ColocVolume_um3']) = OverlapVolumes;
            ColocResults.([chan '_ColocNumObjects']) = OverlapCounts;
        
            % Store overlap voxel indices
            ColocOverlapIndices.(chan) = unique(OverlapIndices);
            clear OverlapIndices mask;
        end
        
        % Multi-Channel Overlaps
        if Q.ColocAllAns == 1
            fprintf('Finding multi-channel overlaps \n');
            numCombos = 0;
            
            for i = 2:numChannels
                combos = nchoosek(1:numChannels, i);
                numCombos = numCombos + size(combos, 1);
            end
            
            comboCount = 0;
            for i = 2:numChannels
                combos = nchoosek(1:numChannels, i);
                
                for comboIdx = 1:size(combos, 1)
                    comboCount = comboCount + 1;
                    idxs = combos(comboIdx, :);
                    comboNames = channels(idxs);
                    label = strjoin(comboNames, '_AND_');
        
                    % Build 'AND' mask
                    andMask = true(size(ChannelMasks.(channels{1})));
                    for j = 1:numel(comboNames)
                        andMask = andMask & ChannelMasks.(comboNames{j});
                    end
        
                    volumes = zeros(numObjMg, 1);
                    OverlapCounts = zeros(numObjMg, 1);
                    OverlapIndices = [];
        
                    for objIdx = 1:numObjMg
                        objVoxels = FullMg{objIdx};
                        overlapVoxels = andMask(objVoxels);
                        overlapCount = nnz(overlapVoxels);
        
                        volumes(objIdx)  = overlapCount * voxscale;
        
                        if overlapCount > 0
                            % Record indices of all overlapping voxels for this combo
                            OverlapIndices = [OverlapIndices; objVoxels(overlapVoxels)];
                            OverlapCounts(objIdx) = 1;  % at least one overlap
                        end
                        clear objVoxels overlapVoxels;
                    end
                    clear andMask;
        
                    % Store results
                    ColocResults.([label '_ColocVol_um3']) = volumes;
                    ColocResults.([label '_ColocNum']) = OverlapCounts;
        
                    % Save overlap voxel indices for this combo
                    ColocOverlapIndices.(label) = unique(OverlapIndices);
                    clear OverlapIndices;
                end
            end
        end
    end
    
    if exist('ColocResults', 'var') && istable(ColocResults)
        CellDataTable = join(CellDataTable, ColocResults, 'Keys', 'CellNum', 'KeepOneCopy', 'CellNum');
    end

    if (Q.ColocAns == 1)
        if exist('BatchColocOrigImgs', 'var') && isstruct(BatchColocOrigImgs)
            
            % Retrieve desired channels for saving original images
            fields1 = fieldnames(BatchColocOrigImgs);
            
            for i = 1:numel(fields1)
                f1 = fields1{i};
                chanName = BatchColocOrigImgs.(f1);
    
                if isfield(BatchColocObjectsSelect, chanName)
                    BatchColocOrigImgs.(f1) = BatchColocObjectsSelect.(chanName);
                else
                    warning('Channel "%s" not found in ColocObjects.', chanName);
                end
            end
        end

        if exist('BatchColocOrigImgs', 'var') && isstruct(BatchColocOrigImgs)
   
            % Retrieve desired channels for saving overlap images
            fields1 = fieldnames(BatchColocOverlapImgs);
            for i = 1:numel(fields1)
                f1 = fields1{i};
                chanName = BatchColocOverlapImgs.(f1);
                if isfield(ColocOverlapIndices, chanName)
                    BatchColocOverlapImgs.(f1) = ColocOverlapIndices.(chanName);
                else
                    warning('Channel "%s" not found in ColocObjects.', chanName);
                end
            end
        end
    end

    % Colocalization Images
    if (Q.ColocAns == 1)
        if (ColocSaveModeOverlay ~= 0)
            MergedMg = zeros(sz);
            for i = 1:length(FullMg)
                MergedMg(FullMg{1, i}) = 1;
            end
            
            channels = fieldnames(BatchColocOrigImgs);
            MergedOrigImgs = struct();
            for f = 1:numel(channels)
                field = channels{f};
                merged = zeros(sz);
                for i = 1:numel(BatchColocOrigImgs)
                    indices = BatchColocOrigImgs(i).(field);
                    for c = 1:numel(indices)
                        inds = indices{c};
                        merged(inds) = 1;
                    end
                    clear indices inds;
                end
                MergedOrigImgs.(field) = merged;
                clear merged;
            end
        end
    end
    
    % 2D All Objects From Other Channels
    if (Q.ColocAns == 1)
        if (ColocSaveModeOrig == 2) || (ColocSaveModeOrig == 1)
            channels = fieldnames(MergedOrigImgs);
            fprintf('Rendering in 2D \n');
            for i = 1:numel(channels)
                field = channels{i};
                vol = MergedOrigImgs.(field);
                mipZ = max(vol, [], 3);
                mipZ = cat(3, mipZ, zeros(size(mipZ)), zeros(size(mipZ)));
                fig = figure('Visible', 'off');
                imagesc(mipZ);
                daspect([1 1 1]);
                filename = fullfile(SaveFolderPath, [name '_step5_' field '_Obj_2D.png']);
                exportgraphics(fig, filename, 'Resolution', 300);
                clear mipZ vol;
            end
        end
    end
    
    % 3D All Objects From Other Channels (long processing time)
    if (Q.ColocAns == 1)
        if (ColocSaveModeOrig == 3) || (ColocSaveModeOrig == 1)
            channels = fieldnames(MergedOrigImgs);
            fprintf('Rendering in 3D (this may take a while) \n');
            for i = 1:numel(channels)
                field = channels{i};
                vol = MergedOrigImgs.(field);
                fv = isosurface(vol, 0.5);
            
                fig = figure('Visible', 'off');
                patch(fv, ...
                    'FaceColor', [1, 0, 0], ...
                    'FaceAlpha', 0.8, ...
                    'EdgeColor', 'none');
            
                camlight; lighting gouraud;
                view(0, -90);
                axis equal tight off;
                daspect([1 1 1]);   
                savefig(fig, fullfile(SaveFolderPath, [name '_step5_' field '_Obj_3D.fig']));
                filename = fullfile(SaveFolderPath, [name '_step5_' field '_Obj_3D.png']);
                exportgraphics(fig, filename, 'Resolution', 300);

                clear vol fv;
            end
        end
    end
    
    % 2D Overlapping Objects Only
    if (Q.ColocAns == 1)
        if (ColocSaveModeOverlay == 2) || (ColocSaveModeOverlay == 1)
            MgZ = double(max(MergedMg, [], 3) > 0);
            channels = fieldnames(BatchColocOverlapImgs);
            for i = 1:numel(channels)
                chan = channels{i};
                inds = BatchColocOverlapImgs.(chan);
    
                if isempty(inds)
                    warning('No overlap indices for channel: %s', chan);
                    continue;
                end
                mask = false(sz);
                mask(inds) = true;
                mip = max(mask, [], 3);
                overlayRGB = zeros([size(mip), 3]);
                overlayRGB(:,:,1) = mip;
                overlayRGB(:,:,2) = MgZ;
                fig = figure('Visible', 'off');
                imagesc(overlayRGB);
                daspect([1 1 1]);
                
                % Legend
                annotation('textbox', [0.27 0.78 0.15 0.15], 'String', ...
                    {['\color{red}' chan], ...
                     '\color{green}Microglia', ...
                     '\color{yellow}Overlap'}, ...
                    'FitBoxToText', 'on', 'EdgeColor', 'none', 'Color', 'w', 'FontSize', 10);
                
                filename = fullfile(SaveFolderPath, [name '_step6_' chan '_Overlap_2D.png']);
                exportgraphics(fig, filename, 'Resolution', 300);
                clear inds mask;
            end
        end
    end
    
    % 3D Overlapping Objects Only
    if (Q.ColocAns == 1)
        if (ColocSaveModeOverlay == 3) || (ColocSaveModeOverlay == 1)
            channels = fieldnames(BatchColocOverlapImgs);
            MgMask = logical(MergedMg);
            fprintf('Rendering overlaps in 3D (this may take a while) \n');
            for i = 1:numel(channels)
                ch = channels{i};
                inds = BatchColocOverlapImgs.(ch);
                    
                % Convert linear indices to 3D logical mask
                merged = false(sz);
                merged(inds) = true;
                merged = logical(merged);

                % Skip empty volumes
                if nnz(merged) == 0
                    warning('Skipping empty channel: %s', ch);
                    continue;
                end
            
                % Get overlapping and exclusive voxels
                fprintf(['Getting overlaps for ' ch '... \n']);
                overlap = merged & MergedMg;
                onlyMerged = merged & ~overlap;
                onlyMergedMg = MergedMg & ~overlap;
            
                % Compute isosurfaces and reduce mesh complexity (speed up processing)
                fvMerged = isosurface(onlyMerged, 0.5);
                if numel(fvMerged.faces) > 1000
                    fvMerged = reducepatch(fvMerged, 0.3);
                end
                fvMerged.vertices = double(fvMerged.vertices);
                fvMergedMg = isosurface(onlyMergedMg, 0.5);
                fvOverlap = isosurface(overlap, 0.5);
            
                % Create figure
                fig = figure('Visible', 'off');
            
                % Render other channel objects in red
                fprintf(['Rendering ' ch ' (Red) \n']);
                if ~isempty(fvMerged.vertices)
                    patch(fvMerged, ...
                        'FaceColor', [1, 0, 0], ...
                        'FaceAlpha', 0.3, ...
                        'EdgeColor', 'none');
                end
            
                % Render microglia in green
                fprintf('Rendering MergedMg (Green) \n');
                if ~isempty(fvMergedMg.vertices)
                    patch(fvMergedMg, ...
                        'FaceColor', [0, 1, 0], ...
                        'FaceAlpha', 0.3, ...
                        'EdgeColor', 'none');
                end
            
                % Render overlap in yellow
                fprintf('Rendering Overlap (Yellow) \n');
                if ~isempty(fvOverlap.vertices)
                    patch(fvOverlap, ...
                        'FaceColor', [1, 1, 0], ...
                        'FaceAlpha', 1, ...
                        'EdgeColor', 'none');
                end
            
                camlight; lighting gouraud;
                axis equal tight off;
                view(0, -90);
                daspect([1 1 1]);
            
                % Legend
                legendEntries = {};
                legendEntries{end+1} = 'Microglia (Green)';
                legendEntries{end+1} = [ch ' (Red)'];
                legendEntries{end+1} = 'Overlap (Yellow)';
            
                if ~isempty(legendEntries)
                    legend(legendEntries, 'TextColor', 'w', 'Location', 'northeastoutside');
                end
                 
                savefig(fig, fullfile(SaveFolderPath, [name '_step6_' ch '_Overlap_3D.fig']));
                filename = fullfile(SaveFolderPath, [name '_step6_' field '_Overlap_3D.png']);
                exportgraphics(fig, filename, 'Resolution', 300);
                
                clear inds merged overlap onlyMerged onlyMergedMg fvMerged fvMergedMg fvOverlap;
                
            end
        end
    end

    % Clean up workspace before skeletonization step
    clear allObjs binImg colocObjs CC ChannelMasks ColocBinaryThresholdedImgs ColocConnectedComponents ...
        ColocImgs ColocComp ColocCompImgs ColocImgsNorm colocIntersect colocUnion ColocObjects ColocObjectsSelect ColocOrigImgs ...
        ColocOverlapImgs ColocOverlapIndices ColocSizeFilteredImgs channelsStruct depthMap DetectedObjs ex ...
        flatex fullimg imgChannels img imgGlia imgGliaNorm img2D img3D mask merged MergedMg MergedOrigImgs MgZ ... 
        Microglia mip NoiseIm onlyMerged onlyMergedMg overlap overlapping ...
        overlappingObjects overlayRGB OverlapIndices ProxCC subStruct thisChannel x y ...
        BatchColocOverlapImgs BatchColocOrigImgs BatchColocObjectsSelect BatchColocExclude BatchColocInclude BatchColocComp BatchProxCC;

    %% SKELETONIZATION

    % must recalculate centroids
    numCells = numel(FullMg);
    centroids = zeros(numCells, 3);
    
    fprintf('Calculating centroids...\n');
    
    for i = 1:numCells
        % Calculate centroid directly from indices, no 3D array needed
        [x, y, z] = ind2sub([sizeXY(1), sizeXY(2), sizeZ], FullMg{1,i});
        centroids(i, 1) = mean(x);
        centroids(i, 2) = mean(y);
        centroids(i, 3) = mean(z);
    end  

    if Q.SkelContinueAns == 1 || Q.SkelContinueAns == 3
    
        kernel(:,:,1) = [1 1 1; 1 1 1; 1 1 1];
        kernel(:,:,2) = [1 1 1; 1 0 1; 1 1 1];
        kernel(:,:,3) = [1 1 1; 1 1 1; 1 1 1]; 
        
        numendpts = zeros(numel(FullMg),1);
        numbranchpts = zeros(numel(FullMg),1);
        MaxBranchLength = zeros(numel(FullMg),1);
        MinBranchLength = zeros(numel(FullMg),1);
        AvgBranchLength = zeros(numel(FullMg),1);
        
        BranchLengthList=cell(1,numel(FullMg));
        
        % Trace skeletons
        % Preallocate to save info from parfor loop
        SkelErrorLog = cell(1, numel(FullMg)); % Error log for debugging later
        Skels(1, numel(FullMg)) = struct('ID', [], 'Mask', [], 'Status', '');
        adjust_pxXY_list = zeros(1, numel(FullMg));
        
        parfor i=1:numel(FullMg)
            try
               % Convert linear indices for this cell to subscripts
                [idxX, idxY, idxZ] = ind2sub([sizeXY sizeZ], FullMg{i});
                
                % Tight bounding box
                minX = min(idxX);  maxX = max(idxX);
                minY = min(idxY);  maxY = max(idxY);
                minZ = min(idxZ);  maxZ = max(idxZ);
                
                % add padding to bounding box
                pad = 5; % change as needed, but 5 works for now
                minX_p = max(minX - pad, 1);
                maxX_p = min(maxX + pad, sizeXY(1));
                minY_p = max(minY - pad, 1);
                maxY_p = min(maxY + pad, sizeXY(2));
                minZ_p = max(minZ - pad, 1);
                maxZ_p = min(maxZ + pad, sizeZ);
                
                % padded size
                bbSize = [maxX_p-minX_p+1, maxY_p-minY_p+1, maxZ_p-minZ_p+1];
                
                % local coordinates
                locX = idxX - minX_p + 1;
                locY = idxY - minY_p + 1;
                locZ = idxZ - minZ_p + 1;
                
                ex = false(bbSize);
                linLoc = sub2ind(bbSize, locX, locY, locZ);
                ex(linLoc) = true;
                ex = single(ex);
                
                % store padded offsets
                Skels(i).Offset = struct( ...
                    'left',   minX_p, ...
                    'right',  maxX_p, ...
                    'top',    minY_p, ...
                    'bottom', maxY_p, ...
                    'zmin',   minZ_p, ...
                    'zmax',   maxZ_p );
        
            if Q.SkelModeAns == 1
                WholeSkel = SlimSkel3D(ex,100);
                DownSampled = 0;
                adjust_pxXY = pxXY;
            end
            
            if Q.SkelModeAns == 2
                if sizeXY(1)>512 %convert large images to 512x512 to speed up skeletonization. The branch lengths are later adjusted to account for this down-sampling.
                    ex = imresize(ex,0.5);
                    adjust_pxXY = 2*pxXY;
                    DownSampled = 1;
                else 
                    DownSampled = 0;
                    adjust_pxXY = pxXY;
                end
        
                adjust_pxXY_list(i) = adjust_pxXY; %store adjusted scaling in case of downsampling, you need this for Sholl downstream        
                SmoothEx = imgaussfilt3(ex); %Smooth the cell so skeleton doesn't pick up many fine hairs
                FastMarchSkel = skeleton(SmoothEx);%Find the skeleton! This uses msfm3d and rk4 files, which have been compiled and the .mexw64 versions included. If errors, re-run compilation of these files (in FastMarching_version3b folder), and add the folder and subfolders to path. 
        
                %Convert cell output of branches into one image for further processing.
                WholeSkel=zeros(size(ex));
                WholeList = round(vertcat(FastMarchSkel{:}));
                SkelIdx = sub2ind(size(ex),WholeList(:,1),WholeList(:,2),WholeList(:,3));
                WholeSkel(SkelIdx)=1;
            end
        
            % Create bounding box 
            BoundedSkel = WholeSkel;
            si = size(WholeSkel);
            Skels(i).Size = si;
        
            % Find centroid inside skeleton (start point for branch tracing)
            globalCentroid = floor(centroids(i,:));
        
            % convert to local bounding-box coordinates
            i2 = [ ...
                globalCentroid(1) - minX_p + 1, ...
                globalCentroid(2) - minY_p + 1, ...
                globalCentroid(3) - minZ_p + 1 ];
                
            % Adjust centroid for downsampling
            if DownSampled == 1
                i2(1) = round(i2(1)/2);
                i2(2) = round(i2(2)/2);
            end
            
            % Snap centroid to nearest skeleton voxel
            i2 = NearestPixel(WholeSkel, i2, pxXY);
        
            % Find endpoints
            endpts = (convn(BoundedSkel, kernel, 'same') == 1) & BoundedSkel;
            Skels(i).EndpointIdx = find(endpts);
        
            EndptList = Skels(i).EndpointIdx;
            [r,c,p] = ind2sub(si, EndptList);
            EndptList = [r c p];
            numendpts(i,:) = length(EndptList);
        
            % trace branches
            ArclenOfEachBranch = zeros(length(EndptList),1);
            fullmask = zeros(si, 'uint8');
        
            for j=1:length(EndptList)
                i1 = EndptList(j,:);
                mask = ConnectPointsAlongPath(BoundedSkel, i1, i2);
                if nnz(mask) < 2
                    % Invalid or degenerate path → skip to next branch
                    ArclenOfEachBranch(j) = 0;
                    continue;
                end
        
                fullmask(mask) = fullmask(mask) + 1;
        
                % Arc length
                pxlist = find(mask);
                distpoint = reorderpixellist(pxlist, si, i1, i2);
                distpoint(:,1) = distpoint(:,1) * adjust_pxXY;
                distpoint(:,2) = distpoint(:,2) * adjust_pxXY;
                distpoint(:,3) = distpoint(:,3) * pxZ;
                ArclenOfEachBranch(j) = arclength(distpoint(:,1), distpoint(:,2), distpoint(:,3));
            end
        
            % Save branch lengths
            BranchLengthList{1,i} = ArclenOfEachBranch;
            if ~isempty(ArclenOfEachBranch)
                AvgBranchLength(i) = mean(ArclenOfEachBranch);
            else
                AvgBranchLength(i) = 0;
            end
        
            % Decode branch levels
            tmp = fullmask; 
            tmp(tmp>3)=4;
        
            pri  = (tmp==4);
            sec  = (tmp==3);
            tert = (tmp==2);
            quat = (tmp==1);
            
            Skels(i).PrimaryIdx     = find(pri);
            Skels(i).SecondaryIdx   = find(sec);
            Skels(i).TertiaryIdx    = find(tert);
            Skels(i).QuaternaryIdx  = find(quat);
        
            % Find branchpoints
            brpts = zeros(si(1), si(2), si(3), 4);
            for kk = 1:3
                temp = fullmask > kk;    % uint8 > double → logical
                tempendpts = (convn(double(temp), kernel, 'same') == 1) & temp;
                brpts(:,:,:,kk+1) = tempendpts;
            end
            
            quatendpts = (convn(double(quat), kernel, 'same') == 1) & quat;
            quatbrpts  = quatendpts - endpts;
            
            % Keep only branchpoints connecting to level-4
            fullrep = double(fullmask);
            fullrep(fullrep < 4) = 0; % keep only primary branches
            qbpts  = fullrep + double(quatbrpts);
            qbpts1 = convn(qbpts, ones([3 3 3]), 'same');
            brpts(:,:,:,1) = (double(quatbrpts) .* qbpts1) >= 5;
            
            % Combine and count
            allbranch = sum(brpts,4);
            BranchptList = find(allbranch);
            numbranchpts(i) = numel(BranchptList);
        
            Skels(i).BranchpointIdx = find(allbranch);
            Skels(i).Status = 'success';
        
            catch ME
                SkelErrorLog{i} = ME;
                Skels(i).Status = 'error';
                BranchLengthList{1,i} = 'failed';
                
                % or leave blank to do nothing if an error is detected, just write in zeros and
                %continue to next loop iteration. 
           end
           
           fprintf(['completed cell ' num2str(i) ' of ' num2str(numel(FullMg)) '\n']);
        end
        
        delete(gcp('nocreate')); %close parallel pool so error isn't generated when program is run again.
        
        % Skeleton Tracing Error Log
        for i = 1:numel(SkelErrorLog)
            if ~isempty(SkelErrorLog{i})
                fprintf('Error in cell %d:\n', i);
                fprintf(SkelErrorLog{i}.getReport('extended', 'hyperlinks', 'off'));
                fprintf('--------------------------------------------------\n');
            end
        end
        
        % Calculate total branch length
        
        nCells = numel(FullMg);
        TotalBranchLengths = zeros(nCells, 1);
        
        for RowNum = 1:nCells
            branch_lengths = BranchLengthList{1, RowNum};
            if ~isempty(branch_lengths)
                TotalBranchLengths(RowNum) = sum(branch_lengths);
            else
                TotalBranchLengths(RowNum) = 0;
            end
        end
        
        % Save images of skeletons
        % this step has been moved outside parfor loop to make customization easier
        for i = 1:length(Skels)
            if ~isfield(Skels(i), 'Status') || ~strcmp(Skels(i).Status, 'success')
                continue;
            end
        
            si = Skels(i).Size;  
            pri = false(si);  pri(Skels(i).PrimaryIdx) = true;
            sec = false(si);  sec(Skels(i).SecondaryIdx) = true;
            tert = false(si); tert(Skels(i).TertiaryIdx) = true;
            quat = false(si); quat(Skels(i).QuaternaryIdx) = true;
            
            branchpts = false(si); branchpts(Skels(i).BranchpointIdx) = true;
            endpts = false(si); endpts(Skels(i).EndpointIdx) = true;
        
            figTitle = [name, '_Cell', num2str(i)];
            figure('Name', figTitle, 'Visible', 'off');
        
            hold on
        
            % Plot primary branches (purple)
            if any(pri(:))
                fv1 = isosurface(pri, 0);
                patch(fv1, 'FaceColor', [1 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                camlight
                lighting gouraud
            end
        
            % Plot secondary branches (yellow)
            if any(sec(:))
                fv1 = isosurface(sec, 0);
                patch(fv1, 'FaceColor', [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                camlight
                lighting gouraud
            end
        
            % Plot tertiary branches (green)
            if any(tert(:))
                fv1 = isosurface(tert, 0);
                patch(fv1, 'FaceColor', [0 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                camlight
                lighting gouraud
            end
        
            % Plot quaternary branches (cyan)
            if any(quat(:))
                fv1 = isosurface(quat, 0);
                patch(fv1, 'FaceColor', [0 1 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                camlight
                lighting gouraud
            end
            
            % Plot Branchpoints (red)
            if any(branchpts(:))
                fv2 = isosurface(branchpts, 0);
                patch(fv2, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
                camlight
                lighting gouraud
            end
        
            % Plot Endpoints (red)
            if any(endpts(:))
                fv2 = isosurface(endpts, 0);
                patch(fv2, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
                camlight
                lighting gouraud
            end
        
            % Final figure settings
            view(0, 270);
            daspect([1 1 1]);
            hold off
        
            % Save image
            SubFolderSkels = fullfile(SaveFolderPath, 'Individual_Skeletons');
            if ~exist(SubFolderSkels, 'dir')
                mkdir(SubFolderSkels);
            end
            filename = sprintf('%s_Skeleton_cell%d.png', name, i);
            exportgraphics(gcf, fullfile(SubFolderSkels, filename), 'Resolution', 300);
            close(gcf);
        end

    %% FINAL DATA OUTPUT

    CellDataTable.AvgBranchLength_um = AvgBranchLength(:);
    CellDataTable.NumEndpoints = numendpts(:);
    CellDataTable.NumBranchpoints = numbranchpts(:);
    CellDataTable.TotalBranchLength_um = TotalBranchLengths;

    end
    
    nCells = height(CellDataTable); 
    CellDataTable.FileName = repmat({name}, nCells, 1);
    cols = CellDataTable.Properties.VariableNames;
    CellDataTable = CellDataTable(:, [{'FileName'}, setdiff(cols, {'FileName'}, 'stable')]);
    csvFilename = fullfile(SaveFolderPath, ['CellData_', name, '.csv']);
    writetable(CellDataTable, csvFilename);

    nCells1 = height(ImageDataTable); 
    ImageDataTable.FileName = repmat({name}, nCells1, 1); % add a column with filename (useful for compiling multiple .csv files from batch processing)
    cols1 = ImageDataTable.Properties.VariableNames;
    ImageDataTable = ImageDataTable(:, [{'FileName'}, setdiff(cols1, {'FileName'}, 'stable')]); % move filename to first column
    csvFilename1 = fullfile(SaveFolderPath, ['ImageData_', name, '.csv']);
    writetable(ImageDataTable, csvFilename1);

    fprintf(['completed file ' num2str(k) ' of ' num2str(numel(imgFiles)) '\n']);
    fprintf('%s \n', datetime('now'));

    fprintf('----------------------------------------------------------------------------\n')

    % Clean up memory for next iteration
    clearvars -except BatchFolder fullScriptPath imgFiles ParametersFile ...
                  ParametersPath SaveBatchPath ScriptFolder FullMgFiles k;
    close all hidden;
    memory
    whos

    fprintf('----------------------------------------------------------------------------\n')

end

fprintf('----------------------------------------------------------------------------\n')
totalTime = toc;
time_duration = seconds(totalTime);
fprintf('Batch completed in %s (hh:mm:ss) \n', string(time_duration, 'hh:mm:ss'));
fprintf('Batch Processing Complete! \n');
diary off; % end log