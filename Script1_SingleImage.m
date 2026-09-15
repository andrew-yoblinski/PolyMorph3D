%% Required Add-Ons!!!
% You need to install the following add-ons if you don't have them already:
    % Image Processing Toolbox (MathWorks)
    % Parallel Computing Toolbox (Mathworks)
    % Statistics and Machine Learning Toolbox (Mathworks)

% Click 'Run' under the 'Editor' tab at the top, select "Add to Path" if prompted

%% Setup
clear; clc; %clear memory and command window log
fullScriptPath = mfilename('fullpath'); %get local path to Script1
ScriptFolder = fileparts(fullScriptPath); %define path to parent folder
addpath(genpath(ScriptFolder)); %add parent folder and all of its subfolders (genpath) to matlab path, you can also do this manually with built-in point and click interface
fprintf('Added %s and all its subfolders to the MATLAB path.\n', ScriptFolder);

%% Welcome Message
app = Welcome_App; %call app
waitfor(app, 'Output');
if strcmp(app.ChoiceContinue, 'Yes')
    delete(app);
elseif strcmp(app.ChoiceContinue, 'No')
    disp('Script aborted until the following required add-ons are installed:');
    disp('Image Processing Toolbox (MathWorks)');
    disp('Parallel Computing Toolbox (Mathworks)');
    disp('Statistics and Machine Learning Toolbox (Mathworks)');
    delete(app);
    return;
end

%% Select File
app = FileDataGUI_App;
waitfor(app, 'Output');
delete(app);

% Create save directory
[filepath, name, ext] = fileparts(file);
folderName = [name, '_data'];
SaveFolderPath = fullfile(SavePath, folderName);
if ~exist(SaveFolderPath, 'dir')
    mkdir(SaveFolderPath);
end

%% Import: Load tiff file
disp('Importing File and Metadata');
progbar = waitbar(0, 'Importing File and Metadata');

tifFile = fullfile(pathname, [name, '.tif']);
    waitbar(0.1, progbar);
if exist(tifFile, 'file') % check that the file exists
    reader = bfGetReader(tifFile); % use bioformats toolbox to read in the tiff file
        waitbar(0.25, progbar);
    reader.setSeries(0);

    % Get image dimensions
    sizeX = reader.getSizeX(); % pixels in X
    sizeY = reader.getSizeY(); % pixels in Y
    sizeZ = reader.getSizeZ(); % number of Z-slices
    sizeCZ = reader.getImageCount();
    sizeC = reader.getSizeC(); % number of channels
    sizeXY = [sizeY, sizeX];
    total_frames = sizeZ * sizeC;
        waitbar(0.9, progbar);

   if reader.getImageCount() ~= total_frames
       error('Mismatch between expected and actual frame count.');
   end

    % Extract voxel dimensions (user will double check in app below)
    omeMeta = reader.getMetadataStore();
    pxXY = omeMeta.getPixelsPhysicalSizeX(0);
    pxZ = omeMeta.getPixelsPhysicalSizeZ(0);

    % Ensure format can be handled by subsequent operations
    if isempty(pxXY)
        pxXY = 'N/A';
    else
        pxXY = double(pxXY.value()); % note that this will be presumed μm going forward
    end
    
    if isempty(pxZ)
        pxZ = 'N/A';
    else
        pxZ = double(pxZ.value()); % note that this will be presumed μm going forward
    end

    img = zeros(sizeY, sizeX, sizeZ, sizeC, 'uint16');
    for c = 1:sizeC
        for z = 1:sizeZ
            index = reader.getIndex(z - 1, c - 1, 0);
            plane = bfGetPlane(reader, index + 1);
            img(:, :, z, c) = plane;
            clear plane;
        end
    end

    % separate channels
    sizeC = size(img, 4);
    imgChannels = cell(1, sizeC);
    for c = 1:sizeC
        imgChannels{c} = squeeze(img(:, :, :, c));
    end
end
        waitbar(1, progbar);
        reader.close();

    if isgraphics(progbar)
        close(progbar);
    end

disp('File Successfully Imported');
clear img; % free memory

%% User check metadata
InfoInput.pxXY = pxXY;
InfoInput.pxZ = pxZ;
app = FileDataOverride_App;
setup(app, InfoInput);
waitfor(app, 'Output');
delete(app);
clear InfoInput;

% reminder that real world unit is presumed to be in μm
voxscale = pxXY*pxXY*pxZ; % define voxel size (volume in μm^3)
voxdim = [pxXY, pxXY, pxZ]; % define voxel dimensions
imgGlia = imgChannels{round(GliaChannel)};
MeanFI_Glia = mean(imgGlia, 'all'); % get mean fluorescent intensity of primary uglia channel

ImageDataTable = table(MeanFI_Glia, ...
    'VariableNames', {'MeanFI_Glia'}); %save MFI to data table containing image level information

% normalization slice by slice with global scaling
globalMin = double(min(imgGlia(:)));
globalMax = double(max(imgGlia(:)));
globalRange = globalMax - globalMin;
imgGliaNorm = zeros(size(imgGlia), 'single');
fprintf('Normalizing glia channel with global scaling...\n');
progbar = waitbar(0, 'Normalizing glia channel...');
for z = 1:size(imgGlia, 3)
    imgGliaNorm(:,:,z) = single((double(imgGlia(:,:,z)) - globalMin) / globalRange);
    if mod(z, 100) == 0
        waitbar(z/size(imgGlia,3), progbar);
    end
end
if isgraphics(progbar)
    close(progbar);
end
fprintf('Normalization complete.\n');

sz = size(imgGliaNorm);

%% Second glia channel if applicable
if SecondGliaChanYN == 1
    % normalize like primary glia channel
    imgGlia2 = imgChannels{round(GliaChannel2)};
    MeanFI_Glia_Ch2 = mean(imgGlia2, 'all'); % get mean fluorescent intensity
    ImageDataTable.MeanFI_Glia_Ch2 = MeanFI_Glia_Ch2(:); %save MFI for glia chan 2
    
    globalMin = double(min(imgGlia2(:)));
    globalMax = double(max(imgGlia2(:)));
    globalRange = globalMax - globalMin;
    imgGliaNorm2 = zeros(size(imgGlia2), 'single');
    fprintf('Normalizing glia channel with global scaling...\n');
    progbar = waitbar(0, 'Normalizing second glia channel...');
    for z = 1:size(imgGlia2, 3)
        imgGliaNorm2(:,:,z) = single((double(imgGlia2(:,:,z)) - globalMin) / globalRange);
        if mod(z, 100) == 0
            waitbar(z/size(imgGlia2,3), progbar);
        end
    end
    if isgraphics(progbar)
        close(progbar);
    end
    fprintf('Normalization of second glia channel complete.\n');

    % user decide how to use second glia channel
    app = GliaChannel2_App;
    setup(app);
    waitfor(app, 'Output');
    delete(app);
end

%% Threshold Glia Structure Channel
ThreshInput.orig = max(imgGlia,[], 3);
ThreshInput.imgGliaNorm = imgGliaNorm;
ThreshInput.pxXY = pxXY;
ThreshInput.pxZ = pxZ;
ThreshInput.length = sizeZ;
ThreshInput.size = sizeXY;
ThreshInput.voxscale = voxscale;
app = ThresholdGUI_App;
setup(app, ThreshInput);
waitfor(app, 'Output');
delete(app);

if SecondGliaChanYN == 1 
    ThreshInput.orig = max(imgGlia2,[], 3);
    ThreshInput.imgGliaNorm2 = imgGliaNorm2;
    app = ThresholdCh2GUI_App;
    setup(app, ThreshInput);
    waitfor(app, 'Output');
    delete(app);
end

clear ThreshInput;

%% Save max projection of original primary glia image
if  SaveOrig == 1
    img2D = max(imgGlia, [], 3);
    figName = [name,'_step0_Glia_OrigMAX'];
    fig = figure('Name',figName);
    imagesc(img2D);
    colormap gray;
    daspect([1 1 1]);
    filename = [name '_step0_Glia_OrigMAX' '.png'];
    exportgraphics(fig, fullfile(SaveFolderPath, filename), 'Resolution', 300);
end

% save max project of secondary glia image if applicable
if SecondGliaChanYN == 1
    if  SaveOrig2 == 1
        img2D = max(imgGlia2, [], 3);
        figName = [name,'_step0_GliaCh2_OrigMAX'];
        fig = figure('Name',figName);
        imagesc(img2D);
        colormap gray;
        daspect([1 1 1]);
        filename = [name '_step0_GliaCh2_OrigMAX' '.png'];
        exportgraphics(fig, fullfile(SaveFolderPath, filename), 'Resolution', 300);
    end
end
clear imgGlia imgGlia2 imgGliaNorm imgGliaNorm2 img2D % Free memory

if SecondGliaChanYN == 1 

    % only keep glia ch1 objects that overlap with at least 1 voxel in glia ch2
    CC = bwconncomp(NoiseIm,26);
    maskKept = false(size(NoiseIm));
    for objIdx = 1:CC.NumObjects
        voxelIdx = CC.PixelIdxList{objIdx};
        if any(NoiseIm2(voxelIdx))
            maskKept(voxelIdx) = true;
        end
    end

    if SaveOverlapImg == 1
        overlap3D = NoiseIm & NoiseIm2;
        
        Q.SaveGlia2 = {'How would you like to see the overlap of your two glia channels?', 
            ['Note: The 2D image may be misleading as voxels in the same XY location from different planes will be rendered as overlapping. ' ...
            'This is a quirk of the 2D visualization and does not affect the data analysis, which is always in 3D.'],
            'Warning: A 3D image will take longer to process!'};
        Q.SaveGlia2Choice = questdlg(Q.SaveGlia2,'Save Image of Cells','3D figure', '2D image', 'Both', '2D image');
        switch Q.SaveGlia2Choice
            case '3D figure'
                Q.SaveGlia2Ans = 3;
            case '2D image'
                Q.SaveGlia2Ans = 2;    
            case 'Both'
                Q.SaveGlia2Ans = 1;
        end
   
        % build and save visualization of glia channels overlap
        if Q.SaveGlia2Ans == 2 || Q.SaveGlia2Ans == 1
            % 2D visual
            R3D = double(NoiseIm2) | double(overlap3D);
            G3D = double(maskKept) | double(overlap3D);
            B3D = zeros(size(maskKept));
            R2D = max(R3D, [], 3);
            G2D = max(G3D, [], 3);
            B2D = max(B3D, [], 3);
            overlapRGB = cat(3, R2D, G2D, B2D);
        
            % save 2D image of overlaps
            figName = [name,'_step0_Glia2_Overlap_2D'];
            fig = figure('Name',figName);
            imagesc(overlapRGB);
            daspect([1 1 1]);
            fileName = [name '_step0_Glia2_Overlap_2D.png'];
            exportgraphics(fig, fullfile(SaveFolderPath, fileName), 'Resolution', 300);
        end
    
        % edit this so that the dimensions are fixed
        if Q.SaveGlia2Ans == 3 || Q.SaveGlia2Ans == 1
            fprintf('Rendering in 3D (this may take a few minutes)...\n');
            
            % 3D visualization
            fig = figure('Visible', 'off', 'Name', '_step0_Glia2_Overlap_3D');
            hold on;
        
            % Primary glia chan only (green)
            NoiseImOnly = maskKept & ~overlap3D;
            if any(NoiseImOnly(:))
                fv = isosurface(double(NoiseImOnly), 0.5);
                p1 = patch(fv);
                set(p1, 'FaceColor', [0 1 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
        
            % Secondary glia chan only (red)
            NoiseIm2Only = NoiseIm2 & ~overlap3D;
            if any(NoiseIm2Only(:))
                fv = isosurface(double(NoiseIm2Only), 0.5);
                p2 = patch(fv);
                set(p2, 'FaceColor', [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
        
            % Overlap (yellow)
            if any(overlap3D(:))
                fv = isosurface(double(overlap3D), 0.5);
                p3 = patch(fv);
                set(p3, 'FaceColor', [1 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
            end
            
            camlight; lighting gouraud;
            axis equal tight off;
            xlim([1 sz(2)]);
            ylim([1 sz(1)]);
            zlim([1 sz(3)]);
            view(0, -90);
            daspect([1 1 1]);
        
            % Legend
            legendEntries = {};
            legendEntries{end+1} = 'Primary Glia Chan (Green)';
            legendEntries{end+1} = 'Secondary Glia Chan (Red)';
            legendEntries{end+1} = 'Overlap (Yellow)';
        
            if ~isempty(legendEntries)
                legend(legendEntries, 'TextColor', 'w', 'Location', 'northeastoutside');
            end
                        
            % save 3D render
            figFilename = fullfile(SaveFolderPath, [name '_step0_Glia2_Overlap_3D.fig']);
            savefig(fig, figFilename);
            pngFilename = fullfile(SaveFolderPath, [name '_step0_Glia2_Overlap_3D.png']);
            exportgraphics(fig, pngFilename, 'Resolution', 300);
            set(fig, 'Visible', 'on');
            hold off;
        end
    end

    if SecondGliaChanChoice == 1
        NoiseIm = maskKept; % keep overlap without merging if previously selected
    
    elseif SecondGliaChanChoice == 2
        % keep overlap and merge glia channels if previously selected
        CC2 = bwconncomp(NoiseIm2, 26);
        maskKept2 = false(size(NoiseIm2));
        for objIdx = 1:CC2.NumObjects
            voxelIdx = CC2.PixelIdxList{objIdx};
            if any(maskKept(voxelIdx))
                maskKept2(voxelIdx) = true;
            end
        end
        NoiseIm = maskKept | maskKept2;
    end
end

progbar = waitbar(0,'Processing thresholded/filtered image...');
ConnectedComponents=bwconncomp(NoiseIm,26); %returns structure with 4 fields. PixelIdxList contains a 1-by-NumObjects cell array where the k-th element in the cell array is a vector containing the linear indices of the pixels in the k-th object. 26 defines connectivity. This looks at cube of connectivity around pixel.
numObj = numel(ConnectedComponents.PixelIdxList); %PixelIdxList is field with list of pixels in each connected component. Find how many connected components there are.

fprintf('Processing %d objects...\n', numObj);

for i = 1:numObj
    waitbar(i/numObj, progbar);
    ex=zeros(sizeXY(1),sizeXY(2),sizeZ);
    ex(ConnectedComponents.PixelIdxList{1,i})=1;%write in only one object to image. Cells are white on black background.
    flatex = sum(ex,3);
    allObjs(:,:,i) = flatex(:,:);
end
clear ex flatex NoiseIm;

if isgraphics(progbar)
    close(progbar);
end
fprintf('Object processing complete.\n');

DetectedObjs = sum(allObjs,3);
cmapCompress = hot(max(DetectedObjs(:)));  
cmapCompress(1,:) = zeros(1,3);

% Save max projection of thresholded glia image
if  ShowImg == 1    
    figName = [file,'_step1_Glia_Thresh'];
    fig = figure('Name', figName);
    imagesc(DetectedObjs);
    colormap(cmapCompress);
    daspect([1 1 1]);
    filename = ([name '_step1_Glia_Thresh' '.png']);
    exportgraphics(fig, fullfile(SaveFolderPath, filename), 'Resolution', 300);
end
clear DetectedObjs;

% Extract list of pixel values and which object they belong to for viewing (in CellSizeCutoffGUI). 
for i = 1:numObj
    ObjectList(i,1) = length(ConnectedComponents.PixelIdxList{1,i}); 
    ObjectList(i,2) = i;  
end

ObjectList = sortrows(ObjectList,-1);%Sort columns by pixel size.
udObjectList = flipud(ObjectList);%ObjectList is large to small, flip upside down so small is plotted first in blue.

%% Cell Segmentation 
% Decreased threshold may cause cells to be inappropriately connected.
% Choose the threshold for a large cell, and segment larger objects into
% separate cells (by identifying number of nuclei and running fitgmdist
% (fit Gaussian mixture distribution). Function is run 3 times to improve
% accuracy and replicability

SegInput.udObjectList = udObjectList;
SegInput.ObjectList = ObjectList;
SegInput.allObjs = allObjs;
SegInput.ConnectedComponents = ConnectedComponents;
SegInput.numObj = numObj;
SegInput.voxscale = voxscale;
SegInput.s = sizeXY;
SegInput.zs = sizeZ;
app = CellSizeCutoffGUI_App;
setup(app, SegInput);
waitfor(app, 'Output');
delete(app);
clear SegInput allObjs; % Free memory

%% Segmentation Loop with Bounding Boxes
col=1;
progbar = waitbar(0,'Segmenting... ');
fprintf('Starting segmentation with bounding box optimization...\n');

for i = 1:numObj
    waitbar(i/numObj, progbar);
    
    if numel(ConnectedComponents.PixelIdxList{1,i}) > CellSizeCutoff
        % use bounding box instead of full image, saves memory/time
        [idxX, idxY, idxZ] = ind2sub([sizeXY(1), sizeXY(2), sizeZ], ConnectedComponents.PixelIdxList{1,i});
        
        minX = max(min(idxX)-5, 1);
        maxX = min(max(idxX)+5, sizeXY(1));
        minY = max(min(idxY)-5, 1);
        maxY = min(max(idxY)+5, sizeXY(2));
        minZ = max(min(idxZ)-5, 1);
        maxZ = min(max(idxZ)+5, sizeZ);
        
        bbSize = [maxX-minX+1, maxY-minY+1, maxZ-minZ+1];
        
        % create small bounding box array instead of full size
        ex = false(bbSize);
        locX = idxX - minX + 1;
        locY = idxY - minY + 1;
        locZ = idxZ - minZ + 1;
        linLoc = sub2ind(bbSize, locX, locY, locZ);
        ex(linLoc) = true;
        
        se=strel('diamond',6);
        nucmask=imerode(ex,se);
        nucsize = round((CellSizeCutoff/50),0);
        nucmask=bwareaopen(nucmask,nucsize);
        indnuc=bwconncomp(nucmask);
        nuc = numel(indnuc.PixelIdxList);
        clear nucmask indnuc; 
        
        if nuc ==0
            nuc = 2;
            warning('Warning: The program finds 0 nuclei to segment your object into. Adjust se=strel(diamond,6) to a lower number.');
        end  
        if nuc ==1
            nuc = 2;
            warning('Warning: The program finds 1 nuclei to segment your object into. Adjust se=strel(diamond,6) to a lower number.');
        end
        
        [x,y,z]=ind2sub(bbSize,find(ex));
        points = [x y z];
        
        GMModel = fitgmdist(points,nuc,'replicates',3); 
        idx = cluster(GMModel,points);
        
        for j=1:nuc
            obj = (idx == j);
            object = points(obj,1:3);
            
            % Create small bounding box for this segment
            ex_seg = false(bbSize);
            locIdx = sub2ind(bbSize, object(:,1), object(:,2), object(:,3));
            ex_seg(locIdx) = true;
            
            individual=bwconncomp(ex_seg,26);
            NumberOfIdentifiedObjects = length(individual.PixelIdxList);
            
            for m = 1:NumberOfIdentifiedObjects
                % Convert local indices back to global
                [lx, ly, lz] = ind2sub(bbSize, individual.PixelIdxList{1,m});
                gx = lx + minX - 1;
                gy = ly + minY - 1;
                gz = lz + minZ - 1;
                globalIdx = sub2ind([sizeXY(1), sizeXY(2), sizeZ], gx, gy, gz);
                
                Microglia{1,col} = globalIdx;
                col=col+1;
            end
        end
        
        clear object individual x y z points ex ex_seg idxX idxY idxZ locIdx;
    else
        Microglia{1,col}=ConnectedComponents.PixelIdxList{1,i};
        col=col+1;
    end
end

if isgraphics(progbar)
    close(progbar);
end
fprintf('Segmentation complete.\n');

numObjSep = numel(Microglia);

for i = 1:numObjSep
    SepObjectList(i,1) = length(Microglia{1,i}); 
    SepObjectList(i,2) = i;  
end
SepObjectList = sortrows(SepObjectList,-1);
udSepObjectList = flipud(SepObjectList);

clear ConnectedComponents; % Free memory

% post-segmentation projection
AllSeparatedObjs = zeros(sizeXY(1), sizeXY(2), numObjSep, 'uint16');
fprintf('Creating separated object projections...\n');

for i = 1:numObjSep
    % create 2D projection where color indicates z-depth
    [x, y, ~] = ind2sub([sizeXY(1), sizeXY(2), sizeZ], Microglia{1,i});
    idx2D = sub2ind([sizeXY(1), sizeXY(2)], x, y);
    counts = accumarray(idx2D, 1, [sizeXY(1)*sizeXY(2), 1]);
    AllSeparatedObjs(:,:,i) = reshape(counts, sizeXY(1), sizeXY(2));
    clear idx2D counts;
end

fprintf('Projections complete.\n');

%% Territorial volume 
%Uses convhulln to create a 3D polygon around the object's external points.
%For the total occupied vs unoccupied volume, don't want to exclude any
%cells/processes. Use Microglia list here, not FullMg.

ConvexVol = zeros(numObjSep,1);
progbar = waitbar(0,'Finding territorial volume...');

for i = 1:numObjSep
    waitbar(i/numObjSep, progbar);
    
    % Get voxel coordinates for this object
    [x, y, z] = ind2sub(sz, [Microglia{1,i}]);
    obj = [y, x, z];
    obj = unique(obj, 'rows');
    
    % Check for degenerate cases
    if size(obj,1) < 4
        warning('Skipping object %d: only %d unique points.', i, size(obj,1));
        ConvexVol(i) = NaN;
        continue;
    end

    try
        [~, v] = convhulln(obj);
        ConvexVol(i) = v * voxscale;
    catch ME
        warning('convhulln failed for object %d: %s', i, ME.message);
        ConvexVol(i) = NaN;
    end
end

if isgraphics(progbar)
    close(progbar);
end

TotMgVol = nansum(ConvexVol);  % Calculate total volume of image covered by microglia (ignores NaNs)
ImageVol  = (sizeXY(1)*sizeXY(2)*sizeZ) * voxscale; %volume of image cube in um^3. 
EmptyVol = ImageVol - TotMgVol; %And the remaining 'empty space'.
PercentMgVol = (TotMgVol / ImageVol) * 100;

ImageDataTable.TotMgVol = TotMgVol(:);
ImageDataTable.ImageVol = ImageVol(:);
ImageDataTable.EmptyVol = EmptyVol(:);
ImageDataTable.PercentMgVol = PercentMgVol(:);


%% Full cells
FullCellInput.numObjSep = numObjSep;
FullCellInput.Microglia = Microglia;
FullCellInput.s = sizeXY;
FullCellInput.zs = sizeZ;
FullCellInput.voxscale = voxscale;
FullCellInput.SepObjectList = SepObjectList;
FullCellInput.udSepObjectList = udSepObjectList;
FullCellInput.AllSeparatedObjs = AllSeparatedObjs;
app = FullCellsGUI_App;
setup(app, FullCellInput);
waitfor(app, 'Output');
delete(app);
clear FullCellInput AllSeparatedObjs; 

% full cells processing
if KeepAllCells == 1
    FullMg = Microglia;
else   
    col=1;
    progbar = waitbar(0,'Finding All Full Cells...');
    
    for i = 1:numObjSep
        waitbar(i/numObjSep, progbar);
        if numel(Microglia{1,i}) >= SmCellCutoff
            if RemoveXY == 1
                % use bounding box
                [idxX, idxY, idxZ] = ind2sub([sizeXY(1), sizeXY(2), sizeZ], Microglia{1,i});
                
                minX = max(min(idxX)-1, 1);
                maxX = min(max(idxX)+1, sizeXY(1));
                minY = max(min(idxY)-1, 1);
                maxY = min(max(idxY)+1, sizeXY(2));
                minZ = max(min(idxZ)-1, 1);
                maxZ = min(max(idxZ)+1, sizeZ);
                
                bbSize = [maxX-minX+1, maxY-minY+1, maxZ-minZ+1];
                
                ex = false(bbSize);
                locX = idxX - minX + 1;
                locY = idxY - minY + 1;
                locZ = idxZ - minZ + 1;
                linLoc = sub2ind(bbSize, locX, locY, locZ);
                ex(linLoc) = true;
                
                antiborder = logical(padarray(ex,[0 0 1],0));
                cleared = bwclearborder(antiborder,26);
                nonedge = max(cleared(:));
                clear ex antiborder cleared linLoc locX locY locZ idxX idxY idxZ;
                
                if nonedge == 1
                   FullMg{1,col} = Microglia{1,i};
                   col = col+1;
                end
            else
                FullMg{1,col} = Microglia{1,i};
                col = col+1;   
            end  
        end
    end
    
    if isgraphics(progbar)
        close(progbar);
    end
end

numObjMg = numel(FullMg);

for i = 1:numObjMg
    MgObjectList(i,1) = length(FullMg{1,i}); 
    MgObjectList(i,2) = i;  
end
MgObjectList = sortrows(MgObjectList,-1);
udMgObjectList = flipud(MgObjectList);
num = 1:3:(3*numObjMg+1);
cmap = jet(max(num));
cmap(1,:) = zeros(1,3);

clear Microglia; % Free memory


%% Quantify 3D masks of cells
% Cell Volumes
NumberOfPixelsPerCell = cellfun(@numel,FullMg);
CellVolume = (NumberOfPixelsPerCell*voxscale)';
CellNum = (1:numel(FullMg))';
CellDataTable = table(CellNum, CellVolume, ... 
    'VariableNames', {'CellNum', 'Volume'});

[biggest,idx] = max(NumberOfPixelsPerCell);
MaxCellVol = biggest*voxscale; 

ImageDataTable.MaxCellVol = MaxCellVol(:);

%% Centroid Calculation
numCells = numel(FullMg);
centroids = zeros(numCells, 3);

fprintf('Calculating centroids...\n');
progbar = waitbar(0, 'Calculating centroids...');

for i = 1:numCells
    waitbar(i / numCells, progbar);
    % Calculate centroid directly from indices, no 3D array needed
    [x, y, z] = ind2sub([sizeXY(1), sizeXY(2), sizeZ], FullMg{1,i});
    centroids(i, 1) = mean(x);
    centroids(i, 2) = mean(y);
    centroids(i, 3) = mean(z);
end   

if isgraphics(progbar)
     close(progbar);
end

fprintf('Centroids complete.\n');

centum = centroids .* [pxXY, pxXY, pxZ];
centdist = pdist2(centum, centum); %Calculate distance from each centroid to all other centroids
centdist = nonzeros(centdist); %Remove all 0s (distance from one centroid to itself)
AvgDist = mean(centdist);
GliaDensity = numCells/ImageVol;

ImageDataTable.AvgDist = AvgDist(:);
ImageDataTable.GliaDensity = GliaDensity(:);

Q.SaveCells = {'How would you like to see your numbered cells?','Warning: A 3D image will take longer to process!'};
Q.SaveCellsChoice = questdlg(Q.SaveCells,'Save Image of Cells','3D figure', '2D image', 'Both', '2D image');
switch Q.SaveCellsChoice
    case '3D figure'
        Q.SaveCellsAns = 1;
    case '2D image'
        Q.SaveCellsAns = 2;    
    case 'Both'
        Q.SaveCellsAns = 3;
end

%% 2D Cell Rendering

if Q.SaveCellsAns == 2 || Q.SaveCellsAns == 3
    fullimg = ones(sizeXY(1), sizeXY(2));
    progbar = waitbar(0, 'Rendering cells in 2D (optimized)');
    for i = 1:numObjMg
        waitbar(i / numObjMg, progbar);
        j = udMgObjectList(i, 2);
        [x, y, ~] = ind2sub([sizeXY(1), sizeXY(2), sizeZ], FullMg{1,j});
        idx2D = sub2ind([sizeXY(1), sizeXY(2)], x, y);
        counts = accumarray(idx2D, 1, [sizeXY(1)*sizeXY(2), 1]);
        flatex = reshape(counts, sizeXY(1), sizeXY(2));
        OutlineImage = zeros(sizeXY(1), sizeXY(2));
        OutlineImage(flatex > 1) = 1;
        se = strel('diamond', 4);
        Outline = imdilate(OutlineImage, se);
        fullimg(Outline == 1) = 1;
        fullimg(flatex > 1) = num(1, i + 1);
        clear idx2D counts flatex OutlineImage Outline;
    end
    if isgraphics(progbar)
        close(progbar);
    end
    
    figName = [file, '_step2__Glia_Cells_2D'];
    fig = figure('Name', figName);
    imagesc(fullimg);
    colormap(cmap);
    colorbar('Ticks', [1, max(num)], 'TickLabels', {'Small', 'Large'});
    daspect([1 1 1]);
    hold on;
    centroidsXY = centroids(:, 1:2);
    plot(centroidsXY(:,2), centroidsXY(:,1), ...
        'w.', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    offsetX = 25;
    for i = 1:numObjMg
        xy = centroidsXY(i, :);
        if all(~isnan(xy))
            text(xy(2) + offsetX, xy(1), num2str(CellNum(i)), ...
                'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'center');
        end
    end
    hold off;
    filename = ([name '_step2_Glia_Cells_2D.png']);
    exportgraphics(fig, fullfile(SaveFolderPath, filename), 'Resolution', 300);
end
clear fullimg;

if Q.SaveCellsAns == 1 || Q.SaveCellsAns == 3
    % Note: 3D rendering still uses full volume but this is necessary for isosurface
    fullvol = zeros(sizeXY(1), sizeXY(2), sizeZ, 'uint16');
    fprintf('Building 3D volume for rendering\n');
    for i = 1:numObjMg
        j = udMgObjectList(i, 2);
        fullvol(FullMg{1, j}) = num(1, i + 1);
    end
    figName = [file, '_step2_Glia_Cells_3D'];
    fig = figure('Name', figName);
    hold on;
    centroidsXY = centroids(:, 1:2);
    progbar = waitbar(0, 'Rendering cells in 3D (this may take a while)');
    for i = 1:numObjMg
        waitbar(i / numObjMg, progbar);
        objMask = (fullvol == num(1, i + 1));
        if any(objMask(:))
            fv = isosurface(objMask, 0.5);
            if ~isempty(fv.vertices)
                patch(fv, ...
                    'FaceColor', cmap(num(1, i + 1), :), ...
                    'FaceAlpha', 0.9, ...
                    'EdgeColor', 'none');
            end
        end
    end
    if isgraphics(progbar)
        close(progbar);
    end
    daspect([1 1 1]);
    axis equal tight off;
    view(0, -90);
    camlight('headlight');
    lighting gouraud;
    colormap(fig, cmap);
    colorbar('Ticks', [min(num), max(num)], ...
        'TickLabels', {'Small', 'Large'}, ...
        'FontSize', 10);
    clim([min(num), max(num)]);
    scatter3(centroids(:,2), centroids(:,1), centroids(:,3), ...
        25, 'w', 'filled');
    offsetX = 25;
    for i = 1:numObjMg
        xy = centroidsXY(i, :);
        if all(~isnan(xy))
            text(xy(2) + offsetX, xy(1), num2str(CellNum(i)), ...
                'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'center');
        end
    end
    figFilename = fullfile(SaveFolderPath, [name '_step2_Glia_Cells_3D.fig']);
    savefig(fig, figFilename);
    pngFilename = fullfile(SaveFolderPath, [name '_step2_Glia_Cells_3D.png']);
    exportgraphics(fig, pngFilename, 'Resolution', 300);
    hold off;
    clear fullvol fv objMask;
end

%% Convex territorial volume of full cells. 
Q.ConvexImg = {'Would you like to see a convex volume image of each full cell?'};
Q.ConvexImgChoice = questdlg(Q.ConvexImg,'Show Convex Volume Images?','Yes please!', 'No thanks','No thanks');
switch Q.ConvexImgChoice 
    case 'Yes please!'
        Q.ConvexImgAns = 1;
    case 'No thanks'
        Q.ConvexImgAns = 0;    
end

FullCellTerritoryVol = zeros(numObjMg,1);
for i = 1:numObjMg
    [x,y,z] = ind2sub(sz,[FullMg{1,i}]); %input: size of array ind values come from, list of values to convert
    obj = [y,x,z]; %concatenate x y z coordinates
    
    try
        [j, v] = convhulln(obj);
        FullCellTerritoryVol(i,:) = v * voxscale;
    
        if Q.ConvexImgAns == 1
            figure;
            trisurf(j, obj(:,1), obj(:,2), obj(:,3));
            axis([0 sizeXY(1) 0 sizeXY(2) 0 sizeZ]);
            daspect([1 1 1]);
        end
    
    catch ME
        warning('Could not compute convex hull for object %d: %s', i, ME.message);
        FullCellTerritoryVol(i,:) = NaN;
    end
    clear obj x y z;
end
CellDataTable.TerritoryVolume = FullCellTerritoryVol(:);

%% Colocalization

Q.Coloc = {'Would you like to colocalize your cells with any other channel(s) in your image?'};
Q.ColocChoice = questdlg(Q.Coloc,'Colocalization?','Yes!', 'No thanks, just cell structure data', 'Yes!');
    switch Q.ColocChoice
        case 'Yes!'
            Q.ColocAns = 1;
        case 'No thanks, just cell structure data'
            Q.ColocAns = 2;    
    end

%% Select Colocalization Channel(s) of Interest
if (Q.ColocAns == 1)
    app = ColocInfo_App;
    waitfor(app, 'Output');
    delete(app);

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
    
    % Normalize for Coloc
    fields = fieldnames(ColocImgs);
    ColocImgsNorm = struct();
    
    for f = 1:numel(fields)
        fieldName = fields{f};
        ColocImgsNorm.(fieldName) = mat2gray(ColocImgs.(fieldName));
    end

elseif (Q.ColocAns == 2)
    % create dummy variables to save in parameters file
    ColocChannelNames = []; ColocChannels = []; ColocShowImg = []; ColocSaveOrig = [];
    ColocThresholdValues = []; ColocSizeCutoffs = []; ColocThreshMode = []; ColocThreshMethod = [];
    ColocFillHoles = []; BatchColocExclude = []; BatchColocInclude = []; BatchColocComp = [];
    BatchColocObjectsSelect = []; BatchColocOrigImgs = []; BatchColocOverlapImgs = [];
    BatchProxCC = []; ProxSkip = []; ColocSaveModeOverlay = []; ColocSaveModeOrig = [];
end
clear imgChannels ColocImgs;

%% Coloc Threshold
if (Q.ColocAns == 1)
    ColocThreshInput.length = sizeZ;
    ColocThreshInput.size = sizeXY;
    ColocThreshInput.voxscale = voxscale;
    ColocThreshInput.ColocImgs = ColocImgsNorm;
    app = ColocThresh_App;
    setup(app, ColocThreshInput);
    waitfor(app, 'Output');
    delete(app);
    clear ColocThreshInput;
end

%% Save Images of Other Channels
if Q.ColocAns == 1

    if  ColocSaveOrig == 1
        fields = fieldnames(ColocImgsNorm);
        for i = 1:numel(fields)
            chName = fields{i};
            img3D = ColocImgsNorm.(chName);
            img2D = max(img3D, [], 3);
            figName = [name, '_step3_', chName, '_OrigMAX'];
            fig = figure('Name', figName);
            imagesc(img2D);
            colormap gray;
            daspect([1 1 1]);
            filename = fullfile(SaveFolderPath, [name, '_step3_', chName, '_OrigMAX.png']);
            exportgraphics(fig, filename, 'Resolution', 300);
            clear img3D img2D; 
        end
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
            fig = figure('Name', ['_step4_', chName , '_Thresh_']);
            imagesc(img2D);
            colormap hot;
            daspect([1 1 1]);
            filename = fullfile(SaveFolderPath, [name, '_step4_', chName, '_Thresh.png']);
            exportgraphics(fig, filename, 'Resolution', 300);
            clear img3D img2D depthMap; 
        end
    end
end
clear ColocImgsNorm; 

%% Explanation of colocalization steps
if (Q.ColocAns == 1)
    app = ColocInstructions_App;
    setup(app);
    waitfor(app, 'Output');
    delete(app);
end

%% Decide how to handle multi-colocalization
if (Q.ColocAns == 1)
    ColocFilterExclInput.chans = ColocSizeFilteredImgs;
    app = ColocFilterExclusive_App;
    setup(app, ColocFilterExclInput);
    waitfor(app, 'Output');
    delete(app);
    clear ColocFilterExclInput;

    % Filter by exclusion criteria
    if ~exist('ColocExclude', 'var')
        BatchColocExclude = [];
    else
        filterNames = fieldnames(ColocExclude);
        for i = 1:numel(filterNames) % loop over named exclusion filtered objects
            fname = filterNames{i};
            % Grab image from each field
            bwKeep = ColocExclude.(fname).ExcludeChan;
            bwIf = ColocExclude.(fname).IfChanExcl;
            mode = ColocExclude.(fname).ExclMode;
        
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
            clear bwIf bwKeep maskKept_obj maskKept_vox props labeledKeep;
        end
    end
    clear ColocExclude;

    % User define inclusion filtering criteria
    ColocFilterInclInput.chans = ColocSizeFilteredImgs;
    app = ColocFilterInclusive_App;
    setup(app, ColocFilterInclInput);
    waitfor(app, 'Output');
    delete(app);
    clear ColocFilterInclInput;

    % Filter by inclusion criteria
    if ~exist('ColocInclude', 'var')
        BatchColocInclude = [];
    else 
        filterNames = fieldnames(ColocInclude);
        for i = 1:numel(filterNames)
            fname = filterNames{i};
            bwKeep = ColocInclude.(fname).KeepChan;
            bwIf = ColocInclude.(fname).IfChan;
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
    clear ColocInclude;

    % Multi-Channel Composite Objects
    ColocCompInput.chans = ColocSizeFilteredImgs;
    app = ColocCompObjects_App;
    setup(app, ColocCompInput);
    waitfor(app, 'Output');
    delete(app);
    clear ColocCompInput;

    if ~exist('ColocComp', 'var')
        BatchColocComp = [];
    else
        ColocCompImgs = struct(); % create empty struct to catch multi-chan composite object images created below
        subNames = fieldnames(ColocComp); % extract names of user-defined multi-chan composites
        
        progbar = waitbar(0,'Filtering and Combining...');
        for i = 1:numel(subNames) % loop over user-defined multi-chan composites
            waitbar(i / numel(subNames), progbar);
            subName = subNames{i};
            subStruct = ColocComp.(subName); % extract user-defined criteria for each multi-chan composite as separate struct
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

        if isgraphics(progbar)
            close(progbar);
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
    clear ColocComp;
end

%% Find objects in other channels and composites
if (Q.ColocAns == 1)
    ColocConnectedComponents = struct();
    channels = fieldnames(ColocSizeFilteredImgs);
    progbar = waitbar(0,'Finding Connected Components...');
    for i = 1:numel(channels)
        waitbar(i / numel(channels), progbar);
        chan = channels{i};
        fprintf('Finding connected components in: %s\n', chan);
        binImg = ColocSizeFilteredImgs.(chan);
        CC = bwconncomp(binImg, 26);
        ColocConnectedComponents.(chan) = CC;
        clear binImg 
    end
    if isgraphics(progbar)
        close(progbar);
    end

    ColocObjects = struct();   
    for i = 1:numel(channels)
        chan = channels{i};
        ColocObjects.(chan) = ColocConnectedComponents.(chan).PixelIdxList;
    end

    % Get image level-data about objects in other channels and add to ImageDataTable
    objCountFields = cellfun(@(f) ['ObjCount_', f], channels, 'UniformOutput', false);
    objDensityFields = cellfun(@(f) ['ObjDensity_', f], channels, 'UniformOutput', false);
    objCounts = cellfun(@(f) numel(ColocObjects.(f)), channels)';
    objDensities = objCounts / ImageDataTable.ImageVol;
    allFields = [objCountFields(:); objDensityFields(:)]';
    ObjCountTable = array2table([objCounts, objDensities], 'VariableNames', allFields);
    ImageDataTable = [ImageDataTable, ObjCountTable];

end

%% Colocalize objects with glia structure
if (Q.ColocAns == 1)
    ColocSelectInput.chans = ColocObjects;
    app = ColocGliaSelect_App;
    setup(app, ColocSelectInput);
    waitfor(app, 'Output');
    delete(app);
    clear ColocSelectInput;
end

if (Q.ColocAns == 1)
    clear ColocObjects;

    ColocResults = table(CellDataTable.CellNum, 'VariableNames', {'CellNum'});
    channels = fieldnames(ColocObjectsSelect);
    numChannels = numel(channels);
    ChannelMasks = struct();
    ColocOverlapIndices = struct();

    for i = 1:numChannels
        chan = channels{i};
        mask = false(size(ColocSizeFilteredImgs.(chan)));
        for j = 1:numel(ColocObjectsSelect.(chan))
            mask(ColocObjectsSelect.(chan){j}) = true;
        end
        ChannelMasks.(chan) = mask;
        clear mask;
    end

    % Single Channel Ovelaps
    progbar = waitbar(0,'Finding single channel overlaps with glia structures...');
    for i = 1:numChannels
        waitbar(i / numChannels, progbar);
        chan = channels{i};
        mask = ChannelMasks.(chan);
        channelObjects = ColocObjectsSelect.(chan);
    
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
            clear overlappingObjects overlapping objVoxels connectedVoxels; % Free memory
        end
    
        % Store overlap summary stats in table
        ColocResults.([chan '_ColocVolume_um3']) = OverlapVolumes;
        ColocResults.([chan '_ColocNumObjects']) = OverlapCounts;
    
        % Store overlap voxel indices
        ColocOverlapIndices.(chan) = unique(OverlapIndices);
        clear OverlapIndices mask;
    end

    if isgraphics(progbar)
        close(progbar);
    end
    
  % Multi-Channel Overlaps
    if Q.ColocAllAns == 1
        progbar = waitbar(0, 'Finding multi-channel overlaps...');
        numCombos = 0;
        
        for k = 2:numChannels
            combos = nchoosek(1:numChannels, k);
            numCombos = numCombos + size(combos, 1);
        end
        
        comboCount = 0;
        for k = 2:numChannels
            combos = nchoosek(1:numChannels, k);
            
            for comboIdx = 1:size(combos, 1)
                comboCount = comboCount + 1;
                waitbar(comboCount / numCombos, progbar);
                
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
    
        if isgraphics(progbar)
            close(progbar);
        end
    end
end

if exist('ColocResults', 'var') && istable(ColocResults)
    CellDataTable = join(CellDataTable, ColocResults, 'Keys', 'CellNum', 'KeepOneCopy', 'CellNum');
end

%% Visualize Colocalization
if (Q.ColocAns == 1)
    ColocSaveInput.OrigChans = ColocObjectsSelect;
    ColocSaveInput.OverlapChans = ColocOverlapIndices;
    app = ColocSaveImages_App;
    setup(app, ColocSaveInput);
    waitfor(app, 'Output');
    delete(app);
    clear ColocSaveInput;
end

if (Q.ColocAns == 1)
    if (ColocSaveModeOverlay ~= 0)
        MergedMg = zeros(sz);
        for i = 1:length(FullMg)
            MergedMg(FullMg{1, i}) = 1;
        end
        
        channels = fieldnames(ColocOrigImgs);
        MergedOrigImgs = struct();
        for f = 1:numel(channels)
            field = channels{f};
            merged = zeros(sz);
            for i = 1:numel(ColocOrigImgs)
                indices = ColocOrigImgs(i).(field);
                for c = 1:numel(indices)
                    inds = indices{c};
                    merged(inds) = 1;
                end
                clear indices inds; % Free memory
            end
            MergedOrigImgs.(field) = merged;
            clear merged; % Free memory
        end
    end
end

%% 2D All Objects From Other Channels
if (Q.ColocAns == 1)
    if (ColocSaveModeOrig == 2) || (ColocSaveModeOrig == 1)
        channels = fieldnames(MergedOrigImgs);
        progbar = waitbar(0, 'Rendering in 2D...');
        for i = 1:numel(channels)
            waitbar(i/numel(channels), progbar);
            field = channels{i};
            vol = MergedOrigImgs.(field);
            mipZ = max(vol, [], 3);
            mipZ = cat(3, mipZ, zeros(size(mipZ)), zeros(size(mipZ)));
            fig = figure('Name', field);
            imagesc(mipZ);
            daspect([1 1 1]);
            filename = fullfile(SaveFolderPath, [name '_step5_' field '_Obj_2D.png']);
            exportgraphics(fig, filename, 'Resolution', 300);
            clear mipZ vol;
        end
        if isgraphics(progbar)
            close(progbar);
        end
    end
end

%% 3D All Objects From Other Channels (long processing time)
if (Q.ColocAns == 1)
    if (ColocSaveModeOrig == 3) || (ColocSaveModeOrig == 1)
        channels = fieldnames(MergedOrigImgs);
        progbar = waitbar(0, 'Rendering in 3D (this may take a while)');
        for i = 1:numel(channels)
            waitbar(i/numel(channels), progbar);
            field = channels{i};
            vol = MergedOrigImgs.(field);
            fv = isosurface(vol, 0.5);
        
            fig = figure('Name', ['_3D_ ' field]);
            hold on;
            patch(fv, ...
                'FaceColor', [1, 0, 0], ...
                'FaceAlpha', 1, ...
                'EdgeColor', 'none');
        
            camlight; lighting gouraud;
            view(0, -90);
            axis equal tight off;
            daspect([1 1 1]);   
            savefig(fig, fullfile(SaveFolderPath, [name '_step5_' field '_Obj_3D.fig']));
            filename = fullfile(SaveFolderPath, [name '_step5_' field '_Obj_3D.png']);
            exportgraphics(fig, filename, 'Resolution', 300);
            
            hold off;
            set(fig, 'Visible', 'on');
            clear vol fv; % Free memory
        end
        if isgraphics(progbar)
            close(progbar);
        end
    end
end

%% 2D Overlapping Objects Only
if (Q.ColocAns == 1)
    if (ColocSaveModeOverlay == 2) || (ColocSaveModeOverlay == 1)
        MgZ = double(max(MergedMg, [], 3) > 0);
        channels = fieldnames(ColocOverlapImgs);
        for i = 1:numel(channels)
            chan = channels{i};
            inds = ColocOverlapImgs.(chan);

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
            fig = figure('Name', ['step6_' chan '_Overlap_2D']);
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

%% 3D Overlapping Objects Only
if (Q.ColocAns == 1)
    if (ColocSaveModeOverlay == 3) || (ColocSaveModeOverlay == 1)
        channels = fieldnames(ColocOverlapImgs);
        MgMask = logical(MergedMg);
        progbar = waitbar(0, 'Rendering overlaps in 3D (this may take a while)');
        for i = 1:numel(channels)
            waitbar(i/numel(channels), progbar);
            ch = channels{i};
            inds = ColocOverlapImgs.(ch);   
           
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
            disp(['Getting overlaps for ' ch '...']);
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
            fig = figure('Visible', 'off', 'Name', ['_Overlap_3D_' ch]);
            hold on;
        
            % Render other channel objects in red
            disp(['Rendering ' ch ' (Red)']);
            if ~isempty(fvMerged.vertices)
                patch(fvMerged, ...
                    'FaceColor', [1, 0, 0], ...
                    'FaceAlpha', 0.1, ...
                    'EdgeColor', 'none');
            end
        
            % Render microglia in green
            disp('Rendering MergedMg (Green)');
            if ~isempty(fvMergedMg.vertices)
                patch(fvMergedMg, ...
                    'FaceColor', [0, 1, 0], ...
                    'FaceAlpha', 0.3, ...
                    'EdgeColor', 'none');
            end
        
            % Render overlap in yellow
            disp('Rendering Overlap (Yellow)');
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
            
            hold off;
            set(fig, 'Visible', 'on');
            clear inds merged overlap onlyMerged onlyMergedMg fvMerged fvMergedMg fvOverlap;
            
        end
        if isgraphics(progbar)
            close(progbar);
        end
    end
end

% Clean up workspace before skeletonization step
clear allObjs binImg colocObjs CC ChannelMasks ColocConnectedComponents ...
    ColocImgs ColocComp ColocCompImgs ColocImgsNorm colocIntersect colocUnion ColocObjects ColocObjectsSelect ColocOrigImgs ...
    ColocOverlapImgs ColocOverlapIndices ColocSizeFilteredImgs channelsStruct depthMap DetectedObjs ex ...
    flatex fullimg imgChannels img imgGlia imgGliaNorm img2D img3D mask merged MergedMg MergedOrigImgs MgZ ... 
    Microglia mip NoiseIm onlyMerged onlyMergedMg overlap overlapping ...
    overlappingObjects overlayRGB OverlapIndices ProxCC subStruct thisChannel x y;

%% Continue to Skeletonization?
Q.SkelContinue = ['The next step is skeletonization of all identified cells, which can be a computationally intensive and time consuming process.' ...
    'If you are not interested in this data, you can safely skip this step.'];
Q.SkelContinueChoice = questdlg(Q.SkelContinue,'Continue to Skeletonization?', ...
    'Yes, continue to skeletonization step', 'No thanks, just output current data file', 'Not right now, but yes for batch process parameters', ...
    'Yes, continue to skeletonization step');
switch Q.SkelContinueChoice 
    case 'Yes, continue to skeletonization step'
        Q.SkelContinueAns = 1;
    case 'No thanks, just output current data file'
        Q.SkelContinueAns = 2;
    case 'Not right now, but yes for batch process parameters'
        Q.SkelContinueAns = 3;
end

if Q.SkelContinueAns == 1 || Q.SkelContinueAns == 3
    Q.SkelMode = ['Would you like to include all small processes or fillipodia in your skeleton analysis? ','Note: only keep small processes if you want all fine fillipodia - this will increase processing time.'];
    Q.SkelModeChoice = questdlg(Q.SkelMode,'Skeletonization Method','Keep small processes', 'Only major branches', 'Only major branches');
    switch Q.SkelModeChoice
        case 'Keep small processes'
            Q.SkelModeAns = 1;
        case 'Only major branches'
            Q.SkelModeAns = 2;    
    end
end

%% 3D Skeleton
% Instead of using the entire image, just create a bounding box (with some padding) around each glia from FullMg, then trace those. 
% This reduces memory load substantially.
% If the program is unable to properly identify a centroid or endpoints, it will output a 0 and move to the next cell.

if Q.SkelContinueAns == 1

    kernel(:,:,1) = [1 1 1; 1 1 1; 1 1 1];
    kernel(:,:,2) = [1 1 1; 1 0 1; 1 1 1];
    kernel(:,:,3) = [1 1 1; 1 1 1; 1 1 1]; 
    
    numendpts = zeros(numel(FullMg),1);
    numbranchpts = zeros(numel(FullMg),1);
    MaxBranchLength = zeros(numel(FullMg),1);
    MinBranchLength = zeros(numel(FullMg),1);
    AvgBranchLength = zeros(numel(FullMg),1);
    
    BranchLengthList=cell(1,numel(FullMg));
    
    %% Trace skeletons
    
    % Preallocate to save info from parfor loop
    SkelErrorLog = cell(1, numel(FullMg)); % Error log for debugging later
    Skels(1, numel(FullMg)) = struct('ID', [], 'Mask', [], 'Status', '');
    adjust_pxXY_list = zeros(1, numel(FullMg));
    
    tic; % timer
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
       
       disp(['completed cell ' num2str(i) ' of ' num2str(numel(FullMg))]); %#ok<PFBNS> %To see which cell we are currently prcoessing.
    end
    
    delete(gcp('nocreate')); %close parallel pool so error isn't generated when program is run again.
    
    totalTime = toc;
    fprintf('parfor completed in %.2f seconds (%.2f minutes)\n', totalTime, totalTime/60);
    
    % Skeleton Tracing Error Log
    for i = 1:numel(SkelErrorLog)
        if ~isempty(SkelErrorLog{i})
            fprintf('Error in cell %d:\n', i);
            disp(SkelErrorLog{i}.getReport('extended', 'hyperlinks', 'off'));
            fprintf('--------------------------------------------------\n');
        end
    end
    
    %% Calculate total branch length
    
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
    
    %% Save images of skeletons
    % this step has been moved outside parfor loop to make customization easier
    progbar = waitbar(0,'Saving Images of Each Skeleton');
    for i = 1:length(Skels)
        waitbar (i/numel(FullMg), progbar);
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
    
        figName = [name, '_Cell', num2str(i)];
        figure('Name', figName, 'Visible', 'off');
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
    if  isgraphics(progbar)
        close(progbar);
    end
    
    % Output Results
    CellDataTable.AvgBranchLength_um = AvgBranchLength(:);
    CellDataTable.NumEndpoints = numendpts(:);
    CellDataTable.NumBranchpoints = numbranchpts(:);
    CellDataTable.TotalBranchLength_um = TotalBranchLengths;
    
end

%% Final Data Output
nCells = height(CellDataTable); 
CellDataTable.FileName = repmat({name}, nCells, 1); % add a column with filename (useful for compiling multiple .csv files from batch processing)
cols = CellDataTable.Properties.VariableNames;
CellDataTable = CellDataTable(:, [{'FileName'}, setdiff(cols, {'FileName'}, 'stable')]); % move filename to first column
csvFilename = fullfile(SaveFolderPath, ['CellData_', name, '.csv']);
writetable(CellDataTable, csvFilename);

nCells1 = height(ImageDataTable); 
ImageDataTable.FileName = repmat({name}, nCells1, 1); % add a column with filename (useful for compiling multiple .csv files from batch processing)
cols1 = ImageDataTable.Properties.VariableNames;
ImageDataTable = ImageDataTable(:, [{'FileName'}, setdiff(cols1, {'FileName'}, 'stable')]); % move filename to first column
csvFilename1 = fullfile(SaveFolderPath, ['ImageData_', name, '.csv']);
writetable(ImageDataTable, csvFilename1);

%% Save FullMg
% this is only relevant if you plan to use Script3
FullMgPath = fullfile(SaveFolderPath, ['FullMg_', name, '.mat']);
save(FullMgPath, 'FullMg')

fprintf('Microglia from this image saved as a .mat file called "FullMg_%s.mat" \n', name);
fprintf('You can find this file in %s \n', SaveFolderPath);

%% Save Parameters File
FullFilePath = fullfile(SaveFolderPath, ['Parameters_', name, '.mat']);
save(FullFilePath, 'Q', 'pxXY', 'pxZ', 'sizeXY', 'sizeZ', 'voxscale', 'GliaChannel', ...
    'SecondGliaChanYN', 'SecondGliaChanChoice', 'SaveOverlapImg', 'GliaChannel2', ...
    'noise', 'ThreshSet', 'ThreshMethod', 'ThreshMode', 'FillGaps', 'FillGapsRadius', ...
    'noise2', 'ThreshSet2', 'ThreshMethod2', 'ThreshMode2', 'FillGaps2', 'FillGapsRadius2', ...
    'CellSizeCutoff', 'SaveFolderPath', 'SmCellCutoff', 'RemoveXY', 'KeepAllCells', ...
    'ColocChannels', 'ColocChannelNames', 'ColocShowImg', 'ColocSaveOrig', ...
    'ColocThresholdValues', 'ColocSizeCutoffs', 'ColocThreshMode', 'ColocThreshMethod', 'ColocFillHoles', ...
    'BatchColocExclude', 'BatchColocInclude', 'BatchColocComp', ...
    'BatchColocObjectsSelect', 'BatchColocOrigImgs', 'BatchColocOverlapImgs', ...
    'ColocSaveModeOverlay', 'ColocSaveModeOrig');

fprintf('Parameters file saved to %s \n', SaveFolderPath);
fprintf('Parameters file name: Parameters_%s \n', name);
fprintf('Single image processing complete. You can now load the parameters file above in the Batch Process script. \n')