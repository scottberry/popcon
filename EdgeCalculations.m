
%% INPUTS NEEDED %%

% 1: Typical Diameter of a cell (Fancy name for how big the GaussianBlur will be), e.g.: TypicalCellDiameter = 150; %
% 2: The size of the actual image, I calculate it from my segmentationimages which is here "DoubleNucleusImage" %
% 3: Centroids of Nucleus for which data should be extracted, e.g.:NucleusCentroidY and NucleusCentroidX %
% (4): In my case the labels for which data should be extracted but just to preallocate a matrix... %
%% OUTPUTS GIVEN %%

% Size calculated in my case from the labels in the image.. guess they need
% to be extracted... I do it like this.. with ObjectLabels being all
% Objects found in the label image (or MetaData in my case)

EdgePerCell = zeros(size(ObjectLabels,1),1); % Indicator whether a cell is at the edge or not
DistanceToEdge = zeros(size(ObjectLabels,1),1); % Well... Distance to Edge


%% Calculations (Of course I don't know how this is handled in the Iterator Pipeline)

% Do Edge Cell Calculations %

% Generate empty Dot Image
DotImage = zeros(size(DoubleNucleusImage));
% Generate linear Index of Nucleus locations
LinearIndex = sub2ind(size(DotImage),NucleusCentroidY,NucleusCentroidX);
% Put the Dots into Dots..
DotImage(LinearIndex) = 1;

% Calculate disk to stretch images
DilationFilter = strel('disk',round(TypicalCellDiameter./2));
% Dilate Image
DotImageEmptySpace = imdilate(DotImage,DilationFilter);

% Make Edge of Images non Edge
DotImageEmptySpace(:,1) = 1;
DotImageEmptySpace(:,end) = 1;
DotImageEmptySpace(1,:) = 1;
DotImageEmptySpace(end,:) = 1;

% Invert Image to mark empty spaces..
ImageEmptySpace = ~DotImageEmptySpace;


% Label all empty spaces and get area count
EmptySpaceLabel = bwlabel(ImageEmptySpace);
SpaceProbs = regionprops(EmptySpaceLabel,'Area');
EmptySpaceSize = cat(1,SpaceProbs.Area);

% Find to small empty regions and exclude those (Probably Optional.. or not)
TooSmallID = find(EmptySpaceSize < 2*((pi*TypicalCellDiameter^2)./2.5));
% Remove to small Areas
if ~isempty(TooSmallID)
    ImageEmptySpace(ismember(EmptySpaceLabel,TooSmallID)) = 0;
end

% Find Edges
EdgeImageEmptySpace = bwperim(ImageEmptySpace);
[EdgeEmptySpaceX,EdgeEmptySpaceY] = find(EdgeImageEmptySpace);

% Calculate Distance to closest Cell for each Edge..
ClosestCellPerEdgePixelID = NaN(size(EdgeEmptySpaceX,1),1);

for CurrentPixel = 1:size(EdgeEmptySpaceX,1)
    [foo,ClosestCellPerEdgePixelID(CurrentPixel)] = min(sqrt((NucleusCentroidY - EdgeEmptySpaceY(CurrentPixel)).^2 + (NucleusCentroidX - EdgeEmptySpaceX(CurrentPixel)).^2));
end

% Determine if a cell is an edge Cell
EdgePerCell = zeros(size(NucleusCentroidX,1),1);
EdgePerCell(unique(ClosestCellPerEdgePixelID)) = 1;

% Calculate closest Edge Pixel for every Cell
DistanceToEdgePerCell = NaN(size(NucleusCentroidX,1),1);

for CurrentCell = 1:size(NucleusCentroidX,1)
    DistanceToEdgePerCell(CurrentCell,1) = min(sqrt((EdgeEmptySpaceY - NucleusCentroidY(CurrentCell,1)).^2 + (EdgeEmptySpaceX - NucleusCentroidX(CurrentCell,1)).^2));
end

   
