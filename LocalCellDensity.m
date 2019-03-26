
%% INPUTS NEEDED %%

% 1: Typical Diameter of a cell (Fancy name for how big the GaussianBlur will be), e.g.: TypicalCellDiameter = 150; %
% 2: Some Scaling Factor for cell type (Can be discarded as it just scales"1"), e.g.: ScalingFactorCellType = 2; %
% 3: Another Scaling Factor for "extended" locality, e.g.: ScalingFactorPara = 5; %
% 4: The size of the actual image, I calculate it from my segmentationimages which is here "DoubleNucleusImage" %
% 5: Centroids of Nucleus for which data should be extracted, e.g.:NucleusCentroidY and NucleusCentroidX %
% (6): In my case the labels for which data should be extracted but just to preallocate a matrix... %
%% OUTPUTS GIVEN %%

% Size calculated in my case from the labels in the image.. guess they need
% to be extracted... I do it like this.. with ObjectLabels being all
% Objects found in the label image (or MetaData in my case)

LocalCDCurrentCells = zeros(size(ObjectLabels,1),1); % LocalCellDensity
ParaCDCurrentCells = zeros(size(ObjectLabels,1),1); % Extended LocalCellDensity
LonerCurrentCells = zeros(size(ObjectLabels,1),1); % Whether a cell is single or not

%% Calculations (Of course I don't know how this is handled in the Iterator Pipeline)

% Local and Para Cell Density (Old Feature with Gaussian on Centroids) %

% Estimate Filter Params and generate gaussian for Local
LocalFilterSize = TypicalCellDiameter*ScalingFactorCellType;
LocalFilterSigma = LocalFilterSize/6; % Denominator has to be tested..
LocalPSF = fspecial('gaussian',LocalFilterSize,LocalFilterSigma);
% Estimate Filter Params and generate gaussian for Parakrine
ParaFilterSize = TypicalCellDiameter*ScalingFactorCellType*ScalingFactorPara;
ParaFilterSigma = ParaFilterSize/6;
ParaPSF = fspecial('gaussian',ParaFilterSize,ParaFilterSigma);

% Generate empty Dot Image
DotImage = zeros(size(DoubleNucleusImage));
% Generate linear Index of Nucleus locations
LinearIndex = sub2ind(size(DotImage),NucleusCentroidY,NucleusCentroidX);
% Put the Dots into Dots..
DotImage(LinearIndex) = 1;
% Blur Image using the filter..
LocalCDImage = imfilter(DotImage,LocalPSF,'symmetric','conv');
ParaCDImage = imfilter(DotImage,ParaPSF,'symmetric','conv');

% Calculate respective Densities for each cell...
LocalCDCurrentCells = LocalCDImage(LinearIndex);
ParaCDCurrentCells = ParaCDImage(LinearIndex);

% Index cells which are loner cells
MaximumLocalPSF = max(LocalPSF(:));
LonerCurrentCells = LocalCDCurrentCells <=(MaximumLocalPSF);
   
