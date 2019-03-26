
%% INPUTS NEEDED %%

% 1: The size of the actual image, I calculate it from my segmentationimages which is here "DoubleNucleusImage" %
% 2: Centroids of Nucleus for which data should be extracted, e.g.:NucleusCentroidY and NucleusCentroidX %
% (3): In my case the labels for which data should be extracted but just to preallocate a matrix... %
%% OUTPUTS GIVEN %%

% Size calculated in my case from the labels in the image.. guess they need
% to be extracted... I do it like this.. with ObjectLabels being all
% Objects found in the label image (or MetaData in my case)

LocalCCCurrentCells = zeros(size(ObjectLabels,1),1); % Local Cell Crowding


%% Calculations (Of course I don't know how this is handled in the Iterator Pipeline)

% Local Cell Crowding (From Gabriele) %

% Quantify the maximal distance in an Image
ImageY = size(DoubleNucleusImage,1);
ImageX = size(DoubleNucleusImage,2);

% Allocate Matrices for calculation
MaximalDistance = sqrt(ImageY^2+ImageX^2);

% Loop over each single cell and calculate scored Distance

SeparationRandom = NaN(size(NucleusCentroidY,1),1);
SeparationTrue = NaN(size(NucleusCentroidY,1),1);


for CurrentCell = 1:size(NucleusCentroidY)
    % Randomize position
    RandomX = randsample(NucleusCentroidX,size(NucleusCentroidX,1)-1,'true');
    RandomY = randsample(NucleusCentroidY,size(NucleusCentroidX,1)-1,'true');
    % Ensure that true cell spots is not included..
    while any((RandomX == NucleusCentroidX(CurrentCell)) & (RandomY == NucleusCentroidY(CurrentCell)))
        RandomX = randsample(NucleusCentroidX,size(NucleusCentroidX,1)-1,'true');
        RandomY = randsample(NucleusCentroidY,size(NucleusCentroidX,1)-1,'true');
    end
    
    % Random distance
    RandomDistances = MaximalDistance./pdist2([NucleusCentroidX(CurrentCell) NucleusCentroidY(CurrentCell)],[RandomX RandomY]);
    % Real Cell
    RealDistances = MaximalDistance./pdist2([NucleusCentroidX(CurrentCell) NucleusCentroidY(CurrentCell)],[NucleusCentroidX NucleusCentroidY]);
    RealDistances(isinf(RealDistances)) = [];
    
    SeparationRandom(CurrentCell) = sum(RandomDistances);
    SeparationTrue(CurrentCell) = sum(RealDistances);
end

AbsoluteDistances = SeparationTrue - SeparationRandom;
DistanceRandom = SeparationRandom;
DistanceReal = SeparationTrue;

LocalCCCurrentCells = AbsoluteDistances;

  
