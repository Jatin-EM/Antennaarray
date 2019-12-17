function [pts] = poissonDiscrandomarray(sizeI,spacing,nPts,showIter,lambda)


%%%%%%% Initial parameters setup
% Parsing inputs and setting default values
if nargin == 3; showIter = 0; end
if nargin == 2; showIter = 0; nPts = 0; end

% Setting properties for iterations
ndim = length(sizeI);   % Number of Dimensions
k = 5;  % Number of 'dart' tries in each grid.
dartFactor = 4; %Select number of sample data in each iterations. Change it to
% reduce run time for code. Have to play around with number. 


%%%%%%% Making Grid read for iterations
%Make grid size such that there is just one pt in each grid
dm = spacing/sqrt(ndim);    % grize cell size [Bridson 2007]
rng('shuffle')
%Make Grid
for i = 1:ndim
    sGrid{1,i} = 1:dm:sizeI(i);
end
[sGrid{:}] = ndgrid(sGrid{:});
sizeGrid = size(sGrid{1});

% Convert Grid points into a nx3 array;
for i = 1:ndim
    sGrid{i} = sGrid{i}(:);
end
sGrid = cell2mat(sGrid);

% Arrays to show eligible grids for dart throws and keeping score of darts
% thrown in a particular grid
emptyGrid = logical(ones(size(sGrid,1),1)); %Eligible Grids
nEmptyGrid = sum(emptyGrid);    %Number of eligible Grids
scoreGrid = zeros(size(emptyGrid)); %Score of darts thrown in Grid

% Darts to be thrown per iterations
% This hugely influences speed of the algorithm. Change dartFactor for it. 
if nPts == 0
    nPts = nEmptyGrid;
    ndarts = round(nEmptyGrid/dartFactor);
end
ndarts = round(nPts/dartFactor);

%%%%%%%%% Iterative process to generate points
% Initialize parameters
ptsCreated = 0;
y_vert=0.25*lambda:0.5*lambda:0.5*sizeI(2);
y_vert=[sizeI(2)/2-y_vert,sizeI(2)/2+y_vert];
x_vert=zeros(1,length(y_vert));
pts_vert=[x_vert'+lambda/4,y_vert'];
x_hor=0.75*lambda:lambda/2:sizeI(1);
% x_hor=[0,x_hor];
y_hor=zeros(1,length(x_hor));
pts_hor=[x_hor',sizeI(2)/2+y_hor'-lambda/4;x_hor',sizeI(2)/2+y_hor'+lambda/4];
pts = [pts_vert;pts_hor];
% pts=[];
iter = 0;

% Start Iterative process
tic
while ptsCreated<nPts & nEmptyGrid >0
    rng('shuffle','mlfg6331_64')
    %Thrown darts in eligible grids
    availGrid = find(emptyGrid == 1);   %Eligible grids for dart throw
    dataPts = min([nEmptyGrid,ndarts]); % Darts to be thrown

    p = datasample(availGrid,dataPts,'Replace',false); %Select grids for darts
    tempPts = sGrid(p,:) + dm*rand(length(p),ndim); %Dart throw!!!
    
    
    % Find good dart throws
    [~,D] = knnsearch([pts;tempPts],tempPts,'k',2); %Finding distance between all darts(pts)
    D = D(:,2);

    withinI = logical(prod(bsxfun(@lt,tempPts,sizeI),2)); %Eligible pts should be withing sizeI 
    eligiblePts = withinI & D>spacing & D<3.5*spacing; %elgible pts should also have minimum separation distance
    
    scorePts = tempPts(~eligiblePts,:); %Keep score from bad dart throws :(
    tempPts = tempPts(eligiblePts,:);   % Save good dart throws :)
    
    
    %Update empty Grid
    emptyPts = floor((tempPts+dm-1)/dm);
    emptyPts = num2cell(emptyPts,1);
    emptyIdx = sub2ind(sizeGrid,emptyPts{:});
    emptyGrid(emptyIdx) = 0;
    
    %Update score pts
    scorePts = floor((scorePts+dm-1)/dm);
    scorePts = num2cell(scorePts,1);
    scoreIdx = sub2ind(sizeGrid,scorePts{:});
    scoreGrid(scoreIdx) = scoreGrid(scoreIdx) + 1;
    
    %Update emptyGrid if scoreGrid has exceeded k dart throws
    emptyGrid = emptyGrid & (scoreGrid<k);
    
    %Update quantities for next iterations
    nEmptyGrid = sum(emptyGrid);
    pts = [pts;tempPts];
    ptsCreated = size(pts,1);
    iter = iter+1;
    ttoc = toc;
    
    %Display iteration details
    if showIter == 1
        disp(sprintf('Iteration: %d    Points Created: %d    EmptyGrid:%d    Total Time: %0.3f',iter,ptsCreated,nEmptyGrid,ttoc));
    end
    
end

% Cut down pts if more points are generated
if size(pts,1)>nPts
    p = 1:size(pts,1);
    p = datasample(p,nPts,'Replace',false);
    pts = pts(p,:);
end

end
