function [bestMatches, errors] = FindMatches(template, sample, G)

% parameters, as used by Efros and Leung
epsilon = 0.1; 
delta = 0.3;

% set the mask variable, need normalizaiton
% validMask is a square mask of width w that is 1 where template is
% filled
% partition sample to blocks (represented by column vectors).
% We can actually do this only once, and pass this representation to this 
% function, but I leave it as is in order not to change function signature 
% that was instructed.

%uncommend to run for color images
for i=1:3
    validMask(:,:,i) = ~isnan(template(:,:,i));
    size(G)
    size(validMask(:,:,i))
    mask(:,:,i) = G.*validMask(:,:,i);
    sum_mask = sum(mask(:,:,i));
    mask(:,:,i) = mask(:,:,i)./sum_mask;    
end
mask = mask(:);
size(mask);

% %uncommend to run for bw images
% validMask = ~isnan(template);
% validMask = validMask(:);
% G = G(:);
% mask = G.*validMask;
% sum_mask = sum(mask);
% mask = mask/sum_mask;


[theight,twidth,nChannels] = size(template);
templateSize = [theight,twidth];
s = [];

for i=1:nChannels
    s = [s ; im2col(sample(:,:,i), templateSize, 'sliding')];
end

nBlocks = size(s,2);
template = template(:);

% some vectorized code that calcualtes SSD(template,sample)*mask for all
% patches
% TODO: find best blocks

errors = nansum(repmat(mask,[1,nBlocks]).*(s - repmat(template,[1,nBlocks])).^2);
errors;
bestBlocks = find(errors<=min(errors)*(1+epsilon) & errors<=delta)
find(errors);
if isempty(bestBlocks)
     bestBlocks = randperm(length(find(errors)),1);
     %bestBlocks = find(randperm(errors),1)    
end


centerPixel = sub2ind(templateSize,round(templateSize(1)/2),round(templateSize(2)/2));
centerPixel(2:nChannels) = centerPixel + prod(templateSize)*[1:nChannels-1];
bestMatches = s(centerPixel,bestBlocks)


