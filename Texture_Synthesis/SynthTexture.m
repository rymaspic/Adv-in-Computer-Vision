function synthIm = SynthTexture(sample, w, s)
% Texture Synthesis by Non-parameteric Sampling / Efros and Leung

% 

SEED_SIZE = 3;

if ~isa(sample,'double')
    sample = im2double(sample);
end


[sheight,swidth,nChannels] = size(sample);
theight = s(1); twidth = s(2);

synthIm = nan(theight,twidth,nChannels); 

% TODO: G = ...
sigma = 6.4;
G = fspecial('gaussian',w,w/sigma);
%G = G(:);
%G = imgaussfilt([])

% G is a 2D zero-mean Gaussian with variance w/6.4 sampled
% on a w x w grid centered about its mean


% Initialization: pick a random 3x3 patch from sample and place in the
% middle of the synthesized image.
% just for convenience, keep some space from the boundary

% i0=32; j0=4;
i0 = round(SEED_SIZE + rand() * (sheight-2*SEED_SIZE));
j0 = round(SEED_SIZE + rand() * (swidth-2*SEED_SIZE));

c = round(.5*s);
synthIm(c(1):c(1)+SEED_SIZE-1,c(2):c(2)+SEED_SIZE-1, :) = ...
                        sample(i0:i0+SEED_SIZE-1, j0:j0+SEED_SIZE-1, :);
    
% bitmap indicating filled pixels
filled = logical(zeros(s));
filled(c(1):c(1)+SEED_SIZE-1,c(2):c(2)+SEED_SIZE-1) = 1;
nFilled = sum(filled(:));
nPixels = prod(s); 

% main loop
nextP = nPixels/10;
while nFilled < nPixels
    % report progress
    if nFilled > nextP
        fprintf('%d%% complete\n', round(100*nFilled/nPixels));
        nextP = nextP+nPixels/10;
    end
        
    % TODO: find [ii, jj], locations of pixels for next round of filling-in
    % dilate current boundary
    
    % imdilate(filled)
    SE = strel('square',3);
    new_filled = imdilate(filled,SE);
    [ii,jj] = find(new_filled - filled);
    filled = new_filled;
    
    % permute (just to insert some random noise, this is not a must)
    perm = randperm(length(ii));  
    ii = ii(perm); jj = jj(perm);
    
    for i=1:length(ii)
        % TODO: set the correct template
        % place window at the center
        dd = (w-1)/2;       
        template = nan(w,w,nChannels);
        x_min = max(ii(i)-dd,1);
        x_max = min(ii(i)+dd,twidth);
        y_min = max(jj(i)-dd,1);
        y_max = min(jj(i)+dd,theight);
        x_min+1-ii(i)+dd
        x_max+twidth-ii(i)-dd
        
    
        
        %template(x_min+1-ii(i)+dd:x_max+w-ii(i)-dd,y_min+1-jj(i)+dd:y_max+w-jj(i)-dd,1:3) = synthIm(x_min:x_max,y_min:y_max,1:3);
        
        template(x_min+1-ii(i)+dd:x_max+w-ii(i)-dd,y_min+1-jj(i)+dd:y_max+w-jj(i)-dd) = synthIm(x_min:x_max,y_min:y_max);
        
        
            
        [bestMatches, errors] = FindMatches(template, sample, G);
        % sample from best matches
        pixelValue = bestMatches(:,randi(size(bestMatches,2),1)); 
        % TODO: fill in the value and update parameters= 
        
         %uncommend to run for color images:
                 for n = 1:3
                 synthIm(ii(i),jj(i),n) = pixelValue(n); 
                 end

        %uncommend to run for color images:
        %synthIm(ii(i),jj(i)) = pixelValue; 
        
        nFilled = nFilled + 1;
                    
    end
        % uncomment to show texture as it is being synthesized 
        imshow(synthIm);
        drawnow;
end
end
