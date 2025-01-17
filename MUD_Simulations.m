% code for multivariate-univariate dependence (MUD) analysis simulations
%
% looks at the relationship between univariate activity in a simulated
% brain region, pattern similarity in that region, and the relationship
% between voxels' contributions to pattern similarity and univariate
% activity
%
% variables of interest are
% 'patternActivationCorrelation','meanOverallActivation','patternSimilarity',
% which have the univariate activity, pattern similarity, and 
% MUD (patternActivationCorrelation) for each of the n iterations
%
% for more details (and results), see Aly M & Turk-Browne NB. (2016). Attention 
% stabilizes representations in the human hippocampus. Cerebral Cortex, 26, 783-796.
%
% Mariam Aly, August 2014
%
% ==========================================================================


clear all

fileName = 'Act-MUD';

for z = 1:1000 % iterations

% -a,+p,-m (simulation 6 in paper)
% parameters for negative overall activity, positive pattern similarity, and negative MUD 
% (see Aly & Turk-Browne, 2016, for all parameters)

    % voxels 1-500 are stable across the 6 patterns (and have more negative values than voxels 501-1000)
        patt1(1:500,1)=-0.25+(randn(500,1));
        patt2(1:500,1)=patt1(1:500,1);
        patt3(1:500,1)=patt1(1:500,1);
        patt4(1:500,1)=patt1(1:500,1);
        patt5(1:500,1)=patt1(1:500,1);
        patt6(1:500,1)=patt1(1:500,1);

    % the other 500 voxels are not stable across the 6 patterns: assigned random values for each pattern, and are not as negative as voxels 1-500
        patt1(501:1000,1)=(randn(500,1));
        patt2(501:1000,1)=(randn(500,1));
        patt3(501:1000,1)=(randn(500,1));
        patt4(501:1000,1)=(randn(500,1));
        patt5(501:1000,1)=(randn(500,1));
        patt6(501:1000,1)=(randn(500,1));


% put all patterns in a single variable, where rows are voxels and columns are different patterns
    patt = [patt1,patt2,patt3,patt4,patt5,patt6];
    pattsim = corrcoef(patt); % correlate each pattern with every other pattern to get pattern similarity matrix

    
% subtract the mean for each pattern (i.e., subtract the mean across rows for each column)
    patt_mean_sub = patt - mean(patt,1);
    
    
 % divide by root SS
     for i = 1:size(patt_mean_sub,2) 
         patt_norm(:,i) = patt_mean_sub(:,i)/(sqrt(sum(patt_mean_sub(:,i).^2)));
     end
 
     
 % calculate correlations to make sure no errors: correlation for each pair of patterns is the sum of the element-wise product between them
     for i = 1:size(patt_norm,2) 
         for j = 1:size(patt_norm,2)
             r(i,j) = sum(patt_norm(:,i).*patt_norm(:,j));

         end
     end
 
    check_answer = round(pattsim - r); % compare the pattern correlations above with the ones calculated earlier; should be all 0s if no mistake made
 
    
% now get voxel contributions to pattern similarity: the element-wise product of the normalized values for each pair of patterns
     for i = 1:size(patt_norm,2)
         for j = 1:size(patt_norm,2)
            patt_prod{i,j} = patt_norm(:,i).*patt_norm(:,j);
         end
     end

     
% get mean voxel contributions to each of the off-diagonal correlations
    voxelContribs = [patt_prod{1,(2:6)}];
    voxelContribs = [patt_prod{1,(2:6)},patt_prod{2,(3:6)},patt_prod{3,(4:6)},patt_prod{4,(5:6)},patt_prod{5,6}];
    meanVoxelContribs = mean(voxelContribs,2);

    
% get mean activation
    meanAct = mean(patt,2); % average across columns to get mean activity per voxel
    meanActivation = mean(meanAct); % mean activity across all voxels

    
% now calculate the correlation between activation and voxels'
% contributions to pattern similarity
    pattAct = [meanAct, meanVoxelContribs];
    pattActCorr = corrcoef(pattAct); % Fisher transform the correlations


% save values from this iteration
    patternActivationCorrelation(z) = pattActCorr(1,2); % the MUD value
    meanOverallActivation(z) = meanActivation; % mean activity across all voxels
    patternSimilarity(z) = mean(nonzeros(triu(pattsim,1))); % mean pattern similarity (upper triangle of the pattern similarity matrix)


% plot  MUD values across the z iterations
    hist(patternActivationCorrelation)
    xlabel('MUD value (r)')
    ylabel(['frequency out of ', num2str(z), ' iterations'])
    
end

save(fileName, 'patternActivationCorrelation','meanOverallActivation','patternSimilarity');

