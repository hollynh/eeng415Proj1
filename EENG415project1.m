% Project 1: Equity in Healthcare
% Holly Hammons
% EENG 415
% Dr. Salman Mohagheghi
% Spring 2022
%% data pre-processing
% read in the raw data
dataIn = readtable('eeng415Proj1Data.xlsx');
% change table to array
dataIn = table2array(dataIn);

% document which rows contain missing values
[rowsSize, colSize] = size(dataIn);
nanMat = isnan(dataIn);
dataRemoved = zeros(rowsSize, 1);
count = 0;
for i = 1:rowsSize
    for j = 1: colSize
        if nanMat(i, j) == 1
            count = count + 1;
            dataRemoved(count) = dataIn(i, 1);
            break;
        end
    end
end

% this vector contains the tracts that were removed from data set
dataRemoved = dataRemoved(1:count);

% convert table to an array, only keep numerical data
rawData = dataIn(:, 3:end);

% remove rows that contain missing data
rawData = rmmissing(rawData);
%% Observe Variances/Remove Features with Low Variability Relative to Others
varMat = var(rawData);
coeffOfVar = zeros((colSize - 2), 1);

% calculate coefficient of variation
for i = 1:(colSize - 2)
    coeffOfVar(i) = (std(rawData(:, i)) / mean(rawData(:, i))) * 100;
end

maxVariation = max(coeffOfVar);
featToRem = zeros((colSize-2), 1);

% compare variations to max variation
for i = 1:size(coeffOfVar)
    if coeffOfVar(i) <= (0.02 * maxVariation)
        featToRem(i) = i;
    end
end
featToRem = sort(featToRem, 'descend');

% remove features that have coefficients of variance that are 2% or less
% of the maximum coefficient of variance
for i = 1:size(featToRem)
    if featToRem(i) ~= 0
        if featToRem(i) == 1
            rawData = rawData(:, 2:end);
        else 
            rawData = rawData(:, [1:(featToRem(i)-1), (featToRem(i)+1):end]);
        end
    else
        break;
    end
end
%% PCA
% normalize the data
rawData = normalize(rawData);

% find the covariance matrix
cVar = cov(rawData);
% find the eigenvectors and eigenvalues of the covariance matrix
[eigVec, temp] = eig(cVar);
eigVals = eig(cVar);

% match the eigenvector to the order of the largest-smallest 
% eigenvalues.
[rows,cols] = size(eigVals);
sortingTemp = zeros(rows, cols + 1);
% create a matrix to keep track of order of eigenvalues
for i = 1:rows
    sortingTemp(i, 1) = eigVals(i);
    sortingTemp(i, 2) = i;
end
% sort eigenvalues based on previously defined order
eigVals = sort(eigVals, 'descend');
sortingTemp = sortrows(sortingTemp, 1, 'descend');
[rowsV,colsV] = size(eigVec);
newEigVec = zeros(rowsV, colsV);
for i = 1:rows
    newEigVec(:, i) = eigVec(:, sortingTemp(i, 2));
end

% find the variances by normalizing eigenvalues
variances = eigVals / sum(eigVals);
variances = sort(variances, 'descend');

figure(1)
bar(variances)
title('Variability in Data Explained vs Principle Components'); 
ylabel('Percentage of Variability in Data Explained');
xlabel('Principle Components');

% sum up variances until at least 0.9
sumVar = 0;
i = 1;
numVar = 0;
while (sumVar <= 0.9000)
    sumVar = sumVar + variances(i);
    i = i + 1;
    numVar = numVar + 1;
end
newEigVec = newEigVec(:, [1:numVar]);

% transform data points into an N x numVar matrix
pcaData = newEigVec' * rawData';
pcaData = pcaData';

fprintf('The first 10 datapoints of the transformed data after PCA are:\n');
for i = 1:10
    fprintf('%f %f \n', pcaData(i,:));
    fprintf('\n');
end
%% Part E: Computing Vulnerability Level
[rowNew, colNew] = size(pcaData);
sumMat = zeros(rowNew, 1);
for i = 1:rowNew
    for j = 1:colNew
        sumMat(i, 1) = sumMat(i, 1) + pcaData(i, j);
    end
end

% plot vulnerability sums to view approximate distribution
figure(2)
histogram(sumMat);
title('Vulnerabilities by Tract');
xlabel('Vulnerability Value');
ylabel('Percentage of Tracts Defined by Vulnerability');

% q q plot of data to check for normality
figure(3)
qqplot(sumMat);
title('QQ Plot of Vulnerability Levels');

% Use Anderson-Darling test for Normality
[adt,p,adstat,cv] = adtest(sumMat);
test = 0;
if adt == 0
    test = "normal";
elseif adt == 1
    test = "not normal";    
end
fprintf('\n');
fprintf('The test statistic for the Anderson-Darling test is: %f. \n', adstat);
fprintf('The critical value for the Anderson-Darling test is: %f. \n', cv);
fprintf('The Anderson-Darling Test result is: %s. \n', test);

% calculate mean and standard deviation to use in vulnerability analysis
meanVul = mean(sumMat);
stdVul = std(sumMat);

tractVul = zeros(rowNew, 2);
for i = 1:rowNew
    tractVul(i, 1) = i;
end

% assign values 1-5 for least to most vulnerable
for i = 1:rowNew
    if sumMat(i) < (meanVul - (2*stdVul))
        tractVul(i, 2) = 1;
    elseif (sumMat(i) > (meanVul - (2*stdVul))) && (sumMat(i) < (meanVul - stdVul))
        tractVul(i, 2) = 2;
    elseif (sumMat(i) > (meanVul - stdVul)) && (sumMat(i) < (meanVul + stdVul))
        tractVul(i, 2) = 3;
    elseif (sumMat(i) > (meanVul + stdVul)) && (sumMat(i) < (meanVul + (2*stdVul)))
        tractVul(i, 2) = 4;
    elseif sumMat(i) > (meanVul + (2*stdVul)) % is 3*std a typo?
        tractVul(i, 2) = 5;
    end
end

sTemp = sortrows(tractVul, 2, 'descend');

[B,I] = maxk(sumMat, 5);

fprintf('\nThe 5 most vulnerable tracts are: \n');
for i = 1:5
    fprintf('Tract number: %d \n', dataIn(tractVul(I(i), 1), 2));
end

% create final table showing vulnerability level for each tract
finalTractVul = zeros(rowNew, 2);
for i = 1:rowNew
    finalTractVul(i, 1) = dataIn(i, 2);
end

% matrix containing vulnerability ranking for each tract
for i = 1:rowNew
    finalTractVul(i, 2) = tractVul(i, 2);
end

writematrix(finalTractVul,'tractsVuln23.csv') 
