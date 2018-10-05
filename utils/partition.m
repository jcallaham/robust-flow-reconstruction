function [Train, Test, mean_flow, mask] = partition(X, frac)
%[Train, Test, mean_flow, mask] = PARTITION(X, frac)
%  Split data X into training and test set (by time)
%  and return both, subtracted by mean of the training set
%  mask will be true for any NaN entries in X (these are set to zero)
%  If frac is an integer greater than 1, will be interpreted as mTrain

% Partition into training and test

if frac < 1
    mTrain = round(frac*size(X, 2));
else
    mTrain = frac;
end

Train = X(:, 1:mTrain);
mean_flow = mean(Train, 2);

Train = bsxfun(@minus, Train, mean_flow);
Test = bsxfun(@minus, X(:, mTrain+1:end), mean_flow);

mask = any(isnan(Train), 2) & any(isnan(Test), 2);
Train(mask, :) = 0;
Test(mask, :) = 0;
mean_flow(mask) = 0;

end

