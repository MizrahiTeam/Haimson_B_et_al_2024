% SVM_CLASSIFIER - A function to perform Support Vector Machine (SVM) classification.
% Inputs:
%   D: Cell array containing two matrices, where each matrix represents a class
%      with neurons as rows and trials as columns.
%   iter: Number of iterations for the outer loop (subsampling in case of different trial number for the two classes).
%   iter2: Number of iterations for the inner loop.
%   frac: Fraction of data to be held out for testing.
% Outputs:
%   Ac: Accuracy of the SVM classifier on original data for each iteration.
%   Ac_shuffled: Accuracy of the SVM classifier on shuffled data for each iteration.

function [Ac, Ac_shuffled] = SVM_classifer(D, iter, iter2, frac)

% Extract data for each class
x1 = D{1};
x2 = D{2};

% Check if the dimensions of the two classes are not equal
if size(x1, 2) ~= size(x2, 2)
    % Initialize arrays to store accuracies
    Ac = zeros(iter, iter2);
    Ac_shuffled = zeros(iter, iter2);
    
    % Determine the minimum number of samples between the two classes
    MIN = min([size(x1, 2), size(x2, 2)]);
    
    % Perform iterations
    for i = 1:iter
        % Randomly select MIN samples from each class
        if size(x1, 2) > size(x2, 2)
            v = randperm(size(x1, 2), MIN);
            x1 = x1(:, v);
        elseif size(x2, 2) > size(x1, 2)
            v = randperm(size(x2, 2), MIN);
            x2 = x2(:, v);
        end
        
        % Concatenate the samples from both classes
        X = [x1, x2];
        % Assign labels (1 for class 1, 0 for class 2)
        Y = [ones(size(x1)), zeros(size(x2))];
        % Convert labels to nominal array
        Y_nom = nominal(Y(1, :));
        
        % Perform inner iterations
        for j = 1:iter2
            % Shuffle the labels
            Y_nom_shuff = nominal(Y(1, randperm(numel(Y(1, :)), numel(Y(1, :))), :));
            
            % Create cross-validation partitions for original and shuffled data
            P = cvpartition(Y_nom, 'Holdout', frac, 'Stratify', false);
            P_sh = cvpartition(Y_nom_shuff, 'Holdout', frac, 'Stratify', false);
            
            % Train SVM model on original data
            SVMModel = compact(fitcsvm(X(:, P.training)', Y_nom(P.training), 'KernelFunction', 'linear', 'KernelScale', 'auto'));
            % Predict labels for test data
            [label, ~] = predict(SVMModel, X(:, P.test)');
            % Compute accuracy
            Ac(i, j) = mean(Y_nom(P.test)' == label);
            
            % Train SVM model on shuffled data
            [label, ~] = predict(SVMModel, X(:, P_sh.test)');
            % Compute accuracy
            Ac_shuffled(i, j) = mean(Y_nom_shuff(P_sh.test)' == label);
        end
    end
    
else
    % Concatenate the samples from both classes
    X = [x1, x2];
    % Assign labels (1 for class 1, 0 for class 2)
    Y = [ones(size(x1)), zeros(size(x2))];
    % Convert labels to nominal array
    Y_nom = nominal(Y(1, :));
    
    % Initialize arrays to store accuracies
    Ac_shuffled = zeros(iter2, 1);
    Ac = zeros(iter2, 1);
    
    % Perform inner iterations
    for i = 1:iter2
        % Shuffle the labels
        Y_nom_shuff = nominal(Y(1, randperm(numel(Y(1, :)), numel(Y(1, :))), :));
        
        % Create cross-validation partitions for original and shuffled data
        P = cvpartition(Y_nom, 'Holdout', frac, 'Stratify', false);
        P_sh = cvpartition(Y_nom_shuff, 'Holdout', frac, 'Stratify', false);
        
        % Train SVM model on original data
        SVMModel = compact(fitcsvm(X(:, P.training)', Y_nom(P.training), 'KernelFunction', 'linear', 'KernelScale', 'auto'));
        % Predict labels for test data
        [label, ~] = predict(SVMModel, X(:, P.test)');
        % Compute accuracy
        Ac(i) = mean(Y_nom(P.test)' == label);
        
        % Train SVM model on shuffled data
        [label, ~] = predict(SVMModel, X(:, P_sh.test)');
        % Compute accuracy
        Ac_shuffled(i) = mean(Y_nom_shuff(P_sh.test)' == label);
    end
end

end
