
Shuffle_num = 1000;

% Loop through each mouse
for m = 1:n_mice
    % Loop through each field for the current mouse
    for f = 1:n_field(m)
        % Extract data for current mouse and field
        X = CR_e_tr_all{m,f};
        Y = Hit_e_tr_all{m,f};
        % Iterate over blocks
        for b = 1:numel(X)
            % Iterate over neurons
            for n = 1:size(Hit_e_tr{m,f},1)
                % Extract data for current trial and neuron
                D1 = Y{b}(n,:);
                D2 = X{b}(n,:);
                % Check if data is not empty
                if ~isempty(D1) && ~isempty(D2)
                    % Check if number of data points in two conditions is not equal
                    if numel(D1) ~= numel(D2)
                        % Determine minimum number of data points
                        MIN = min([numel(D1),numel(D2)]);
                        % Initialize array to store AUC values
                        x = zeros(1,1000);
                        % Perform shuffling and subsampling
                        for i = 1:1000
                            if numel(D1) > numel(D2)
                                v = randperm(numel(D1),MIN);
                                D1 = D1(:,v);
                            elseif numel(D2) > numel(D1)
                                v = randperm(numel(D2),MIN);
                                D2 = D2(:,v);
                            end
                            % Generate binary index for AUC calculation
                            index = [ones(1,numel(D1)),zeros(1,numel(D2))];
                            r = [D1,D2];
                            [~,~,~,x(i)] = perfcurve(logical(index),r,'true');
                        end
                        % Compute mean AUC
                        AUC_e_b{m,f,b}(n) = mean(0.5 + abs(x - 0.5));
                        % Compute mean AUC for shuffled data
                        AUC_shuffled_trail = zeros(1,Shuffle_num);
                        for k = 1:Shuffle_num
                            shuffled_index = index(randperm(length(index)));
                            [~,~,~,AUC_shuffled_trail(k)] = perfcurve(logical(shuffled_index),r,'true');
                        end
                        % Compute mean shuffled AUC
                        AUC_e_shuff_b{m,f,b}(n) = mean(AUC_shuffled_trail);
                    else
                        % Generate binary index for AUC calculation
                        index = [ones(1,numel(D1)),zeros(1,numel(D2))];
                        r = [D1,D2];
                        % Compute AUC
                        [~,~,~,AUC_e_b{m,f,b}(n)] = perfcurve(logical(index),r,'true');
                        % Compute mean AUC for shuffled data
                        AUC_shuffled_trail = zeros(1,Shuffle_num);
                        for k = 1:Shuffle_num
                            shuffled_index = index(randperm(length(index)));
                            [~,~,~,AUC_shuffled_trail(k)] = perfcurve(logical(shuffled_index),r,'true');
                        end
                        % Compute mean shuffled AUC
                        AUC_e_shuff_b{m,f,b}(n) = mean(AUC_shuffled_trail);                        
                    end
                end
            end
        end
        % Repeat the process for 'CR_h_tr_all' and 'Hit_h_tr_all'
        X = CR_h_tr_all{m,f};
        Y = Hit_h_tr_all{m,f};
        for b = 1:numel(X)
            for n = 1:size(Hit_e_tr{m,f},1)
                D1 = Y{b}(n,:);
                D2 = X{b}(n,:);
                if ~isempty(D1) && ~isempty(D2)
                    if numel(D1) ~= numel(D2)
                        MIN = min([numel(D1),numel(D2)]);
                        x = zeros(1,1000);
                        for i = 1:1000
                            if numel(D1) > numel(D2)
                                v = randperm(numel(D1),MIN);
                                D1 = D1(:,v);
                            elseif numel(D2) > numel(D1)
                                v = randperm(numel(D2),MIN);
                                D2 = D2(:,v);
                            end
                            index = [ones(1,numel(D1)),zeros(1,numel(D2))];
                            r = [D1,D2];
                            [~,~,~,x(i)] = perfcurve(logical(index),r,'true');
                        end
                        AUC_h_b{m,f,b}(n) = mean(0.5 + abs(x - 0.5));
                        AUC_shuffled_trail = zeros(1,Shuffle_num);
                        for k = 1:Shuffle_num
                            shuffled_index = index(randperm(length(index)));
                            [~,~,~,AUC_shuffled_trail(k)] = perfcurve(logical(shuffled_index),r,'true');
                        end
                        AUC_h_shuff_b{m,f,b}(n) = mean(AUC_shuffled_trail);
                    else
                        index = [ones(1,numel(D1)),zeros(1,numel(D2))];
                        r = [D1,D2];
                        [~,~,~,AUC_h_b{m,f,b}(n)] = perfcurve(logical(index),r,'true');
                        AUC_shuffled_trail = zeros(1,Shuffle_num);
                        for k = 1:Shuffle_num
                            shuffled_index = index(randperm(length(index)));
                            [~,~,~,AUC_shuffled_trail(k)] = perfcurve(logical(shuffled_index),r,'true');
                        end
                        AUC_h_shuff_b{m,f,b}(n) = mean(AUC_shuffled_trail);                       
                    end                   
                end
            end
        end
    end 
end

% Adjust AUC values
AUC_e_b = cellfun(@(x) 0.5 + abs(x - 0.5), AUC_e_b, 'UniformOutput', false);
AUC_e_shuff_b = cellfun(@(x) 0.5 + abs(x - 0.5), AUC_e_shuff_b, 'UniformOutput', false);
AUC_h_b = cellfun(@(x) 0.5 + abs(x - 0.5), AUC_h_b, 'UniformOutput', false);
AUC_h_shuff_b = cellfun(@(x) 0.5 + abs(x - 0.5), AUC_h_shuff_b, 'UniformOutput', false);

% Compute average AUC values
for m = 1:n_mice
    for f = 1:n_field(m)
        AUC_e_avg{m,f} = nanmean(cat(1,AUC_e_b{m,f,:}));
        AUC_h_avg{m,f} = nanmean(cat(1,AUC_h_b{m,f,:}));
        AUC_e_shuff_avg{m,f} = nanmean(cat(1,AUC_e_shuff_b{m,f,:}));
        AUC_h_shuff_avg{m,f} = nanmean(cat(1,AUC_h_shuff_b{m,f,:}));
    end
end
