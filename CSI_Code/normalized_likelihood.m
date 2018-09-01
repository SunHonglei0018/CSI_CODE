function [output_top_aoas] = normalized_likelihood(tof_packet_data, aoa_packet_data, num_packets,data_name)
    if nargin < 4
        data_name = ' - ';
    end
    % Find the number of elements that will be in the full_measurement_matrix
    % The value must be computed since each AoA may have a different number of ToF peaks
    full_measurement_matrix_size = 0;
    % Packet Loop
    for packet_index = 1:num_packets
        tof_matrix = tof_packet_data{packet_index};
        aoa_matrix = aoa_packet_data{packet_index};
        % AoA Loop
        for j = 1:size(aoa_matrix, 1)
            % ToF Loop
            for k = 1:size(tof_matrix(j, :), 2)
                % Break once padding is hit
                if tof_matrix(j, k) < 0
                    break
                end
                full_measurement_matrix_size = full_measurement_matrix_size + 1;
            end
        end
    end

    % Construct the full measurement matrix
    full_measurement_matrix = zeros(full_measurement_matrix_size, 2);%�������ݾ���
    full_measurement_matrix_index = 1;
    % Packet Loop
    for packet_index = 1:num_packets
        tof_matrix = tof_packet_data{packet_index};
        aoa_matrix = aoa_packet_data{packet_index};
        % AoA Loop
        for j = 1:size(aoa_matrix, 1)
            % ToF Loop
            for k = 1:size(tof_matrix(j, :), 2)
                % Break once padding is hit
                if tof_matrix(j, k) < 0
                    break
                end
                full_measurement_matrix(full_measurement_matrix_index, 1) = aoa_matrix(j, 1);
                full_measurement_matrix(full_measurement_matrix_index, 2) = tof_matrix(j, k);
                full_measurement_matrix_index = full_measurement_matrix_index + 1;
            end
        end
    end

    % Normalize AoA & ToF
%     aoa_max = max(abs(full_measurement_matrix(:, 1)));
%     tof_max = max(abs(full_measurement_matrix(:, 2)));
%     full_measurement_matrix(:, 1) = full_measurement_matrix(:, 1) / aoa_max;
%     full_measurement_matrix(:, 2) = full_measurement_matrix(:, 2) / tof_max;

    %% ����
     %���������֮��ľ���
     Y = pdist(full_measurement_matrix, 'seuclidean');
     %���������
     linkage_tree = linkage(Y, 'average');
     %�������࣬���ص���ÿ�����������������
     cluster_indices_vector = cluster(linkage_tree,'CutOff', 2.0);

    % Cluster AoA and ToF for each packet
    % Worked Pretty Well
%     linkage_tree = linkage(full_measurement_matrix, 'ward');
    % cluster_indices_vector = cluster(linkage_tree, 'CutOff', 0.45, 'criterion', 'distance');
    % cluster_indices_vector = cluster(linkage_tree, 'CutOff', 0.85, 'criterion', 'distance');
%     cluster_indices_vector = cluster(linkage_tree, 'CutOff', 1.0, 'criterion', 'distance');
    cluster_count_vector = zeros(0, 1);%��ı��
    num_clusters = 0;%һ���ж�����
    for ii = 1:size(cluster_indices_vector, 1)
        if ~ismember(cluster_indices_vector(ii), cluster_count_vector)
            cluster_count_vector(size(cluster_count_vector, 1) + 1, 1) = cluster_indices_vector(ii);
            num_clusters = num_clusters + 1;
        end
    end

    % Collect data and indices into cluster-specific cell arrays
    clusters = cell(num_clusters, 1);%ÿ���������
    cluster_indices = cell(num_clusters, 1);
    for ii = 1:size(cluster_indices_vector, 1)
        % Save off the data
        tail_index = size(clusters{cluster_indices_vector(ii, 1)}, 1) + 1;
        clusters{cluster_indices_vector(ii, 1)}(tail_index, :) = full_measurement_matrix(ii, :);
        % Save off the indexes for the data
        cluster_index_tail_index = size(cluster_indices{cluster_indices_vector(ii, 1)}, 1) + 1;
        cluster_indices{cluster_indices_vector(ii, 1)}(cluster_index_tail_index, 1) = ii;
    end

    % Delete outliers from each cluster
    for ii = 1:size(clusters, 1)
        % Delete clusters that are < 5% of the size of the number of packets
        if size(clusters{ii}, 1) < (0.05 * num_packets)
            clusters{ii} = [];
            cluster_indices{ii} = [];
            continue;
        end
        alpha = 0.05;
        [~, outlier_indices, ~] = deleteoutliers(clusters{ii}(:, 1), alpha);
        cluster_indices{ii}(outlier_indices(:), :) = [];
        clusters{ii}(outlier_indices(:), :) = [];

        alpha = 0.05;
        [~, outlier_indices, ~] = deleteoutliers(clusters{ii}(:, 2), alpha);
        cluster_indices{ii}(outlier_indices(:), :) = [];
        clusters{ii}(outlier_indices(:), :) = [];
    end

    cluster_plot_style = {'bo', 'go', 'ro', 'ko', ...
                        'bs', 'gs', 'rs', 'ks', ...
                        'b^', 'g^', 'r^', 'k^', ... 
                        'bp', 'gp', 'rp', 'kp', ... 
                        'b*', 'g*', 'r*', 'k*', ... 
                        'bh', 'gh', 'rh', 'kh', ... 
                        'bx', 'gx', 'rx', 'kx', ... 
                        'b<', 'g<', 'r<', 'k<', ... 
                        'b>', 'g>', 'r>', 'k>', ... 
                        'b+', 'g+', 'r+', 'k+', ... 
                        'bd', 'gd', 'rd', 'kd', ... 
                        'bv', 'gv', 'rv', 'kv', ... 
                        'b.', 'g.', 'r.', 'k.', ... 
                        'co', 'mo', 'yo', 'wo', ...
                        'cs', 'ms', 'ys', ...
                        'c^', 'm^', 'y^', ... 
                        'cp', 'mp', 'yp', ... 
                        'c*', 'm*', 'y*', ... 
                        'ch', 'mh', 'yh', ... 
                        'cx', 'mx', 'yx', ... 
                        'c<', 'm<', 'y<', ... 
                        'c>', 'm>', 'y>', ... 
                        'c+', 'm+', 'y+', ... 
                        'cd', 'md', 'yd', ... 
                        'cv', 'mv', 'yv', ... 
                        'c.', 'm.', 'y.', ... 
    };

    %% TODO: Tune parameters
    %% TODO: Tuning parameters using SVM results
    % Good base: 5, 10000, 75000, 0 (in order)
    % Likelihood parameters
    weight_num_cluster_points = 0.0;
    weight_aoa_variance = 0.0004;
    weight_tof_variance = -0.0016;
    weight_tof_mean = -0.0000;
    constant_offset = -1;
    %{
    weight_num_cluster_points = 0.0;
    weight_aoa_variance = -0.0010;
    weight_tof_variance = -0.0079;
    weight_tof_mean = -0.0003;
    constant_offset = -0.9997;
    %}
    %{
    % Old Likelihood parameters
    weight_num_cluster_points = 0.0001 * 10^-3;
    weight_aoa_variance = -0.7498 * 10^-3;
    weight_tof_variance = 0.0441 * 10^-3;
    weight_tof_mean = -0.0474 * 10^-3;
    constant_offset = -1;
    %}
    %constant_offset = 300;
    % Compute likelihoods
    likelihood = zeros(length(clusters), 1);
    cluster_aoa = zeros(length(clusters), 1);
    max_likelihood_index = -1;
    top_likelihood_indices = [-1; -1; -1; -1; -1;];
    for ii = 1:length(clusters)
        % Ignore clusters of size 1
        if size(clusters{ii}, 1) == 0
            continue
        end
        % Initialize variables
        num_cluster_points = size(clusters{ii}, 1);
        aoa_mean = 0;
        tof_mean = 0;
        aoa_variance = 0;
        tof_variance = 0;
        % Compute Means
        for jj = 1:num_cluster_points
            aoa_mean = aoa_mean + clusters{ii}(jj, 1);
            tof_mean = tof_mean + clusters{ii}(jj, 2);
        end
        aoa_mean = aoa_mean / num_cluster_points;
        tof_mean = tof_mean / num_cluster_points;
        % Compute Variances
        for jj = 1:num_cluster_points
            aoa_variance = aoa_variance + (clusters{ii}(jj, 1) - aoa_mean)^2;
            tof_variance = tof_variance + (clusters{ii}(jj, 2) - tof_mean)^2;
        end
        aoa_variance = aoa_variance / (num_cluster_points - 1);
        tof_variance = tof_variance / (num_cluster_points - 1);
        % Compute Likelihood
        %% TODO: Trying result from SVM
        %{
        exp_body = weight_num_cluster_points * num_cluster_points ...
                - weight_aoa_variance * aoa_variance ...
                - weight_tof_variance * tof_variance ...
                - weight_tof_mean * tof_mean ...
                - constant_offset;
        %}
        exp_body = weight_num_cluster_points * num_cluster_points ...
                + weight_aoa_variance * aoa_variance ...
                + weight_tof_variance * tof_variance ...
                + weight_tof_mean * tof_mean ...
                + constant_offset;
        likelihood(ii, 1) = exp_body;%exp(exp_body);
        % Compute Cluster Average AoA
        for jj = 1:size(clusters{ii}, 1)
            cluster_aoa(ii, 1) = cluster_aoa(ii, 1) + clusters{ii}(jj, 1);
        end
        cluster_aoa(ii, 1) = cluster_aoa(ii, 1) / size(clusters{ii}, 1);
		
		% Output
		fprintf('\nCluster Properties for cluster %d\n', ii)
		fprintf('Num Cluster Points: %d, Weighted Num Cluster Points: %d\n', ...
				num_cluster_points, (weight_num_cluster_points * num_cluster_points))
		fprintf('AoA Variance: %.9f, Weighted AoA Variance: %.9f\n', ...
				aoa_variance, (weight_aoa_variance * aoa_variance))
		fprintf('ToF Variance: %.9f, Weighted ToF Variance: %.9f\n', ...
				tof_variance, (weight_tof_variance * tof_variance))
		fprintf('AoA Mean %.9f\n', cluster_aoa(ii, 1))
		fprintf('ToF Mean: %.9f, Weighted ToF Mean: %.9f\n', ...
				tof_mean, (weight_tof_mean * tof_mean))
		fprintf('Exponential Body: %.9f\n', exp_body)
		fprintf('Cluster %d has formatting: %s, AoA: %f, and likelihood: %f\n', ...
				ii, cluster_plot_style{ii}, cluster_aoa(ii, 1), likelihood(ii, 1))
				
        % Check for maximum likelihood
        if max_likelihood_index == -1 ...
                || likelihood(ii, 1) > likelihood(max_likelihood_index, 1)
            max_likelihood_index = ii;
        end
        % Record the top maximum likelihoods
        for jj = 1:size(top_likelihood_indices, 1)
            % Replace empty slot
            if top_likelihood_indices(jj, 1) == -1
                top_likelihood_indices(jj, 1) = ii;
                break;
            % Add somewhere in the list
            elseif likelihood(ii, 1) > likelihood(top_likelihood_indices(jj, 1), 1)
                % Shift indices down
                for kk = size(top_likelihood_indices, 1):-1:(jj + 1)
                    top_likelihood_indices(kk, 1) = top_likelihood_indices(kk - 1, 1);
                end
                top_likelihood_indices(jj, 1) = ii;
                break;
            % Add an extra item to the list because the likelihoods are all equal...
            elseif likelihood(ii, 1) == likelihood(top_likelihood_indices(jj, 1), 1) ...
                    && jj == size(top_likelihood_indices, 1)
                top_likelihood_indices(jj + 1, 1) = ii;
                break;
            % TODO: Make sure I want to keep this
            elseif jj == size(top_likelihood_indices, 1) 
                top_likelihood_indices(jj + 1, 1) = ii;
                break;
            end
        end
    end
	
	fprintf('\nThe cluster with the maximum likelihood is cluster %d\n', max_likelihood_index)
	fprintf('The top clusters are: ')
	for jj = 1:size(top_likelihood_indices, 1)
		if top_likelihood_indices(jj, 1) ~= -1
			fprintf('%d, ', top_likelihood_indices(jj, 1));
		end
	end
	fprintf('\n')
	fprintf('The clusters have AoAs: ')
	for jj = 1:size(top_likelihood_indices, 1)
		if top_likelihood_indices(jj, 1) ~= -1
			fprintf('%g, ', cluster_aoa(top_likelihood_indices(jj, 1), 1));
		end
	end
	fprintf('\n')
	fprintf('The clusters have plot point designations: ')
	for jj = 1:size(top_likelihood_indices, 1)
		if top_likelihood_indices(jj, 1) ~= -1
			fprintf('%s, ', ...
					cluster_plot_style{top_likelihood_indices(jj, 1)})
		end
	end
	fprintf('\n')
	fprintf('The clusters have likelihoods: ')
	for jj = 1:size(top_likelihood_indices, 1)
		if top_likelihood_indices(jj, 1) ~= -1
			fprintf('%g, ', ...
					likelihood(top_likelihood_indices(jj, 1), 1))
		end
	end
	fprintf('\n')

    %Plot AoA & ToF
	figure_name_string = sprintf('%s: AoA vs ToF Plot', data_name);
	figure('Name', figure_name_string, 'NumberTitle', 'off')
	hold on
	% Plot the data from each cluster and draw a circle around each cluster 
	for ii = 1:size(cluster_indices, 1)
		% Some clusters may all be considered outliers...we'll see how this affects the 
		% results, but it shouldn't result in a crash!!
		if isempty(cluster_indices{ii})
			continue
		end
		% If there's an index out of bound error, 'cluster_plot_style' is to blame
		plot(full_measurement_matrix(cluster_indices{ii, 1}, 1), ...
				full_measurement_matrix(cluster_indices{ii}, 2), ...
				cluster_plot_style{ii}, ...
				'MarkerSize', 8, ...
				'LineWidth', 2.5)
		% Compute Means
		num_cluster_points = size(cluster_indices{ii}, 1);
		aoa_mean = 0;
		tof_mean = 0;
		for jj = 1:num_cluster_points
			% Denormalize AoA for presentation
			aoa_mean = aoa_mean + ( clusters{ii}(jj, 1));
			tof_mean = tof_mean + clusters{ii}(jj, 2);
		end
		aoa_mean = aoa_mean / num_cluster_points;
		tof_mean = tof_mean / num_cluster_points;
		cluster_text = sprintf('Size: %d', num_cluster_points);
		text(aoa_mean + 5, tof_mean, cluster_text);
		drawnow
    end
    
    % plot cycle for top 5 cluster
    rank_index=0;
    for jj = 1:size(top_likelihood_indices, 1)
		if top_likelihood_indices(jj, 1) ~= -1
            top_index=top_likelihood_indices(jj, 1);
			[ellipse_x, ellipse_y] = compute_ellipse(...
			clusters{top_index}(:, 1) , ...
			clusters{top_index}(:, 2));
            rank_index=rank_index+1;
            rank_text=sprintf('Rank: %d, likelihood: %g',rank_index,likelihood(top_index, 1));
            plot(ellipse_x, ellipse_y, 'c-', 'LineWidth', 3);
            text(ellipse_x(1,size(ellipse_x,2)*3/4), ellipse_y(1,size(ellipse_y,2)*3/4), rank_text);
            drawnow;
		end
    end
    
	xlabel('Angle of Arrival (AoA)')
	ylabel('(Normalized) Time of Flight (ToF)')
	title('Angle of Arrival vs. (Normalized) Time of Flight')
	grid on
	hold off

    % Select AoA
    max_likelihood_average_aoa = cluster_aoa(max_likelihood_index, 1);
	fprintf('The Estimated Angle of Arrival for data set %s is %f\n', ...
			data_name, max_likelihood_average_aoa)
    % Profit
    
    % Trim remaining -1 indices from the end
    ii = size(top_likelihood_indices, 1);
    while ii > 0
        if top_likelihood_indices(ii, 1) == -1
            top_likelihood_indices(ii, :) = [];
            ii = ii - 1;
        else
            break;
        end
    end
%     top_likelihood_indices
    output_top_aoas = cluster_aoa(top_likelihood_indices);
end

%% Computes ellipse totally enclosing the points defined by x and y, includes room for markers, etc.
% x         -- the x coordinates of the points that the ellipse should enclose
% y         -- the y coordinates of the points that the ellipse should enclose
% Return:
% ellipse_x -- the x coordinates of the enclosing ellipse
% ellipse_y -- the y coordinates of the enclosing ellipse
function [ellipse_x, ellipse_y] = compute_ellipse(x, y)
    % Buffer room for each dimension
    marker_adjustment_quantity_x = 4;
    marker_adjustment_quantity_y = 0.05;
    % Find centroid
    centroid_x = sum(x) / length(x);
    centroid_y = sum(y) / length(y);
    % Find max difference between points in each dimension (diameter)
    % x
    diameter_x = 0;
    for ii = 1:length(x)
        for jj = 1:length(x)
            if abs(x(ii) - x(jj)) > diameter_x
                diameter_x = abs(x(ii) - x(jj));
            end
        end
    end
    radius_x = diameter_x / 2;
    radius_x = radius_x + marker_adjustment_quantity_x; 
    % y
    diameter_y = 0;
    for ii = 1:length(y)
        for jj = 1:length(y)
            if abs(y(ii) - y(jj)) > diameter_y
                diameter_y = abs(y(ii) - y(jj));
            end
        end
    end
    radius_y = diameter_y / 2;
    radius_y = radius_y + marker_adjustment_quantity_y;
    % Generate points of ellipse
    t = 0:0.001:(2 * pi);
    ellipse_x = radius_x * cos(t) + centroid_x;
    ellipse_y = radius_y * sin(t) + centroid_y;
end
