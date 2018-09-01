%% MUSIC algorithm with SpotFi method including ToF and AoA
function [estimated_aoas, estimated_tofs] = aoa_tof_music(x, ...
        antenna_distance, frequency, sub_freq_delta, plot_switch)
    if nargin == 4
        plot_switch = 1;
    end
    
    % Data covarivance matrix
    R = x * x'; 
    % Find the eigenvalues and eigenvectors of the covariance matrix
    [eigenvectors, eigenvalue_matrix] = eig(R);
%     eigenvalue_matrix对角矩阵，对角线上是特征值
    % eig������õ�����ֵ��һ���ǰ��մ�С�����˳�����еģ����԰��մ�С�����˳����������ֵ�Ͷ�Ӧ����������
    eigenvalues = diag(eigenvalue_matrix);
    [eigenvalues, eigenvectors] = sort_eigenvectors(eigenvalues, eigenvectors);
    eigenvalues_descend = sort(eigenvalues, 'descend');
    
    
    
    
    num_computed_paths = aic_new(30, 32, eigenvalues_descend);
%     num_computed_paths = 1;
%     % Find max eigenvalue for normalization
%     max_eigenvalue = max(diag(eigenvalue_matrix));
%     eigenvalue_matrix = eigenvalue_matrix / max_eigenvalue;
%     
%     % Find the largest decrease ratio that occurs between the last 10 elements (largest 10 elements)
%     % and is not the first decrease (from the largest eigenvalue to the next largest)
%     % Compute the decrease factors between each adjacent pair of elements, except the first decrease
%     start_index = size(eigenvalue_matrix, 1) - 2;
%     end_index = start_index - 10;
%     decrease_ratios = zeros(start_index - end_index + 1, 1);
%     k = 1;
%     for ii = start_index:-1:end_index
%         temp_decrease_ratio = eigenvalue_matrix(ii + 1, ii + 1) / eigenvalue_matrix(ii, ii);
%         decrease_ratios(k, 1) = temp_decrease_ratio;
%         k = k + 1;
%     end
%     [max_decrease_ratio, max_decrease_ratio_index] = max(decrease_ratios);
% 
%     index_in_eigenvalues = size(eigenvalue_matrix, 1) - max_decrease_ratio_index;
%     num_computed_paths = size(eigenvalue_matrix, 1) - index_in_eigenvalues + 1;
    %num_computed_paths=1;
    
    % Estimate noise subspace
    column_indices = 1:(size(eigenvalue_matrix, 1) - num_computed_paths);
    eigenvectors_noise = eigenvectors(:, column_indices); 
    % �ź��ӿռ�
    column_sig_indices = (size(eigenvalue_matrix, 1) - num_computed_paths+1):size(eigenvalue_matrix, 1);
    eigenvectors_sig = eigenvectors(:,column_sig_indices);
    % Peak search
    % Angle in degrees (converts to radians in phase calculations)
    %% TODO: Tuning theta too??
    theta = -90:1:90; 
    % time in milliseconds
    %% TODO: Tuning tau....
    %tau = 0:(1.0 * 10^-9):(50 * 10^-9);
    tau = 0:(2.0 * 10^-9):(400 * 10^-9);
    Pmusic = zeros(length(theta), length(tau));
    % Angle of Arrival Loop (AoA)
    for ii = 1:length(theta)
        % Time of Flight Loop (ToF)
        for jj = 1:length(tau)
            steering_vector = compute_steering_vector(theta(ii), tau(jj), ...
                    frequency, sub_freq_delta, antenna_distance);
            PP = steering_vector' * (eigenvectors_noise * eigenvectors_noise') * steering_vector;
            Pmusic(ii, jj) = abs(1 /  PP);
            Pmusic(ii, jj) = 10 * log10(Pmusic(ii, jj));% / max(Pmusic(:, jj))); 
            %Pmusic(ii, jj) = abs(Pmusic(ii, jj));
        end
    end
    %��һ��Pmusic
    if plot_switch == 1
        p_max = max(max(Pmusic));
        p_min = min(min(Pmusic));
        Pmusic = (Pmusic - p_min)./(p_max - p_min);
        % Theta (AoA) & Tau (ToF) 3D Plot
        figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off')
        mesh(tau, theta, Pmusic)
        view(2);
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Spectrum Peaks')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end

% 	% Theta (AoA)
% 	figure_name_string = sprintf('%s: Number of Paths: %d', data_name, num_computed_paths);
% 	figure('Name', figure_name_string, 'NumberTitle', 'off')
% 	plot(theta, Pmusic(:, 1), '-k')
% 	xlabel('Angle, \theta')
% 	ylabel('Spectrum function P(\theta, \tau)  / dB')
% 	title('AoA Estimation as a function of theta')
% 	grid on

%     % Find AoA peaks
%     [~, aoa_peak_indices] = findpeaks(Pmusic(:, 1));
%     estimated_aoas = theta(aoa_peak_indices);
%     
% % 	% Tau (ToF)
% %     figure_name_string = sprintf('ToF Estimation as a Function of Tau w/ AoA: %f', ...
% %             estimated_aoas(1));
% %     figure('Name', figure_name_string, 'NumberTitle', 'off')
% %     plot(tau, Pmusic(1, :), '-k')
% %     xlabel('Time of Flight \tau / degree')
% %     ylabel('Spectrum function P(\theta, \tau)  / dB')
% %     title(figure_name_string)
% %     grid on
%     
%     % Find ToF peaks
%     time_peak_indices = zeros(length(aoa_peak_indices), length(tau));
%     % AoA loop (only looping over peaks in AoA found above)
%     for ii = 1:length(aoa_peak_indices)
%         aoa_index = aoa_peak_indices(ii);
%         % For each AoA, find ToF peaks
%         [peak_values, tof_peak_indices] = findpeaks(Pmusic(aoa_index, :));
%         if isempty(tof_peak_indices)
% 
%             tof_peak_indices = 1;
%         end
%         % Pad result with -1 so we don't have a jagged matrix (and so we can do < 0 testing)
%         negative_ones_for_padding = -1 * ones(1, length(tau) - length(tof_peak_indices));
%         time_peak_indices(ii, :) = horzcat(tau(tof_peak_indices), negative_ones_for_padding);
%     end
% 
%     % Set return values
%     % AoA is now a column vector
%     estimated_aoas = transpose(estimated_aoas);
%     % ToF is now a length(estimated_aoas) x length(tau) matrix, with -1 padding for unused cells
%     estimated_tofs = time_peak_indices;
    %�Ҳ���
    [pks,locs,xymin,smin] = extrema2(Pmusic,0.5);
    estimated_aoas = theta(locs(:,1));
    estimated_tofs = tau(locs(:,2));
end

%% Computes the steering vector for SpotFi. 
function steering_vector = compute_steering_vector(theta, tau, freq, sub_freq_delta, ant_dist)
    steering_vector = zeros(30, 1);
    k = 1;
    base_element = 1;
    for ii = 1:2
        for jj = 1:15
            steering_vector(k, 1) = base_element * omega_tof_phase(tau, sub_freq_delta)^(jj - 1);
            k = k + 1;
        end
        base_element = base_element * phi_aoa_phase(theta, freq, ant_dist);
    end
end

%% Compute the phase shifts across subcarriers as a function of ToF
function time_phase = omega_tof_phase(tau, sub_freq_delta)
    time_phase = exp(-1i * 2 * pi * sub_freq_delta * tau);
end

%% Compute the phase shifts across the antennas as a function of AoA
function angle_phase = phi_aoa_phase(theta, frequency, d)
    % Speed of light (in m/s)
    c = 3.0 * 10^8;
    % Convert to radians
    theta = theta / 180 * pi;
    angle_phase = exp(-1i * 2 * pi * d * sin(theta) * (frequency / c));
end

%% ��С������������ֵ��������������Ӧ������������˳��
function [sorted_eigenvalues,sorted_eigenvectors] = sort_eigenvectors(eigenvalues, eigenvectors)
    [sorted_eigenvalues, b] = sort(eigenvalues);
    sorted_eigenvectors = eigenvectors(:,b);
end
%GET_TOTAL_RSS Calculates the Received Signal Strength (RSS) in dBm from
% a CSI struct.
%
% (c) 2011 Daniel Halperin <dhalperi@cs.washington.edu>
%
function ret = get_total_rss(csi_st)
    error(nargchk(1,1,nargin));

    % Careful here: rssis could be zero
    rssi_mag = 0;
    if csi_st.rssi_a ~= 0
        rssi_mag = rssi_mag + dbinv(csi_st.rssi_a);
    end
    if csi_st.rssi_b ~= 0
        rssi_mag = rssi_mag + dbinv(csi_st.rssi_b);
    end
    if csi_st.rssi_c ~= 0
        rssi_mag = rssi_mag + dbinv(csi_st.rssi_c);
    end
    
    ret = db(rssi_mag, 'pow') - 44 - csi_st.agc;
end