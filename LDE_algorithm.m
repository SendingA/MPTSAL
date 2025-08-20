% function [fp_index, fp_index_fractional] = LDE_algorithm(cir, window_size, lde_thresh, required_consecutive)
% 
% % 输入：
% % cir: 复数CIR数据（长度N）
% % window_size: 滑动窗口大小 (建议4~8)
% % lde_thresh: 阈值系数 (建议3~6)
% % required_consecutive: 连续超阈tap数 (建议4~8)
% 
% % 输出：
% % fp_index：整数tap（LDE检测结果）
% % fp_index_fractional：亚tap精度结果（通过抛物线拟合）
% 
%     % 1. 计算CIR能量（功率）
%     cir_power = abs(cir).^2;
% 
%     % 2. 滑动均值滤波 (平滑噪声)
%     cir_power_smooth = movmean(cir_power, window_size);
% 
%     % 3. 估计噪声阈值
%     noise_floor = mean(cir_power_smooth(1:100)); % 取前100个tap为噪声
%     threshold = noise_floor * lde_thresh;
% 
%     % 4. First Path检测
%     consecutive_count = 0;
%     fp_index = -1;
% 
%     for i = 1:length(cir_power_smooth)
%         if cir_power_smooth(i) > threshold
%             consecutive_count = consecutive_count + 1;
%             if consecutive_count >= required_consecutive
%                 fp_index = i - required_consecutive + 1; % 取第一次连续超阈的起点
%                 break;
%             end
%         else
%             consecutive_count = 0;
%         end
%     end
% 
%     if fp_index == -1
%         warning('No First Path detected!');
%         fp_index_fractional = NaN;
%         return;
%     end
% 
%     % 5. 抛物线拟合（亚tap位置计算）
%     if fp_index > 1 && fp_index < length(cir_power)-1
%         y1 = cir_power(fp_index - 1);
%         y2 = cir_power(fp_index);
%         y3 = cir_power(fp_index + 1);
% 
%         delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
%         fp_index_fractional = fp_index + delta;
%     else
%         fp_index_fractional = fp_index; % 边界点，无插值
%     end
% 
%     % 6. 画图显示（可选）
%     figure;
%     plot(10*log10(cir_power),'b'); hold on;
%     plot(fp_index,10*log10(cir_power(fp_index)),'ro','MarkerSize',8,'LineWidth',2);
%     plot(fp_index_fractional,10*log10(cir_power(round(fp_index_fractional))),'gx','MarkerSize',8,'LineWidth',2);
%     title('LDE Algorithm First Path Detection');
%     xlabel('CIR Tap Index'); ylabel('Power (dB)');
%     legend('CIR Power','Integer FP','Fractional FP');
%     grid on;
% 
% end
function [path_indices, path_indices_fractional] = LDE_algorithm(cir, window_size, lde_thresh, required_consecutive, N_paths, min_separation)
    cir_power = abs(cir).^2;
    cir_power_smooth = movmean(cir_power, window_size);
    
    noise_floor = mean(cir_power_smooth(1:100)); % 可调区域
    threshold = noise_floor * lde_thresh;

    % --- 找出所有可能路径位置 ---
    candidates = find(cir_power_smooth > threshold);

    % --- 提取局部最大值点 ---
    peaks = [];
    for i = 2:length(cir_power_smooth)-1
        if cir_power_smooth(i) > threshold && ...
           cir_power_smooth(i) > cir_power_smooth(i-1) && ...
           cir_power_smooth(i) > cir_power_smooth(i+1)
            peaks(end+1) = i; %#ok<AGROW>
        end
    end

    % --- 按功率排序（强->弱） ---
    peak_powers = cir_power(peaks);
    [~, sort_idx] = sort(peak_powers, 'descend');
    peaks_sorted = peaks(sort_idx);

    % --- 选择前N个且间隔不小于min_separation ---
    path_indices = [];
    for i = 1:length(peaks_sorted)
        idx = peaks_sorted(i);
        if isempty(path_indices) || all(abs(path_indices - idx) >= min_separation)
            path_indices(end+1) = idx; %#ok<AGROW>
        end
        if length(path_indices) >= N_paths
            break;
        end
    end

    % --- 抛物线拟合提升精度 ---
    path_indices_fractional = zeros(size(path_indices));
    for k = 1:length(path_indices)
        idx = path_indices(k);
        if idx > 1 && idx < length(cir_power)-1
            y1 = cir_power(idx - 1);
            y2 = cir_power(idx);
            y3 = cir_power(idx + 1);
            delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
            path_indices_fractional(k) = idx + delta;
        else
            path_indices_fractional(k) = idx; % 无法拟合
        end
    end

    % --- 可选画图 ---
    figure;
    plot(10*log10(cir_power), 'b'); hold on;
    plot(path_indices, 10*log10(cir_power(path_indices)), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    plot(path_indices_fractional, 10*log10(cir_power(round(path_indices_fractional))), 'gx', 'MarkerSize', 8, 'LineWidth', 2);
    title('LDE Multi-path Detection');
    xlabel('CIR Tap Index'); ylabel('Power (dB)');
    legend('CIR Power','Integer Delay','Fractional Delay');
    grid on;
end
