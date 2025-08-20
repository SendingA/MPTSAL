function [aod_est_rad, slope_series, cleaned_diff] = extract_aod_from_pdoa(aod_pdoa_diff_wrapped, window_size, residual_thresh)
% 提取多轮 PDOA 差值中的真实 AoD 相位
% 
% 输入:
%   aod_pdoa_diff_wrapped: wrap 后的相位差序列 (rad), 长度 N
%   window_size: 滑动窗口大小 (推荐 10~20)
%   residual_thresh: 残差剔除阈值 (rad), 比如 0.3
%
% 输出:
%   aod_est_rad: 估计出的 AoD 相位差 (rad)
%   slope_series: 每个滑窗中的拟合斜率
%   cleaned_diff: 去除异常值后的 unwrap 相位差序列

    % === Step 1: 相位展开 ===
    aod_unwrapped = unwrap(aod_pdoa_diff_wrapped(:));

    N = length(aod_unwrapped);
    slope_series = nan(1, N);
    residual_series = nan(1, N);
    
    % === Step 2: 滑动线性拟合 ===
    for i = 1:N - window_size
        idx = i:(i+window_size-1);
        x = (0:window_size-1)';
        y = aod_unwrapped(idx);
        p = polyfit(x, y, 1); % 一阶拟合
        slope_series(i + floor(window_size/2)) = p(1); % 中心点保存
        y_fit = polyval(p, x);
        residual = y - y_fit;
        residual_series(i + floor(window_size/2)) = std(residual);
    end

    % === Step 3: 异常点剔除 ===
    cleaned_idx = residual_series < residual_thresh;
    cleaned_diff = aod_unwrapped(cleaned_idx & ~isnan(residual_series));

    % === Step 4: 复数平均提取最终相位 ===
    z = exp(1j * cleaned_diff);
    aod_est_rad = angle(mean(z));

    % === 可视化（可选）===
    figure;
    subplot(3,1,1);
    plot(aod_unwrapped, '-o'); hold on;
    plot(find(~cleaned_idx), aod_unwrapped(~cleaned_idx), 'rx');
    title('Unwrapped PDOA Phase Diff (剔除异常点)');
    ylabel('Phase (rad)');
    grid on;

    subplot(3,1,2);
    plot(slope_series, 'LineWidth', 1.5);
    title('Local Linear Slope (Phase Drift)');
    ylabel('dϕ/dn');
    grid on;

    subplot(3,1,3);
    plot(residual_series, 'LineWidth', 1.2);
    yline(residual_thresh, 'r--', 'Threshold');
    title('Residual (拟合残差)');
    xlabel('Packet Index');
    ylabel('Residual Std (rad)');
    grid on;
end
