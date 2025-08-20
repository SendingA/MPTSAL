% function [ToA_est, powers, p_vec, scan_axis] = estimate_ToA_SAMV3(CIR, N_path, step)
% % 使用 SAMV-3 算法估计 CIR 中多个路径的 ToA（支持亚采样精度）
% %
% % 输入：
% %   CIR     ：N×T 矩阵（N 为 CIR 长度，T 为快照数，可为1）
% %   N_path  ：希望估计的路径数（例如直射+反射 = 2）
% %   step    ：亚采样精度步长（如0.1表示十分之一采样点）
% %
% % 输出：
% %   ToA_est ：1×N_path，亚采样级别的 ToA 索引（单位：采样点）
% %   powers  ：每个路径对应的功率
% %   p_vec   ：估计出的功率谱向量（length(scan_axis)×1）
% %   scan_axis ：扫描的延迟坐标（单位：采样点）
% 
% % 默认值
% if nargin < 3
%     step = 0.1;
% end
% 
% % === 1. 预处理 ===
% CIR = double(CIR);
% if isvector(CIR)
%     CIR = CIR(:);  % 单条 CIR，转列向量
% end
% 
% [N, T] = size(CIR);
% 
% % === 2. 构建亚采样延迟模板 ===
% scan_axis = 0:step:(N-1);  % 扫描位置（单位：采样点）
% K = length(scan_axis);
% A = zeros(N, K);
% 
% for k = 1:K
%     delay = scan_axis(k);
%     A(:,k) = sinc((0:N-1)' - delay);  % Sinc 插值模板
% end
% 
% % === 3. 初始估计（DAS） ===
% DAS_init = mean(abs(CIR).^2, 2);  % 平均能量
% DAS_init_interp = A' * DAS_init; % 投影到扫描位置
% 
% % === 4. 构建自相关矩阵 ===
% R_N = (CIR * CIR') / T;
% sigma = mean(abs(CIR(:)).^2);   % 初始噪声功率
% p_vec_Old = abs(DAS_init_interp).^2;
% 
% % === 5. SAMV-3 迭代 ===
% threshold = 1e-6;
% maxIter = 30;
% 
% for iter = 1:maxIter
%     R = A * diag(p_vec_Old) * A' + sigma * eye(N);
%     Rinv = inv(R);
%     RinvA = Rinv * A;
% 
%     num = sum(conj(A) .* (Rinv * R_N * RinvA), 1).';
%     den = sum(conj(A) .* RinvA, 1).';
% 
%     p_vec = p_vec_Old .* (num ./ (den + 1e-12));
%     sigma = real(trace(Rinv * Rinv * R_N)) / real(trace(Rinv * Rinv));
% 
%     if norm(p_vec - p_vec_Old) / norm(p_vec_Old) < threshold
%         break;
%     end
%     p_vec_Old = p_vec;
% end
% 
% p_vec = real(p_vec);
% 
% % === 6. 选择最强 N_path 个峰 ===
% [~, idxs] = findpeaks(p_vec, 'SortStr', 'descend');
% 
% if length(idxs) < N_path
%     warning('未检测到足够路径，结果可能不可靠');
%     idxs = sort(idxs, 'ascend');
% else
%     idxs = sort(idxs(1:N_path));  % 保留前 N_path 个峰，升序排列
% end
% 
% ToA_est = scan_axis(idxs);
% powers = p_vec(idxs);
% 
% end
function [ToA_est, powers, p_vec, scan_axis] = estimate_ToA_SAMV3(CIR, N_path, step)
% 使用 SAMV-3 算法估计 CIR 中多个路径的 ToA（支持亚采样精度），使用KMeans聚类消除主峰多重效应
%
% 输入：
%   CIR     ：N×T 矩阵（N 为 CIR 长度，T 为快照数，可为1）
%   N_path  ：希望估计的路径数（例如直射+反射 = 2）
%   step    ：亚采样精度步长（如0.1表示十分之一采样点）
%
% 输出：
%   ToA_est ：1×N_path，亚采样级别的 ToA 索引（单位：采样点）
%   powers  ：每个路径对应的功率
%   p_vec   ：估计出的功率谱向量（length(scan_axis)×1）
%   scan_axis ：扫描的延迟坐标（单位：采样点）

if nargin < 3
    step = 0.1;
end

CIR = double(CIR);
if isvector(CIR)
    CIR = CIR(:);
end

[N, T] = size(CIR);
scan_axis = 0:step:(N-1);
K = length(scan_axis);
A = zeros(N, K);

for k = 1:K
    delay = scan_axis(k);
    A(:,k) = sinc((0:N-1)' - delay);
end

DAS_init = mean(abs(CIR).^2, 2);
DAS_init_interp = A' * DAS_init;

R_N = (CIR * CIR') / T;
sigma = mean(abs(CIR(:)).^2);
p_vec_Old = abs(DAS_init_interp).^2;

threshold = 1e-6;
maxIter = 30;

for iter = 1:maxIter
    R = A * diag(p_vec_Old) * A' + sigma * eye(N);
    Rinv = inv(R);
    RinvA = Rinv * A;

    num = sum(conj(A) .* (Rinv * R_N * RinvA), 1).';
    den = sum(conj(A) .* RinvA, 1).';

    p_vec = p_vec_Old .* (num ./ (den + 1e-12));
    sigma = real(trace(Rinv * Rinv * R_N)) / real(trace(Rinv * Rinv));

    if norm(p_vec - p_vec_Old) / norm(p_vec_Old) < threshold
        break;
    end
    p_vec_Old = p_vec;
end

p_vec = real(p_vec);

% 查找所有局部峰
[pks_all, idxs_all] = findpeaks(p_vec, 'MinPeakProminence', max(p_vec) * 0.01);

if length(idxs_all) < N_path
    warning('检测到的峰数量少于路径数');
    ToA_est = nan(1, N_path);
    powers = nan(1, N_path);
    return;
end

% % KMeans 聚类聚集峰
% toa_candidates = scan_axis(idxs_all)';
% [pks_sorted, idxs_sort] = sort(pks_all, 'descend');
% toa_top = toa_candidates(idxs_sort);
% 
% opts = statset('Display','off');
% [idx_class, C] = kmeans(toa_top(1:min(10,end)), N_path, 'Replicates', 5, 'Options', opts);
% 
% ToA_est = zeros(1, N_path);
% powers = zeros(1, N_path);
% for k = 1:N_path
%     cluster_points = find(idx_class == k);
%     if isempty(cluster_points)
%         ToA_est(k) = nan;
%         powers(k) = nan;
%         continue;
%     end
%     [~, max_local] = max(pks_sorted(cluster_points));
%     idx_real = cluster_points(max_local);
%     ToA_est(k) = toa_top(idx_real);
%     powers(k) = pks_sorted(idx_real);
% end
% 
% [ToA_est, sort_idx] = sort(ToA_est);
% powers = powers(sort_idx);
% 找出前若干个峰的候选 ToA
toa_candidates = scan_axis(idxs_all)';
[pks_sorted, idxs_sort] = sort(pks_all, 'descend');
toa_top = toa_candidates(idxs_sort(1:min(15,end)));  % 取前15个峰

% DBSCAN 聚类（自适应）
epsilon = 5;  % 峰值最大间距阈值（根据实际采样间隔调整）
minpts = 1;   % 至少1个点为簇
idx_class = dbscan(toa_top, epsilon, minpts);

% 保留聚类标签 > 0 的
unique_clusters = unique(idx_class(idx_class > 0));

ToA_est = zeros(1, length(unique_clusters));
powers = zeros(1, length(unique_clusters));

for k = 1:length(unique_clusters)
    cl = unique_clusters(k);
    cluster_points = find(idx_class == cl);
    [~, max_idx] = max(pks_sorted(cluster_points));
    ToA_est(k) = toa_top(cluster_points(max_idx));
    powers(k) = pks_sorted(cluster_points(max_idx));
end

% 排序
[ToA_est, sort_idx] = sort(ToA_est);
powers = powers(sort_idx);


end
