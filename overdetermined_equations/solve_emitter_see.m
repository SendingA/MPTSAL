function [x, y, beta] = solve_emitter_see(x0, y0, alpha, phi1, phi2, D, K, weights)
% 解算发射器位置和方向的超定方程组
% 输入参数修正说明:
%   x0, y0: 接收器位置（根据图片应为(0,2)）
%   alpha: 接收器朝向角度(弧度)
%   phi1, phi2: 到达角(弧度) - 注意符号：φ1为正, φ2为负
%   D: sinθ1 - sinθ2
%   K: 路径长度差 b+c-a
%   weights: 各方程的权重系数 [w1, w2, w3, w4]

% 设置默认权重(如果未提供)
if nargin < 8
    weights = [1, 1, 1, 1];
end

% 根据新约定调整角度：法线朝左为正，需要将角度转换为标准数学约定
alpha_std = alpha + pi; % 将新约定转换为标准约定
phi1_std = phi1; % 到达角保持不变
phi2_std = phi2;

% 计算线性化后的系数矩阵A和常数向量b
A = zeros(4, 3);
b = zeros(4, 1);

% 方程1: 来自phi1的方程（使用调整后的角度）
A(1, 1) = -tan(alpha_std + phi1_std);
A(1, 2) = 1;
b(1) = y0 - tan(alpha_std + phi1_std)*x0;

% 方程2: 来自phi2的方程（使用调整后的角度）
A(2, 1) = -tan(alpha_std + phi2_std);
A(2, 2) = -1;
b(2) = -y0 - tan(alpha_std + phi2_std)*x0;

% 方程3: 路径长度差(需要线性化近似)
[x_init, y_init] = initial_estimate(x0, y0, alpha_std, phi1_std, phi2_std);
dx = x_init - x0;
dy_direct = y_init - y0;
dy_reflect = y_init + y0;

% 计算雅可比矩阵项
d_direct = sqrt(dx^2 + dy_direct^2);
d_reflect = sqrt(dx^2 + dy_reflect^2);

A(3, 1) = (dx/d_reflect) - (dx/d_direct);
A(3, 2) = (dy_reflect/d_reflect) - (dy_direct/d_direct);
b(3) = K + d_direct - d_reflect + A(3,1)*x_init + A(3,2)*y_init;

% 方程4: 离开角正弦差(需要线性化近似)
theta1_global = atan2(y0 - y_init, x0 - x_init);
theta2_global = atan2(-y_init, (x0 - x_init)*y_init/(y_init + y0));

% 计算雅可比矩阵项
A(4, 1) = (y0 - y_init)/(d_direct^2) - (y_init*(y0 + y_init))/((x0 - x_init)^2*y_init^2 + y_init^2*(y0 + y_init)^2);
A(4, 2) = -(x0 - x_init)/(d_direct^2) + ((x0 - x_init)*y0*(y0 + y_init))/((x0 - x_init)^2*y_init^2 + y_init^2*(y0 + y_init)^2);
A(4, 3) = -cos(theta1_global - 0) + cos(theta2_global - 0); % 初始beta=0
b(4) = D - (sin(theta1_global - 0) - sin(theta2_global - 0)) + A(4,1)*x_init + A(4,2)*y_init + A(4,3)*0;

% 应用权重
W = diag(sqrt(weights));
A_weighted = W * A;
b_weighted = W * b;

% 使用SVD求解代替直接求逆，提高数值稳定性
[U, S, V] = svd(A_weighted, 'econ');
s = diag(S);
% 设置奇异值阈值，避免除以极小值
threshold = max(size(A_weighted)) * eps(norm(A_weighted));
s_inv = zeros(size(s));
s_inv(s > threshold) = 1./s(s > threshold);
solution = V * diag(s_inv) * U' * b_weighted;

% 提取结果并转换回新约定
x = solution(1);
y = solution(2);
beta_std = solution(3);
beta = beta_std - pi; % 将标准约定转换回新约定
end

function [x_init, y_init] = initial_estimate(x0, y0, alpha_std, phi1_std, phi2_std)
% 基于前两个方程提供初始估计
m1 = tan(alpha_std + phi1_std);
m2 = tan(alpha_std + phi2_std);

if abs(m1 - m2) < 1e-10
    x_init = x0;
    y_init = y0;
else
    x_init = x0 - (2*y0)/(m1 + m2);
    y_init = (y0*(m2 - m1))/(m1 + m2);
end
end