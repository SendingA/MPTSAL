function [x, y, beta, residual] = solve_equations(x0, y0, alpha, phi1, phi2, K, D, weights)
    % 求解超定方程：最小化加权残差平方和
    % 输入参数:
    %   x0, y0: 发射端坐标（距离单位）
    %   alpha: 接收端朝向角（弧度）
    %   phi1: 直接路径到达角（弧度）
    %   phi2: 反射路径到达角（弧度）
    %   K: 路径长度差（距离单位）
    %   D: 发射端离开角正弦差（无量纲）
    %   weights: 可选，残差权重向量 [w1, w2, w3, w4]，默认全1
    % 输出:
    %   x, y, beta: 求解得到的未知变量（x,y为距离，beta为弧度）
    %   residual: 最终残差向量，用于评估拟合质量

    % 处理可选权重参数
    if nargin < 8
        weights = [1, 1, 1, 1]; % 默认权重
    end

    % 初始猜测：避免分母为零，根据物理意义选择
    initial_guess = [x0 + 1, y0 + 1, 0]; % [x, y, beta]

    % 定义残差函数（内部加权）
    function res = residual_function(sol)
        x_val = sol(1);
        y_val = sol(2);
        beta_val = sol(3);
        
        % 计算各方程残差
        r1 = (y_val - y0) / (x_val - x0) - tan(alpha + phi1);
        r2 = (-(y_val + y0)) / (x_val - x0) - tan(alpha + phi2);
        r3 = sqrt((x_val - x0)^2 + (y_val + y0)^2) - sqrt((x_val - x0)^2 + (y_val - y0)^2) - K;
        
        % 处理第四个方程中的arctan2和分母
        term1 = atan2(y0 - y_val, x0 - x_val) - beta_val;
        term2_denom = y_val + y0;
        if abs(term2_denom) < 1e-8 % 避免除以零
            term2 = atan2(-y_val, sign(x0 - x_val) * 1e8); % 近似处理
        else
            term2 = atan2(-y_val, (x0 - x_val) * y_val / term2_denom) - beta_val;
        end
        r4 = sin(term1) - sin(term2) - D;
        
        % 应用权重：返回加权残差向量
        res = [sqrt(weights(1)) * r1, sqrt(weights(2)) * r2, sqrt(weights(3)) * r3, sqrt(weights(4)) * r4];
    end

    % 使用lsqnonlin求解最小二乘问题
    options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
    [sol, ~, residual] = lsqnonlin(@residual_function, initial_guess, [], [], options);
    
    % 提取解
    x = sol(1);
    y = sol(2);
    beta = sol(3);
end