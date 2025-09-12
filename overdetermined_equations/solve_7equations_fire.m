function [best_solution, best_res] = solve_7equations_fire(x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 参数设置
    initial_temp = 10000000;    % 初始温度
    cooling_rate = 0.97;    % 冷却速率
    num_iterations = 300000;  % 每个温度下的迭代次数,注意控制迭代次数
    min_temp = 1e-6;        % 最小温度
    rng('shuffle');
    % 定义变量的上下界
    bounds = [
        -5, 5;     % x
        -10, 10;      % y
        -pi, pi;    % beta
        -2, 2;      % a
        -5, 5;      % b
        -2, 2;      % p
        -5, 5       % q
    ];
    
    % 随机初始化当前解
    current_solution = zeros(7, 1);
    for i = 1:7
        current_solution(i) = bounds(i, 1) + (bounds(i, 2) - bounds(i, 1)) * rand();
    end
    
    % 计算当前解的残差
    current_res = compute_residuals(current_solution, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2);
    best_solution = current_solution;
    best_res = current_res;
    
    % 初始化温度
    temperature = initial_temp;
    
    % 模拟退火主循环
    while temperature > min_temp
        for i = 1:num_iterations
            % 生成新解
            new_solution = current_solution;
            for j = 1:7
                % 在当前解附近随机扰动
                perturbation = (bounds(j, 2) - bounds(j, 1)) * 0.1 * (2 * rand() - 1);
                new_solution(j) = new_solution(j) + perturbation;
                % 确保新解在边界内
                new_solution(j) = max(bounds(j, 1), min(bounds(j, 2), new_solution(j)));
            end
            
            % 计算新解的残差
            new_res = compute_residuals(new_solution, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2);
            if new_res == inf
                continue;
            end
            % 计算能量差
            delta_energy = new_res - current_res;
            
            % 决定是否接受新解
            if delta_energy < 0 
                % 新解更好，总是接受
                current_solution = new_solution;
                current_res = new_res;
               
                % 更新最佳解
                if new_res < best_res 
                    best_solution = new_solution;
                    best_res = new_res;
                end
            else
                % 新解更差，以一定概率接受（避免陷入局部最优）
                acceptance_prob = exp(-delta_energy / temperature);
                if rand() < acceptance_prob
                    current_solution = new_solution;
                    current_res = new_res;
                end
            end
        end
        
        % 降低温度
        temperature = temperature * cooling_rate;
        
        % 显示进度
        fprintf('温度: %.6f, 当前残差: %.6f, 最佳残差: %.6f\n', temperature, current_res, best_res);
        
    end
    
    % 返回最佳解
    fprintf('模拟退火完成。最佳残差: %.6f\n', best_res);
end

function res = compute_residuals(params, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 解包参数
    x = params(1);
    y = params(2);
    beta = params(3);
    a = params(4);
    b = params(5);
    p = params(6);
    q = params(7);

    % 计算所有残差
    [err1, err2, err3, err4, err5, err6, err7] = ...
        compute_all_residuals(x, y, beta, a, b, p, q, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2);
    
    % 如果有任何方程出现数值问题，返回一个大残差值
    if any(isnan([err1, err2, err3, err4, err5, err6, err7])) || ...
       any(isinf([err1, err2, err3, err4, err5, err6, err7]))
        res = inf;
        return;
    end
    
    
    % 计算总残差（加权平方和）
    % 如果这里alpha+phi2是无穷大，直接把这个方程弄成0如何
    weights = [1.515, 1, 1, 0.546, 10, 10, 10];
    res = weights(1)*err1^2 + weights(2)*err2^2 + weights(3)*err3^2 + ...
          weights(4)*err4^2 + weights(5)*err5^2 + weights(6)*err6^2 + ...
          weights(7)*err7^2;
   % fprintf("There are errors %d %d %d %d %d %d %d\n",err1,err2,err3,err4,err5,err6,err7);
end
function [err1, err2, err3, err4, err5, err6, err7] = compute_all_residuals(x, y, beta, a, b, p, q, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 初始化所有残差
    err1 = NaN; err2 = NaN; err3 = NaN; err4 = NaN;
    err5 = NaN; err6 = NaN; err7 = NaN;
    
    % 计算方程1的残差：直接路径到达角
    if abs(x - x0) > 1e-10
        %err1 = (y - y0) / (x - x0) - tan(alpha + phi1);
        err1 = (y-y0)*cos(alpha+phi1) - (x-x0)*sin(alpha+phi1);
    else
        % 处理x = x0的情况
        err1 = 1e6; % 给予一个大残差值
    end
    
    % 计算第一个反射点S1
    k1 = (y0 - a*x0 - b) / (1 + a^2);
    x_m1 = x0 + 2*a*k1;
    y_m1 = y0 - 2*k1;
    
    denominator = (y_m1 - y) - a*(x_m1 - x);
    if abs(denominator) < 1e-10
        return; % 跳过分母为零的情况
    end
    t = (a*x + b - y) / denominator;
    x_s1 = x + t*(x_m1 - x);
    y_s1 = y + t*(y_m1 - y);
    
    % 计算方程2的残差：第一个反射路径到达角
    %fprintf("There are something:%d %d %d",y0-y_s1,x0-x_s1,tan(alpha+phi2))
    if abs(x0 - x_s1) > 1e-10
        %fprintf("坐标:%d %d %d %d %d\n",y0,y_s1,x0,x_s1,tan(alpha+phi2));
        tmp = tan(alpha + phi2);
        if(abs(alpha + phi2) == pi/2)
            tmp = 0;
        end
        %err2 = (y0 - y_s1) / (x0 - x_s1) - tmp;
        err2 = (y0-y_s1)*cos(alpha+phi2) - (x0-x_s1)*sin(alpha+phi2);
    else
        err2 = inf; % 给予一个大残差值
        return;
    end
    
    % 计算方程3的残差：第一个反射路径长度差
    direct_path = sqrt((x0 - x)^2 + (y0 - y)^2);
    reflect_path1 = sqrt((x_s1 - x)^2 + (y_s1 - y)^2) + sqrt((x0 - x_s1)^2 + (y0 - y_s1)^2);
    err3 = reflect_path1 - direct_path - K1;
    
    % 计算方程4的残差：第一个反射路径发射端离开角正弦差
    A_angle = atan2(y0 - y, x0 - x);
    B_angle = atan2(y_s1 - y, x_s1 - x);
    err4 = sin(A_angle - beta) - sin(B_angle - beta) - D1;
    
    % 计算第二个反射点S2
    k2 = (y0 - p*x0 - q) / (1 + p^2);
    x_m2 = x0 + 2*p*k2;
    y_m2 = y0 - 2*k2;
    
    denominator2 = (y_m2 - y) - p*(x_m2 - x);
    if abs(denominator2) < 1e-10
        return; % 跳过分母为零的情况
    end
    s = (p*x + q - y) / denominator2;
    x_s2 = x + s*(x_m2 - x);
    y_s2 = y + s*(y_m2 - y);
    
    % 计算方程5的残差：第二个反射路径到达角
    if abs(x0 - x_s2) > 1e-10
       % err5 = (y0 - y_s2) / (x0 - x_s2) - tan(alpha + phi3);
        err5 = (y0 - y_s2)*cos(alpha+phi3) - (x0 - x_s2)*sin(alpha + phi3);
    else
        err5 = inf; % 给予一个大残差值
    end
    %问题在于tan(alpha+phi3)这里直接干到tan(90°)了，所以直接炸了！
    
    % 计算方程6的残差：第二个反射路径长度差
    reflect_path2 = sqrt((x_s2 - x)^2 + (y_s2 - y)^2) + sqrt((x0 - x_s2)^2 + (y0 - y_s2)^2);
    err6 = reflect_path2 - direct_path - K2;
    
    % 计算方程7的残差：第二个反射路径发射端离开角正弦差
    C_angle = atan2(y_s2 - y, x_s2 - x);
    err7 = sin(A_angle - beta) - sin(C_angle - beta) - D2;
end