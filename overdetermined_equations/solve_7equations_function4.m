function [best_x, best_y, best_beta, best_a, best_b, best_p, best_q, min_res] = solve_7equations_function4(x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 定义合理的搜索范围和步长
    % 注意：由于有7个变量，我们需要使用非常粗糙的网格
    x_values = -5:0.5:5;        % 发射端x坐标
    y_values = 0:0.5:10;        % 发射端y坐标（假设在地面以上）
    beta_values = -pi:0.3:pi;   % 发射端朝向
    
    % 反射面参数范围（根据实际情况调整）
    a_values = -1:0.2:1;        % 第一个反射面斜率
    b_values = -2:0.5:2;        % 第一个反射面截距
    p_values = -1:0.2:1;        % 第二个反射面斜率
    q_values = -2:0.5:2;        % 第二个反射面截距
    
    % 初始化最小残差和最佳解
    min_res = inf;
    best_x = NaN;
    best_y = NaN;
    best_beta = NaN;
    best_a = NaN;
    best_b = NaN;
    best_p = NaN;
    best_q = NaN;
    
    % 计算总迭代次数（用于显示进度）
    total_iterations = length(x_values) * length(y_values) * length(beta_values) * ...
                      length(a_values) * length(b_values) * length(p_values) * length(q_values);
    current_iteration = 0;
    
    % 循环所有变量
    for i = 1:length(x_values)
        x = x_values(i);
        for j = 1:length(y_values)
            y = y_values(j);
            for k = 1:length(beta_values)
                beta = beta_values(k);
                for l = 1:length(a_values)
                    a = a_values(l);
                    for m = 1:length(b_values)
                        b = b_values(m);
                        for n = 1:length(p_values)
                            p = p_values(n);
                            for o = 1:length(q_values)
                                q = q_values(o);
                                
                                current_iteration = current_iteration + 1;
                                
                                % 每10000次迭代显示一次进度
                                if mod(current_iteration, 100000) == 0
                                    fprintf('进度: %.2f%%\n', 100 * current_iteration / total_iterations);
                                end
                                
                                % 计算所有七个方程的残差
                                [err1, err2, err3, err4, err5, err6, err7] = ...
                                    compute_all_residuals(x, y, beta, a, b, p, q, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2);
                                
                                % 如果任何方程出现数值问题，跳过这个点
                                if any(isnan([err1, err2, err3, err4, err5, err6, err7])) || ...
                                   any(isinf([err1, err2, err3, err4, err5, err6, err7]))
                                    continue;
                                end
                                
                                % 计算总残差（加权平方和）
                                weights = [1.515, 0.658, 1, 0.546, 1, 1, 1];
                                res = weights(1)*err1^2 + weights(2)*err2^2 + weights(3)*err3^2 + ...
                                      weights(4)*err4^2 + weights(5)*err5^2 + weights(6)*err6^2 + ...
                                      weights(7)*err7^2;
                                
                                % 更新最佳解
                                if res < min_res
                                    min_res = res;
                                    best_x = x;
                                    best_y = y;
                                    best_beta = beta;
                                    best_a = a;
                                    best_b = b;
                                    best_p = p;
                                    best_q = q;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function [err1, err2, err3, err4, err5, err6, err7] = compute_all_residuals(x, y, beta, a, b, p, q, x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 初始化所有残差
    err1 = NaN; err2 = NaN; err3 = NaN; err4 = NaN;
    err5 = NaN; err6 = NaN; err7 = NaN;
    
    % 计算方程1的残差：直接路径到达角
    if abs(x - x0) > 1e-10
        err1 = (y - y0) / (x - x0) - tan(alpha + phi1);
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
    if abs(x0 - x_s1) > 1e-10
        err2 = (y0 - y_s1) / (x0 - x_s1) - tan(alpha + phi2);
    else
        err2 = 1e6; % 给予一个大残差值
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
        err5 = (y0 - y_s2) / (x0 - x_s2) - tan(alpha + phi3);
    else
        err5 = 1e6; % 给予一个大残差值
    end
    
    % 计算方程6的残差：第二个反射路径长度差
    reflect_path2 = sqrt((x_s2 - x)^2 + (y_s2 - y)^2) + sqrt((x0 - x_s2)^2 + (y0 - y_s2)^2);
    err6 = reflect_path2 - direct_path - K2;
    
    % 计算方程7的残差：第二个反射路径发射端离开角正弦差
    C_angle = atan2(y_s2 - y, x_s2 - x);
    err7 = sin(A_angle - beta) - sin(C_angle - beta) - D2;
end