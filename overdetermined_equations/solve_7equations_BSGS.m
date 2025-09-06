function [best_solution, best_res] = solve_7equations_BSGS(x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 第一层：粗网格搜索
    fprintf('开始第一层（粗网格）搜索...\n');
    
    % 粗网格参数
    x_coarse = -5:0.5:5;          % 发射端x坐标
    y_coarse = -10:0.5:10;          % 发射端y坐标
    beta_coarse = -pi:0.5:pi;   % 发射端朝向
    a_coarse = -1:0.5:1;        % 第一个反射面斜率
    b_coarse = -5:0.5:5;          % 第一个反射面截距
    p_coarse = -1:0.5:1;        % 第二个反射面斜率
    q_coarse = -5:0.5:5;          % 第二个反射面截距
    
    % 初始化最小残差和最佳解
    min_res = inf;
    best_solution = zeros(7, 1);
    
    % 粗网格搜索
    for i = 1:length(x_coarse)
        x = x_coarse(i);
        for j = 1:length(y_coarse)
            y = y_coarse(j);
            for k = 1:length(beta_coarse)
                beta = beta_coarse(k);
                for l = 1:length(a_coarse)
                    a = a_coarse(l);
                    for m = 1:length(b_coarse)
                        b = b_coarse(m);
                        for n = 1:length(p_coarse)
                            p = p_coarse(n);
                            for o = 1:length(q_coarse)
                                q = q_coarse(o);
                                
                                % 计算残差
                                res = compute_residuals([x; y; beta; a; b; p; q], x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2);
                                
                                % 更新最佳解
                                if res < min_res
                                    min_res = res;
                                    best_solution = [x; y; beta; a; b; p; q];
                                    fprintf('新最佳残差: %.6f\n', min_res);
                                    fprintf("There are solutions:x:%d y:%d beta:%d a:%d b:%d p%d q%d",x,y,beta,a,b,p,q);
                                end
                            end
                        end
                    end
                end
            end
        end
        fprintf('x进度: %.1f%%\n', 100*i/length(x_coarse));
    end
    
    fprintf('第一层搜索完成。最佳残差: %.6f\n', min_res);
    fprintf('最佳解: x=%.3f, y=%.3f, beta=%.3f, a=%.3f, b=%.3f, p=%.3f, q=%.3f\n', ...
            best_solution(1), best_solution(2), best_solution(3), ...
            best_solution(4), best_solution(5), best_solution(6), best_solution(7));
    
    % 第二层：细网格搜索（在最佳解附近）
    fprintf('开始第二层（细网格）搜索...\n');
    
    % 定义细网格范围和步长
    search_range = 1.0;  % 搜索范围（围绕最佳解）
    step_size = 0.1;     % 步长
    
    x_fine = best_solution(1) - search_range : step_size : best_solution(1) + search_range;
    y_fine = best_solution(2) - search_range : step_size : best_solution(2) + search_range;
    beta_fine = best_solution(3) - 0.5 : 0.1 : best_solution(3) + 0.5;
    a_fine = best_solution(4) - 0.2 : 0.05 : best_solution(4) + 0.2;
    b_fine = best_solution(5) - 0.5 : 0.1 : best_solution(5) + 0.5;
    p_fine = best_solution(6) - 0.2 : 0.05 : best_solution(6) + 0.2;
    q_fine = best_solution(7) - 0.5 : 0.1 : best_solution(7) + 0.5;
    
    % 确保细网格在合理范围内
    x_fine = max(-5, min(5, x_fine));
    y_fine = max(0, min(10, y_fine));
    beta_fine = max(-pi, min(pi, beta_fine));
    a_fine = max(-1, min(1, a_fine));
    b_fine = max(-2, min(2, b_fine));
    p_fine = max(-1, min(1, p_fine));
    q_fine = max(-2, min(2, q_fine));
    
    % 细网格搜索
    for i = 1:length(x_fine)
        x = x_fine(i);
        for j = 1:length(y_fine)
            y = y_fine(j);
            for k = 1:length(beta_fine)
                beta = beta_fine(k);
                for l = 1:length(a_fine)
                    a = a_fine(l);
                    for m = 1:length(b_fine)
                        b = b_fine(m);
                        for n = 1:length(p_fine)
                            p = p_fine(n);
                            for o = 1:length(q_fine)
                                q = q_fine(o);
                                
                                % 计算残差
                                res = compute_residuals([x; y; beta; a; b; p; q], x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2);
                                if res == inf
                                    continue;
                                end
                                % 更新最佳解
                                if res < min_res
                                    min_res = res;
                                    best_solution = [x; y; beta; a; b; p; q];
                                    fprintf('新最佳残差: %.6f\n', min_res);
                                end
                            end
                        end
                    end
                end
            end
        end
        fprintf('细网格x进度: %.1f%%\n', 100*i/length(x_fine));
    end
    
    fprintf('分层搜索完成。最终残差: %.6f\n', min_res);
    fprintf('最终解: x=%.3f, y=%.3f, beta=%.3f, a=%.3f, b=%.3f, p=%.3f, q=%.3f\n', ...
            best_solution(1), best_solution(2), best_solution(3), ...
            best_solution(4), best_solution(5), best_solution(6), best_solution(7));
    
    best_res = min_res;
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