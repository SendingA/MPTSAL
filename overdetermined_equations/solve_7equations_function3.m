function [best_x, best_y, best_beta, best_a, best_b, best_p, best_q, min_res] = solve_7equations_function3(x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 定义网格范围：x, y, beta
    x_values = -10:0.05:10;  % 步长增大以减少计算量
    y_values = -10:0.05:10;
    beta_values = -pi:0.2:pi; %确实2pi的方向均搜索过了
    % 耐心看完所有的内容
    % 初始化最小残差和最佳解
    min_res = inf;
    best_x = NaN;
    best_y = NaN;
    best_beta = NaN;
    best_a = NaN;
    best_b = NaN;
    best_p = NaN;
    best_q = NaN;
    % 需要把当前计算的所有内容和
    % 循环所有x, y和beta值
    for i = 1:length(x_values)
        x = x_values(i);
        for j = 1:length(y_values)
            y = y_values(j);
            
            % 跳过x等于x0的情况（避免除以零）
            if abs(x - x0) < 1e-10
                continue;
            end
            
            for k = 1:length(beta_values)
                beta = beta_values(k);
                
                % 计算方程1的残差：直接路径到达角
                err1 = (y - y0) / (x - x0) - tan(alpha + phi1);
                
                % 使用方程1和方程2计算第一个反射面参数a和b
                % 从方程1和方程2可以推导出a和b的表达式
                % 这里使用最小二乘法求解a和b
                A = [x0, 1; x, 1];
                b_vec = [y0 - tan(alpha + phi2)*(x0 - x); y];
                ab = A \ b_vec;
                a = ab(1);
                b = ab(2); %这种方法假设反射面必须经过发射端了，大错特错
                %y = ax + b大概不适用于这种方案
                
                % 计算第一个反射点S1
                k1 = (y0 - a*x0 - b) / (1 + a^2);
                x_m1 = x0 + 2*a*k1;
                y_m1 = y0 - 2*k1;
                
                % 计算参数t
                denominator = (y_m1 - y) - a*(x_m1 - x);
                if abs(denominator) < 1e-10
                    continue; % 跳过分母为零的情况
                end
                t = (a*x + b - y) / denominator;
                x_s1 = x + t*(x_m1 - x);
                y_s1 = y + t*(y_m1 - y);
                
                % 计算方程2的残差：第一个反射路径到达角
                err2 = (y0 - y_s1) / (x0 - x_s1) - tan(alpha + phi2);
                
                % 计算方程3的残差：第一个反射路径长度差
                direct_path = sqrt((x0 - x)^2 + (y0 - y)^2);
                reflect_path1 = sqrt((x_s1 - x)^2 + (y_s1 - y)^2) + sqrt((x0 - x_s1)^2 + (y0 - y_s1)^2);
                err3 = reflect_path1 - direct_path - K1;
                
                % 计算方程4的残差：第一个反射路径发射端离开角正弦差
                A_angle = atan2(y0 - y, x0 - x);
                B_angle = atan2(y_s1 - y, x_s1 - x);
                err4 = sin(A_angle - beta) - sin(B_angle - beta) - D1;
                
                % 使用方程5计算第二个反射面参数p和q
                % 从方程5可以推导出p和q的表达式
                % 这里使用最小二乘法求解p和q
                A2 = [x0, 1; x, 1];
                b_vec2 = [y0 - tan(alpha + phi3)*(x0 - x); y];
                pq = A2 \ b_vec2;
                p = pq(1);
                q = pq(2);
                
                % 计算第二个反射点S2
                k2 = (y0 - p*x0 - q) / (1 + p^2);
                x_m2 = x0 + 2*p*k2;
                y_m2 = y0 - 2*k2;
                
                % 计算参数s
                denominator2 = (y_m2 - y) - p*(x_m2 - x);
                if abs(denominator2) < 1e-10
                    continue; % 跳过分母为零的情况
                end
                s = (p*x + q - y) / denominator2;
                x_s2 = x + s*(x_m2 - x);
                y_s2 = y + s*(y_m2 - y);
                
                % 计算方程5的残差：第二个反射路径到达角
                err5 = (y0 - y_s2) / (x0 - x_s2) - tan(alpha + phi3);
                
                % 计算方程6的残差：第二个反射路径长度差
                reflect_path2 = sqrt((x_s2 - x)^2 + (y_s2 - y)^2) + sqrt((x0 - x_s2)^2 + (y0 - y_s2)^2);
                err6 = reflect_path2 - direct_path - K2;
                
                % 计算方程7的残差：第二个反射路径发射端离开角正弦差
                C_angle = atan2(y_s2 - y, x_s2 - x);
                err7 = sin(A_angle - beta) - sin(C_angle - beta) - D2;
                
                % 计算总残差（加权平方和）
                weights = [1.515, 0.658, 1, 0.546, 1, 1, 1]; % 可以根据需要调整权重
                res = weights(1)*err1^2 + weights(2)*err2^2 + weights(3)*err3^2 + ...
                      weights(4)*err4^2 + weights(5)*err5^2 + weights(6)*err6^2 + ...
                      weights(7)*err7^2;
                if x == 0 && y == 0
                    fprintf("当前res = %d",res)
                end
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
        
        % 显示进度
        if mod(i, 5) == 0
            fprintf('已完成 %.1f%%\n', 100*i/length(x_values));
        end
    end
end