function [best_x, best_y, best_beta, min_res] = find_solution(x0, y0, alpha, phi1, phi2, K, D)
    % 定义网格：x和y从0到10米，步长0.1米
    x_values = -10:0.01:10;
    y_values = -10:0.01:10;
    
    % 初始化最小残差和最佳解
    min_res = inf;
    best_x = NaN;
    best_y = NaN;
    best_beta = NaN;
    
    % 循环所有x和y值
    for i = 1:length(x_values)
        x = x_values(i);
        for j = 1:length(y_values)
            y = y_values(j);
            
            % 跳过x等于x0的情况（避免除以零）
            if x == x0
                continue;
            end
            
            % 计算方程1的残差：直接路径到达角
            err1 = (y - y0) / (x - x0) - tan(alpha + phi1);
            
            % 计算方程2的残差：反射路径到达角
            err2 = (-(y + y0)) / (x - x0) - tan(alpha + phi2);
            
            % 计算方程3的残差：路径长度差
            term1 = sqrt((x - x0)^2 + (y + y0)^2);
            term2 = sqrt((x - x0)^2 + (y - y0)^2);
            err3 = term1 - term2 - K;
            
            % 计算方程4的残差和最佳beta
            A = atan2(y0 - y, x0 - x); % 计算A
            if y + y0 == 0
                % 处理y+y0=0的情况（分母为零），设置大残差
                err4 = 1e9;
                beta_candidate = 0;
            else
                B = atan2(-y, (x0 - x) * y / (y + y0)); % 计算B
                C = A - B;
                sin_halfC = sin(C / 2);
                if abs(sin_halfC) < 1e-8 % 处理sin(C/2)≈0的情况
                    err4 = 1e9;
                    beta_candidate = 0;
                else
                    value = D / (2 * sin_halfC);
                    if abs(value) > 1 % 处理无解情况
                        err4 = 1e9;
                        beta_candidate = 0;
                    else
                        % 计算两个可能的beta值
                        acos_value = acos(value);
                        beta1 = (A + B)/2 - acos_value;
                        beta2 = (A + B)/2 + acos_value;
                        % 选择使方程4残差较小的beta
                        err4_1 = sin(A - beta1) - sin(B - beta1) - D;
                        err4_2 = sin(A - beta2) - sin(B - beta2) - D;
                        if abs(err4_1) < abs(err4_2)
                            err4 = err4_1;
                            beta_candidate = beta1;
                        else
                            err4 = err4_2;
                            beta_candidate = beta2;
                        end
                    end
                end
            end
            
            % 计算总残差（平方和）
            %weights = [1.75,2.94,1,1.83];
            %res = 1.75*err1^2 + 2.94*err2^2 + err3^2 + 1.83*err4^2;
            res = 1.515*err1^2 + 0.658*err2^2 + err3^2 + 0.546*err4^2;
            %res = err1^2 + err2^2 + err3^2 + err4^2;
            if best_x == -4 && best_y == 2
                fprintf(res)
            end
            % 更新最佳解
            if res < min_res
                min_res = res;
                best_x = x;
                best_y = y;
                best_beta = beta_candidate;
            end
        end
    end
end