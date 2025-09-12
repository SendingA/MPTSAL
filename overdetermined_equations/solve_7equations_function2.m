function [best_x, best_y, min_res] = solve_7equations_function2(x0, y0, alpha, phi1, phi2, phi3, D1, D2, K1, K2)
    % 定义网格范围：只搜索x和y
    x_values = -10:0.1:10;
    y_values = -10:0.1:10;
    
    % 初始化最小残差和最佳解
    min_res = inf;
    best_x = NaN;
    best_y = NaN;
    
    % 循环所有x和y值
    for i = 1:length(x_values)
        x = x_values(i);
        for j = 1:length(y_values)
            y = y_values(j);
            
            % 跳过x等于x0的情况（避免除以零）
            if abs(x - x0) < 1e-10
                continue;
            end
            
            % 计算方程1的残差：直接路径到达角
            err1 = (y - y0) / (x - x0) - tan(alpha + phi1);
            
            % 使用几何关系估计第一个反射面参数
            % 假设反射点位于发射端和接收端的垂直平分线上,该假设正确
            mid_x = (x + x0) / 2;
            mid_y = (y + y0) / 2;
            
            % 计算第一个反射面的法线方向
            if abs(y - y0) > 1e-10
                a1 = -(x - x0) / (y - y0); % 法线斜率
                % 反射面方程: y = a1*x + b1
                b1 = mid_y - a1 * mid_x;
            else
                % 处理水平情况
                a1 = 0;
                b1 = mid_y;
            end
            
            % 计算第一个反射点S1
            k1 = (y0 - a1*x0 - b1) / (1 + a1^2);
            x_m1 = x0 + 2*a1*k1;
            y_m1 = y0 - 2*k1;
            
            % 计算参数t
            denominator = (y_m1 - y) - a1*(x_m1 - x);
            if abs(denominator) < 1e-10
                continue; % 跳过分母为零的情况
            end
            t = (a1*x + b1 - y) / denominator;
            x_s1 = x + t*(x_m1 - x);
            y_s1 = y + t*(y_m1 - y);
            
            % 计算方程2的残差：第一个反射路径到达角
            if abs(x0 - x_s1) > 1e-10
                err2 = (y0 - y_s1) / (x0 - x_s1) - tan(alpha + phi2);
            else
                err2 = 0; % 特殊处理
            end
            
            % 计算方程3的残差：第一个反射路径长度差
            direct_path = sqrt((x0 - x)^2 + (y0 - y)^2);
            reflect_path1 = sqrt((x_s1 - x)^2 + (y_s1 - y)^2) + sqrt((x0 - x_s1)^2 + (y0 - y_s1)^2);
            err3 = reflect_path1 - direct_path - K1;
            
            % 使用几何关系估计第二个反射面参数
            % 假设第二个反射面与第一个反射面垂直
            a2 = -1/a1; % 垂直面的斜率
            b2 = mid_y - a2 * mid_x;
            
            % 计算第二个反射点S2
            k2 = (y0 - a2*x0 - b2) / (1 + a2^2);
            x_m2 = x0 + 2*a2*k2;
            y_m2 = y0 - 2*k2;
            
            % 计算参数s
            denominator2 = (y_m2 - y) - a2*(x_m2 - x);
            if abs(denominator2) < 1e-10
                continue; % 跳过分母为零的情况
            end
            s = (a2*x + b2 - y) / denominator2;
            x_s2 = x + s*(x_m2 - x);
            y_s2 = y + s*(y_m2 - y);
            
            % 计算方程5的残差：第二个反射路径到达角
            if abs(x0 - x_s2) > 1e-10
                err5 = (y0 - y_s2) / (x0 - x_s2) - tan(alpha + phi3);
            else
                err5 = 0; % 特殊处理
            end
            
            % 计算方程6的残差：第二个反射路径长度差
            reflect_path2 = sqrt((x_s2 - x)^2 + (y_s2 - y)^2) + sqrt((x0 - x_s2)^2 + (y0 - y_s2)^2);
            err6 = reflect_path2 - direct_path - K2;
            
            % 计算总残差（加权平方和）
            % 给直接路径方程和路径差方程更高的权重
            weights = [2, 1, 2, 1, 1, 2];
            res = weights(1)*err1^2 + weights(2)*err2^2 + weights(3)*err3^2 + ...
                  weights(4)*err5^2 + weights(5)*err6^2;
            
            % 更新最佳解
            fprintf(res,minres,x,y);
            if res < min_res
                min_res = res;
                best_x = x
                best_y = y
            end
        end
        
        % 显示进度
        if mod(i, 10) == 0
            fprintf('已完成 %.1f%%\n', 100*i/length(x_values));
        end
    end
    
    % 使用优化算法进一步细化结果
    if ~isnan(best_x) && ~isnan(best_y)
        options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6);
        [refined, fval] = fminsearch(@(xy) calculate_residuals(xy, x0, y0, alpha, phi1, phi2, phi3, K1, K2), [best_x, best_y], options);
        best_x = refined(1);
        best_y = refined(2);
        min_res = fval;
    end
end

function res = calculate_residuals(xy, x0, y0, alpha, phi1, phi2, phi3, K1, K2)
    x = xy(1);
    y = xy(2);
    
    % 计算方程1的残差：直接路径到达角
    if abs(x - x0) > 1e-10
        err1 = (y - y0) / (x - x0) - tan(alpha + phi1);
    else
        err1 = 0;
    end
    
    % 使用几何关系估计反射面参数
    mid_x = (x + x0) / 2;
    mid_y = (y + y0) / 2;
    
    if abs(y - y0) > 1e-10
        a1 = -(x - x0) / (y - y0);
        b1 = mid_y - a1 * mid_x;
    else
        a1 = 0;
        b1 = mid_y;
    end
    
    % 计算第一个反射点S1
    k1 = (y0 - a1*x0 - b1) / (1 + a1^2);
    x_m1 = x0 + 2*a1*k1;
    y_m1 = y0 - 2*k1;
    
    denominator = (y_m1 - y) - a1*(x_m1 - x);
    if abs(denominator) < 1e-10
        res = inf;
        return;
    end
    t = (a1*x + b1 - y) / denominator;
    x_s1 = x + t*(x_m1 - x);
    y_s1 = y + t*(y_m1 - y);
    
    % 计算方程2的残差：第一个反射路径到达角
    if abs(x0 - x_s1) > 1e-10
        err2 = (y0 - y_s1) / (x0 - x_s1) - tan(alpha + phi2);
    else
        err2 = 0;
    end
    
    % 计算方程3的残差：第一个反射路径长度差
    direct_path = sqrt((x0 - x)^2 + (y0 - y)^2);
    reflect_path1 = sqrt((x_s1 - x)^2 + (y_s1 - y)^2) + sqrt((x0 - x_s1)^2 + (y0 - y_s1)^2);
    err3 = reflect_path1 - direct_path - K1;
    
    % 计算第二个反射面参数
    a2 = -1/a1;
    b2 = mid_y - a2 * mid_x;
    
    % 计算第二个反射点S2
    k2 = (y0 - a2*x0 - b2) / (1 + a2^2);
    x_m2 = x0 + 2*a2*k2;
    y_m2 = y0 - 2*k2;
    
    denominator2 = (y_m2 - y) - a2*(x_m2 - x);
    if abs(denominator2) < 1e-10
        res = inf;
        return;
    end
    s = (a2*x + b2 - y) / denominator2;
    x_s2 = x + s*(x_m2 - x);
    y_s2 = y + s*(y_m2 - y);
    
    % 计算方程5的残差：第二个反射路径到达角
    if abs(x0 - x_s2) > 1e-10
        err5 = (y0 - y_s2) / (x0 - x_s2) - tan(alpha + phi3);
    else
        err5 = 0;
    end
    
    % 计算方程6的残差：第二个反射路径长度差
    reflect_path2 = sqrt((x_s2 - x)^2 + (y_s2 - y)^2) + sqrt((x0 - x_s2)^2 + (y0 - y_s2)^2);
    err6 = reflect_path2 - direct_path - K2;
    
    % 计算总残差
    weights = [2, 1, 2, 1, 2];
    res = weights(1)*err1^2 + weights(2)*err2^2 + weights(3)*err3^2 + ...
          weights(4)*err5^2 + weights(5)*err6^2;
end