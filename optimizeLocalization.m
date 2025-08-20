function xOpt = optimizeLocalization(x0, theta1, theta2, phi1, phi2, rTOF)
    % 转换角度为弧度
    theta1 = deg2rad(theta1);
    theta2 = deg2rad(theta2);
    phi1 = deg2rad(phi1);
    phi2 = deg2rad(phi2);

    % 迭代搜索参数
    searchParams = {
        struct('range', 1, 'step', 0.03) % 中等搜索范围和步长
    };

    % 逐级迭代搜索
    for i = 1:length(searchParams)
        param = searchParams{i};
        fineRange = param.range;
        fineStep = param.step;

        % 动态生成精细的网格点坐标
        x1Range = (x0(1)-fineRange):fineStep:(x0(1)+fineRange);
        y1Range = (x0(2)-fineRange):fineStep:(x0(2)+fineRange);
        x2Range = (x0(3)-fineRange):fineStep:(x0(3)+fineRange);
        y2Range = (x0(4)-fineRange):fineStep:(x0(4)+fineRange);

        [X1, Y1, X2, Y2] = ndgrid(x1Range, y1Range, x2Range, y2Range);
        fineGrid = [X1(:), Y1(:), X2(:), Y2(:)];

        % 计算目标函数值
        fineObjectiveValues = objectiveFunction(fineGrid, theta1, theta2, phi1, phi2, rTOF);

        % 找到最小值
        [~, fineIdx] = min(fineObjectiveValues);
        x0 = fineGrid(fineIdx, :);
        % disp(['Iteration ', num2str(i), ' Best x0: ', mat2str(x0)]);
    end

    xOpt = x0;
end


function fval = objectiveFunction(x, theta1, theta2, phi1, phi2, rTOF)
    x1 = x(:, 1);
    y1 = x(:, 2);
    x2 = x(:, 3);
    y2 = x(:, 4);

    % 计算目标函数值
    fval = (tan(theta1) - x1./y1).^2 + (tan(theta2) - x2./y2).^2 + (tan(phi1 + theta1 - phi2) - (x1 - x2)./(y1 - y2)).^2 + ((rTOF * physconst('Lightspeed') * 1e-9) - (-sqrt(x1.^2 + y1.^2) + sqrt(x2.^2 + y2.^2) + sqrt((x1 - x2).^2 + (y1 - y2).^2))).^2;
end