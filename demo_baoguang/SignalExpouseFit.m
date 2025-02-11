function [betaResult,roResult,tcResult,VnoiseResult] = SignalExpouseFit(ave_skData,T_values)
% 单曝光散斑衬比成像，拟合求解（每一个 T 与对应的 ks 拟合出一张 tc 图，15 个 T 和 ks 拟合出 15 个 tc 图）
% 输入参数：
%     ave_skData = 存储散斑衬比数据的结构体，前置使用 Sum_calcSpeckleContrast.m 文件
%     T_values = [1; 2; 3] | 1 （可输入多组数据同时计算，也可输入一组数据算一个）
% 输出参数：
%     [betaResult,roResult,tcResult,VnoiseResult] = 存储拟合出结果的结构体


% 创建结构体数组，用于存储每组数据的拟合参数
betaResult = struct('result', cell(1, numel(T_values)));
roResult = struct('result', cell(1, numel(T_values)));
tcResult = struct('result', cell(1, numel(T_values)));
VnoiseResult = struct('result', cell(1, numel(T_values)));

[row,col]=size(ave_skData(1).average_sk);

% 创建数组来存储结果
betaMatrix = zeros(row ,col);
roMatrix = zeros(row,col);
tcMatrix = zeros(row,col);
VnoiseMatrix = zeros(row,col);

% 创建进度条
h = waitbar(0, 'Fitting Progress...');
for tIndex = 1:numel(T_values)
    T = T_values(tIndex);
    
    % 获取对应 T 值下的 ave_skData 矩阵
    ave_skMatrix = ave_skData(tIndex).average_sk; % 假设结构体中包含名为 average_sk 的矩阵字段

    % 在内部循环中遍历 ave_skData 矩阵
    parfor l = 1:row
        for w = 1:col
            ks_value = ave_skMatrix(l,w) ^ 2;
            model = @(params, T) sqrt((params(1).*(params(2).^2).*((exp(-2*(T./params(3)))-1+2*(T./params(3)))./(2*(T./params(3)).^2)) +...
                4.*params(1).*params(2).*(1-params(2)).*((exp(-(T./params(3)))-1+(T./params(3)))./((T./params(3)).^2)) + ...
                params(1).*(1-params(2)).^2) +...
                params(4));
            
            % 创建 Options 结构体，设置容差值和关闭显示迭代
            options = optimset('TolFun', 1e-6, 'Display', 'off');
            
            % 非线性拟合
            lower_limit = [1e-9, 0, 0, 0]; % 下限，可以根据需要调整
            upper_limit = [Inf, Inf, 1, Inf]; % 上限，可以根据需要调整
            Initial_guess = [1e-6,0.5,0.9,0];
            
            % 拟合
            estimated_params = lsqcurvefit(model, Initial_guess, double(T), double(ks_value), lower_limit, upper_limit, options);

            % 提取拟合后的参数
            beta = estimated_params(1);
            ro = estimated_params(2);
            tc = estimated_params(3);
            Vnoise = estimated_params(4);
            
            % 存储结果到相应的结果矩阵中
            betaMatrix(l,w) = beta;
            roMatrix(l,w) = ro;
            tcMatrix(l,w) = tc;
            VnoiseMatrix(l,w) = Vnoise;
        end
    end
    
    % 存储结果到相应的结构体
    betaResult(tIndex).result = betaMatrix;
    roResult(tIndex).result = roMatrix;
    tcResult(tIndex).result = tcMatrix;
    VnoiseResult(tIndex).result = VnoiseMatrix;
    
    % 更新进度条
    waitbar(tIndex / numel(T_values), h);
end

% 关闭进度条
close(h);
end
