function [betaMatrix,roMatrix,tcMatrix,VnoiseMatrix] = MutiExpouse(ave_skData,T_values)
% 多曝光散斑衬比成像
%   此处显示详细说明
% 给定 T 值和 Ave_sk 结构体
% T_values = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 12; 14; 16; 18; 20];

[L,W] = size(ave_skData(1).average_sk);

% 创建数组来存储结果
betaMatrix = zeros(L,W);
roMatrix = zeros(L,W);
tcMatrix = zeros(L,W);
VnoiseMatrix = zeros(L,W);

parpool('local', 18);
parfor row = 1:L
    for col = 1:W
        Ks_value = zeros(numel(T_values),1);
        
        % 提取每个曝光时间下 ave_sk 在 (row,col) 处的值
        for i = 1:numel(T_values)
            ave_skMatrix = ave_skData(i).average_sk; % 获取第i个矩阵
            Ks_value(i) = ave_skMatrix(row,col); % 获取每个矩阵的 (row,col) 位置的值
        end
        
        model = @(params, T) sqrt((params(1).*(params(2).^2).*((exp(-2*(T./params(3)))-1+2*(T./params(3)))./(2*(T./params(3)).^2)) +...
            4.*params(1).*params(2).*(1-params(2)).*((exp(-(T./params(3)))-1+(T./params(3)))./((T./params(3)).^2)) + ...
            params(1).*(1-params(2)).^2) +...
            params(4))
        
        % 创建 Options 结构体，设置容差值和关闭显示迭代
        options = optimset('TolFun', 1e-6, 'Display', 'off');
        
        % 非线性拟合
        lower_limit = [1e-9, 0, 0, 0]; % 下限，可以根据需要调整
        upper_limit = [Inf, Inf, 1, Inf]; % 上限，可以根据需要调整
        Initial_guess = [1e-6,0.5,0.9,0];
        
        % 拟合
        estimated_params = lsqcurvefit(model, Initial_guess, double(T_values), double(Ks_value), lower_limit, upper_limit, options);

        % 提取拟合后的参数
        beta = estimated_params(1);
        ro = estimated_params(2);
        tc = estimated_params(3);
        Vnoise = estimated_params(4);
        
        % 存储结果到相应的结果矩阵中
        betaMatrix(row,col) = beta;
        roMatrix(row,col) = ro;
        tcMatrix(row,col) = tc;
        VnoiseMatrix(row,col) = Vnoise;
        
    end
end
delete(gcp('nocreate')); 
end

