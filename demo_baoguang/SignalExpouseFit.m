function [betaResult,roResult,tcResult,VnoiseResult] = SignalExpouseFit(ave_skData,T_values)
% ���ع�ɢ�߳ıȳ��������⣨ÿһ�� T ���Ӧ�� ks ��ϳ�һ�� tc ͼ��15 �� T �� ks ��ϳ� 15 �� tc ͼ��
% ���������
%     ave_skData = �洢ɢ�߳ı����ݵĽṹ�壬ǰ��ʹ�� Sum_calcSpeckleContrast.m �ļ�
%     T_values = [1; 2; 3] | 1 ���������������ͬʱ���㣬Ҳ������һ��������һ����
% ���������
%     [betaResult,roResult,tcResult,VnoiseResult] = �洢��ϳ�����Ľṹ��


% �����ṹ�����飬���ڴ洢ÿ�����ݵ���ϲ���
betaResult = struct('result', cell(1, numel(T_values)));
roResult = struct('result', cell(1, numel(T_values)));
tcResult = struct('result', cell(1, numel(T_values)));
VnoiseResult = struct('result', cell(1, numel(T_values)));

[row,col]=size(ave_skData(1).average_sk);

% �����������洢���
betaMatrix = zeros(row ,col);
roMatrix = zeros(row,col);
tcMatrix = zeros(row,col);
VnoiseMatrix = zeros(row,col);

% ����������
h = waitbar(0, 'Fitting Progress...');
for tIndex = 1:numel(T_values)
    T = T_values(tIndex);
    
    % ��ȡ��Ӧ T ֵ�µ� ave_skData ����
    ave_skMatrix = ave_skData(tIndex).average_sk; % ����ṹ���а�����Ϊ average_sk �ľ����ֶ�

    % ���ڲ�ѭ���б��� ave_skData ����
    parfor l = 1:row
        for w = 1:col
            ks_value = ave_skMatrix(l,w) ^ 2;
            model = @(params, T) sqrt((params(1).*(params(2).^2).*((exp(-2*(T./params(3)))-1+2*(T./params(3)))./(2*(T./params(3)).^2)) +...
                4.*params(1).*params(2).*(1-params(2)).*((exp(-(T./params(3)))-1+(T./params(3)))./((T./params(3)).^2)) + ...
                params(1).*(1-params(2)).^2) +...
                params(4));
            
            % ���� Options �ṹ�壬�����ݲ�ֵ�͹ر���ʾ����
            options = optimset('TolFun', 1e-6, 'Display', 'off');
            
            % ���������
            lower_limit = [1e-9, 0, 0, 0]; % ���ޣ����Ը�����Ҫ����
            upper_limit = [Inf, Inf, 1, Inf]; % ���ޣ����Ը�����Ҫ����
            Initial_guess = [1e-6,0.5,0.9,0];
            
            % ���
            estimated_params = lsqcurvefit(model, Initial_guess, double(T), double(ks_value), lower_limit, upper_limit, options);

            % ��ȡ��Ϻ�Ĳ���
            beta = estimated_params(1);
            ro = estimated_params(2);
            tc = estimated_params(3);
            Vnoise = estimated_params(4);
            
            % �洢�������Ӧ�Ľ��������
            betaMatrix(l,w) = beta;
            roMatrix(l,w) = ro;
            tcMatrix(l,w) = tc;
            VnoiseMatrix(l,w) = Vnoise;
        end
    end
    
    % �洢�������Ӧ�Ľṹ��
    betaResult(tIndex).result = betaMatrix;
    roResult(tIndex).result = roMatrix;
    tcResult(tIndex).result = tcMatrix;
    VnoiseResult(tIndex).result = VnoiseMatrix;
    
    % ���½�����
    waitbar(tIndex / numel(T_values), h);
end

% �رս�����
close(h);
end
