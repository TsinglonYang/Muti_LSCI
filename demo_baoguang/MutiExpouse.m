function [betaMatrix,roMatrix,tcMatrix,VnoiseMatrix] = MutiExpouse(ave_skData,T_values)
% ���ع�ɢ�߳ıȳ���
%   �˴���ʾ��ϸ˵��
% ���� T ֵ�� Ave_sk �ṹ��
% T_values = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 12; 14; 16; 18; 20];

[L,W] = size(ave_skData(1).average_sk);

% �����������洢���
betaMatrix = zeros(L,W);
roMatrix = zeros(L,W);
tcMatrix = zeros(L,W);
VnoiseMatrix = zeros(L,W);

parpool('local', 18);
parfor row = 1:L
    for col = 1:W
        Ks_value = zeros(numel(T_values),1);
        
        % ��ȡÿ���ع�ʱ���� ave_sk �� (row,col) ����ֵ
        for i = 1:numel(T_values)
            ave_skMatrix = ave_skData(i).average_sk; % ��ȡ��i������
            Ks_value(i) = ave_skMatrix(row,col); % ��ȡÿ������� (row,col) λ�õ�ֵ
        end
        
        model = @(params, T) sqrt((params(1).*(params(2).^2).*((exp(-2*(T./params(3)))-1+2*(T./params(3)))./(2*(T./params(3)).^2)) +...
            4.*params(1).*params(2).*(1-params(2)).*((exp(-(T./params(3)))-1+(T./params(3)))./((T./params(3)).^2)) + ...
            params(1).*(1-params(2)).^2) +...
            params(4))
        
        % ���� Options �ṹ�壬�����ݲ�ֵ�͹ر���ʾ����
        options = optimset('TolFun', 1e-6, 'Display', 'off');
        
        % ���������
        lower_limit = [1e-9, 0, 0, 0]; % ���ޣ����Ը�����Ҫ����
        upper_limit = [Inf, Inf, 1, Inf]; % ���ޣ����Ը�����Ҫ����
        Initial_guess = [1e-6,0.5,0.9,0];
        
        % ���
        estimated_params = lsqcurvefit(model, Initial_guess, double(T_values), double(Ks_value), lower_limit, upper_limit, options);

        % ��ȡ��Ϻ�Ĳ���
        beta = estimated_params(1);
        ro = estimated_params(2);
        tc = estimated_params(3);
        Vnoise = estimated_params(4);
        
        % �洢�������Ӧ�Ľ��������
        betaMatrix(row,col) = beta;
        roMatrix(row,col) = ro;
        tcMatrix(row,col) = tc;
        VnoiseMatrix(row,col) = Vnoise;
        
    end
end
delete(gcp('nocreate')); 
end

