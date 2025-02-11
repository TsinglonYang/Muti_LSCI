function [roMatrix,tcMatrix] = SignalExpouse(skData,tkData,Tvalue,Beta)
% ���ع�ɢ�߳ıȳ���
%   �˴���ʾ��ϸ˵��
% ���� T ֵ�� Ave_sk �ṹ��
% Tvalue = 1 ms;

% ����ƥ��
crop_rowStart = floor((size(tkData,1) - size(skData,1))/2)+1;
crop_colStart = floor((size(tkData,2) - size(skData,2))/2)+1;

cropped_tk = tkData(crop_rowStart:(crop_rowStart+size(skData,1)-1),crop_colStart:(crop_colStart+size(skData,2)-1));

% ���� ro
roMatrix = 1-sqrt(abs(skData.^2 - cropped_tk.^2) / Beta);

% �� tau ����Ԥ����ռ�
[row,col] = size(roMatrix);
tcMatrix = zeros(row,col);

% ʹ�� ������� ��� tc
parpool('local', 18);
parfor l = 1:row
    for w = 1:col
        ro = roMatrix(l,w);
        ks = skData(l,w);
%         ks=ks^2
        model = @(tc,T) ((Beta*(ro^2)*((exp(-2*(Tvalue/tc)) - 1 + 2*(Tvalue/tc)) / (2*(Tvalue/tc)^2)) + ...
            4*Beta*ro*(1-ro)*((exp(-(Tvalue/tc)) - 1 + (Tvalue/tc)) / ((Tvalue/tc)^2)) + ...
            Beta*(1-ro)^2));
        tc = lsqcurvefit(model, ones(1,1), double(Tvalue), double(ks));
        tcMatrix(l,w) = tc;
    end
end
delete(gcp('nocreate')); 
end
