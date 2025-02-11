function [tkData, skData, ave_skData] = Sum_calcSpeckleContrast(data_struct, Method, Average, Parameter)
% ���������� .tiff �ļ��еĶ���ԭʼɢ��ͼ��ļ��㷽������
% ���������
%     data_struct = �洢���ݵĽṹ�壬ǰ��ʹ�� ReaderTiffToMatrix.m �ļ�
%     Method = "sk" | "awk" | "sdk" | "awsdk" | "tk"
%     Average = "average_sk"
%     Parameter = �������������Ϊһ��ʾ��

%{
para.sk.window = 3;
para.sk.method = 'ConvFilter';
para.sdk.window = 9;
para.awk.filter_window = [11,11];
para.awk.window = 5;
para.awk.cluster = 3;
para.awsdk.window = 9;
para.awsdk.cluster = 3;
para.tk.window = 30;
para.tk.method = 'Discrete';
%}

% ���������
%     [tkData, skData] = �洢 sk��tk ���ݵĽṹ��

skData = struct();
tkData = struct();
len = length(data_struct); % tiff �ļ���
ave_skData = struct('average_sk', cell(len, 1)); % ����һ���µĽṹ�����飬ÿ���ṹ�����һ�� average_sk �ֶ�

% sk��tk ���о���ת����������
if Method == "sk" || Method == "tk"
    % �������������Ѿ������� imageData �ṹ���У�������ѭ������ÿ��cell
    for cellIndex = 1:len
        % �洢�ڽṹ���е�֡��
        RawSpeckleData2D = data_struct(cellIndex).frames;
        numFrames = length(RawSpeckleData2D);
        % ԭʼɢ��ͼ��� row��col
        [matrixRowSize,  matrixColSize] = size(RawSpeckleData2D{1});

        % ����һ���յ�3D���飬����ά��Ӧ��ͬ��֡��
        RawSpeckleData = zeros([matrixRowSize, matrixColSize, numFrames]);

        % �� cell �����еľ����Ƶ� 3D ������
        for i = 1:numFrames
            RawSpeckleData(:, :, i) = cell2mat(RawSpeckleData2D(i));
        end

        % ���ڣ�RawSpeckleData3D ��һ�� [matrixRowSizexmatri,ColSize,numFrames] �� 3D ����
        
        switch Method
            case "sk"
                sk = SpatialContrastMethod(RawSpeckleData, Parameter.sk.window, Parameter.sk.method);
                for frameIndex = 1:numFrames
                    skData(cellIndex).skValues{frameIndex} = sk(:, :, frameIndex);
                end
            case "tk"
                tk = TemporalContrastMethod(RawSpeckleData, Parameter.tk.window, Parameter.tk.method);
                tkData(cellIndex).tkValues = tk;
        end
    end
    
% �����ֿռ�ıȳ��񷽷�ֱ�Ӽ���
elseif Method == "sdk" || Method == "awk" || Method == "awsdk"
    switch Method
        case "sdk"
            sdk_window = Parameter.sdk.window;
            % ÿ�� TIFF �ļ���������
            parpool('local', 6);
            parfor cellIndex = 1:len
                % ѭ������ÿ�� TIFF �ļ��� ��֡ �ĳı�ͼ��
                for subFrameNum = 1:length(data_struct(cellIndex).frames)
                    skData(cellIndex).skValues{subFrameNum} = SpaceDirectionalContrastMethod_paper_mex(data_struct(cellIndex).frames{subFrameNum},sdk_window);
                end
            end
            delete(gcp('nocreate')); 
            
        case "awk"
            sk_window = Parameter.sk.window;
            sk_method = Parameter.sk.method;
            
            awk_filter_window = Parameter.awk.filter_window;
            
            awk_window = Parameter.awk.window;
            awk_cluster = Parameter.awk.cluster;
            
            % ÿ�� TIFF �ļ���������
            parpool('local', 6);
            parfor cellIndex = 1:len
                % 1������ ave_sk
                [row, col] = size(data_struct(cellIndex).frames{1});
                subFrameNums = length(data_struct(cellIndex).frames);
                ave_sk = zeros(row - sk_window + 1, col - sk_window + 1);
                for subFrameNum = 1:subFrameNums
                    ave_sk = ave_sk + SpatialContrastMethod(data_struct(cellIndex).frames{subFrameNum}, sk_window, sk_method);
                end
                ave_sk = ave_sk ./ subFrameNums;
                
                % 2����ֵ�˲�
                filter_ave_sk = medfilt2(ave_sk, awk_filter_window);
                
                % 3���� awk ���� ��֡ ֱ�� ��ֵ������
                for subFrameNum = 1:subFrameNums
                    skData(cellIndex).skValues{subFrameNum} = AdaptiveWindowContrastMethod_mex(data_struct(cellIndex).frames{subFrameNum}, filter_ave_sk, awk_window, awk_cluster);
                end
                
            end
            delete(gcp('nocreate')); 
            
        case "awsdk"
            sk_window = Parameter.sk.window;
            sk_method = Parameter.sk.method;
            
            awsdk_window = Parameter.awsdk.window;
            awsdk_cluster = Parameter.awsdk.cluster;
            
            % ÿ�� TIFF �ļ���������
            parpool('local', 6);
            parfor cellIndex = 1:len
                % ѭ������ÿ�� TIFF �ļ��� ��֡ �ĳı�ͼ��
                for subFrameNum = 1:length(data_struct(cellIndex).frames)
                    % 1���� sk ���� ��֡
                    temp_sk = SpatialContrastMethod(data_struct(cellIndex).frames{subFrameNum}, sk_window, sk_method);
                    % 2������������ɢ�˲� ���� ��֡
                    skFilter = AnisotropicDiffusionFilter_mex(temp_sk,50,0.25,"self adaption");
                    % 3���� awsdk ���� ��֡ ֱ�� ��ֵ������
                    skData(cellIndex).skValues{subFrameNum} = idea2_mex(data_struct(cellIndex).frames{subFrameNum},skFilter, awsdk_window, awsdk_cluster);
                end
            end
            delete(gcp('nocreate')); 
    end
    
end

if Average == "average_sk"
    for cellIndex = 1:len
        % ���� skData(cellIndex).skValues ���� 60 �� 1496x796 �ľ���
        skValues = skData(cellIndex).skValues;
    
        % ��ʼ��һ��������С��ͬ�ľ������ڴ洢ƽ��ֵ
        averageMatrix = zeros(size(skValues{1}));

        % ���� 60 �������ƽ��ֵ
        for frameIndex = 1:length(skData(cellIndex).skValues)
            averageMatrix = averageMatrix + skValues{frameIndex};
        end
        averageMatrix = averageMatrix / length(skData(cellIndex).skValues); % ��ƽ��ֵ

        % ��ƽ������洢��Ave_sk�ṹ����
        ave_skData(cellIndex).average_sk = averageMatrix;
    end
    fprintf('\nsk average !');
else
    fprintf('\nno average !');
    return

end

end