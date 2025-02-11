function [tkData, skData, ave_skData] = Sum_calcSpeckleContrast(data_struct, Method, Average, Parameter)
% 批量计算多个 .tiff 文件中的多张原始散斑图像的计算方法汇总
% 输入参数：
%     data_struct = 存储数据的结构体，前置使用 ReaderTiffToMatrix.m 文件
%     Method = "sk" | "awk" | "sdk" | "awsdk" | "tk"
%     Average = "average_sk"
%     Parameter = 输入参数，以下为一个示例

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

% 输出参数：
%     [tkData, skData] = 存储 sk、tk 数据的结构体

skData = struct();
tkData = struct();
len = length(data_struct); % tiff 文件数
ave_skData = struct('average_sk', cell(len, 1)); % 创建一个新的结构体数组，每个结构体包含一个 average_sk 字段

% sk、tk 进行矩阵转换后计算更快
if Method == "sk" || Method == "tk"
    % 假设您的数据已经存在于 imageData 结构体中，您可以循环遍历每个cell
    for cellIndex = 1:len
        % 存储在结构体中的帧数
        RawSpeckleData2D = data_struct(cellIndex).frames;
        numFrames = length(RawSpeckleData2D);
        % 原始散斑图像的 row，col
        [matrixRowSize,  matrixColSize] = size(RawSpeckleData2D{1});

        % 创建一个空的3D数组，第三维对应不同的帧数
        RawSpeckleData = zeros([matrixRowSize, matrixColSize, numFrames]);

        % 将 cell 数组中的矩阵复制到 3D 数组中
        for i = 1:numFrames
            RawSpeckleData(:, :, i) = cell2mat(RawSpeckleData2D(i));
        end

        % 现在，RawSpeckleData3D 是一个 [matrixRowSizexmatri,ColSize,numFrames] 的 3D 数组
        
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
    
% 这三种空间衬比成像方法直接计算
elseif Method == "sdk" || Method == "awk" || Method == "awsdk"
    switch Method
        case "sdk"
            sdk_window = Parameter.sdk.window;
            % 每个 TIFF 文件独立计算
            parpool('local', 6);
            parfor cellIndex = 1:len
                % 循环计算每个 TIFF 文件下 子帧 的衬比图像
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
            
            % 每个 TIFF 文件独立计算
            parpool('local', 6);
            parfor cellIndex = 1:len
                % 1、先算 ave_sk
                [row, col] = size(data_struct(cellIndex).frames{1});
                subFrameNums = length(data_struct(cellIndex).frames);
                ave_sk = zeros(row - sk_window + 1, col - sk_window + 1);
                for subFrameNum = 1:subFrameNums
                    ave_sk = ave_sk + SpatialContrastMethod(data_struct(cellIndex).frames{subFrameNum}, sk_window, sk_method);
                end
                ave_sk = ave_sk ./ subFrameNums;
                
                % 2、中值滤波
                filter_ave_sk = medfilt2(ave_sk, awk_filter_window);
                
                % 3、算 awk —— 单帧 直接 赋值给矩阵
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
            
            % 每个 TIFF 文件独立计算
            parpool('local', 6);
            parfor cellIndex = 1:len
                % 循环计算每个 TIFF 文件下 子帧 的衬比图像
                for subFrameNum = 1:length(data_struct(cellIndex).frames)
                    % 1、算 sk —— 单帧
                    temp_sk = SpatialContrastMethod(data_struct(cellIndex).frames{subFrameNum}, sk_window, sk_method);
                    % 2、各向异性扩散滤波 —— 单帧
                    skFilter = AnisotropicDiffusionFilter_mex(temp_sk,50,0.25,"self adaption");
                    % 3、算 awsdk —— 单帧 直接 赋值给矩阵
                    skData(cellIndex).skValues{subFrameNum} = idea2_mex(data_struct(cellIndex).frames{subFrameNum},skFilter, awsdk_window, awsdk_cluster);
                end
            end
            delete(gcp('nocreate')); 
    end
    
end

if Average == "average_sk"
    for cellIndex = 1:len
        % 假设 skData(cellIndex).skValues 包含 60 个 1496x796 的矩阵
        skValues = skData(cellIndex).skValues;
    
        % 初始化一个与矩阵大小相同的矩阵，用于存储平均值
        averageMatrix = zeros(size(skValues{1}));

        % 计算 60 个矩阵的平均值
        for frameIndex = 1:length(skData(cellIndex).skValues)
            averageMatrix = averageMatrix + skValues{frameIndex};
        end
        averageMatrix = averageMatrix / length(skData(cellIndex).skValues); % 求平均值

        % 将平均矩阵存储在Ave_sk结构体中
        ave_skData(cellIndex).average_sk = averageMatrix;
    end
    fprintf('\nsk average !');
else
    fprintf('\nno average !');
    return

end

end