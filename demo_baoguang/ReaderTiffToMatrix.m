function data_struct = ReaderTiffToMatrix(InputFiles, StartFrame, EndFrame, OutputDataFormat)
% 批量读取 TIFF 文件到结构体矩阵中
% 输入参数：
%     InputFiles = 处理文件夹中的所有文件（前置使用 ChooseFile.m）
%     StartFrame = int 类型的整数，批量读取的开始帧
%     EndFrame = int 类型的整数，批量读取的结束帧
%     OutputDataFormat =  'double' | 'uint8' (8 bits depth) | uint16 (16 bits depth)
% 输出参数：
%     data_struct = 存储数据的结构体

tiff_num = length(InputFiles.tiff); % 求结构体长度

data_struct = struct('frames', cell(1, tiff_num));

% 循环遍历每个 TIFF 文件
for fileIndex = 1:tiff_num
    fileFullName = fullfile(InputFiles.tiff(fileIndex).folder, InputFiles.tiff(fileIndex).name);
    
    % 读取 TIFF 文件
    tiffInfo = imfinfo(fileFullName);
    tureFrameNum = length(tiffInfo);
    
    % 设置有效的开始帧和结束帧（怕输入参数不正确）
    [StartFrame, EndFrame] = setStartEndFrames(StartFrame, EndFrame, tureFrameNum);
    
    % 初始化当前文件的图像帧
    framesData = cell(1, EndFrame - StartFrame + 1);
    
    % 循环遍历帧图像
    for frameIndex = 1:(EndFrame - StartFrame + 1)
        % 读取当前帧
        frame = imread(fileFullName, frameIndex);
        
        % 存储到当前文件的图像帧数据中
        framesData{frameIndex} = setMatrixDataFormat(frame, OutputDataFormat);
    end
    
    % 将当前文件的图像帧数据存储到结构体中
    data_struct(fileIndex).frames = framesData;
    
end

end


%% 数据类型转换，转换为给定的输出格式
function InputMatrix = setMatrixDataFormat(InputMatrix, OutputDataFormat)
% Convert to the given output foramt（转换为给定的输出格式）

switch(OutputDataFormat)
    case 'uint8'
        cInt8 = 255; % coefficient to convert to 8 bit integer（将系数转换为8位整数）
        InputMatrix = cInt8.*(InputMatrix./max(InputMatrix, [], 'all'));
        InputMatrix = uint8(InputMatrix);
    case 'uint16'
        InputMatrix = uint16(InputMatrix);
    case 'double'
        InputMatrix = double(InputMatrix);
    otherwise
        fprintf('\n\nUnsupported output data type --> Data Type = %s\n', OutputDataFormat);
        error('Exit due to the error above!');
end

end


%% 设定有效的开始帧和结束帧
function [StartFrame, EndFrame] = setStartEndFrames(StartFrame, EndFrame, Frames)
% Check validity and set start and end frames（检查有效性，设置开始帧和结束帧）

if EndFrame < 1 || EndFrame > Frames
    EndFrame = Frames;
end

if StartFrame < 1 || StartFrame > Frames
    StartFrame = 1;
end

if StartFrame > EndFrame
    fprintf('\nSatrt frame index is bigger than End frame index --> StartFrame = %d, EndFrame = %d\n', StartFrame, EndFrame);
    error('Exit due to above error');
end

end