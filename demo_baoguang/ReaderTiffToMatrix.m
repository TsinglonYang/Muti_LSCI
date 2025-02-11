function data_struct = ReaderTiffToMatrix(InputFiles, StartFrame, EndFrame, OutputDataFormat)
% ������ȡ TIFF �ļ����ṹ�������
% ���������
%     InputFiles = �����ļ����е������ļ���ǰ��ʹ�� ChooseFile.m��
%     StartFrame = int ���͵�������������ȡ�Ŀ�ʼ֡
%     EndFrame = int ���͵�������������ȡ�Ľ���֡
%     OutputDataFormat =  'double' | 'uint8' (8 bits depth) | uint16 (16 bits depth)
% ���������
%     data_struct = �洢���ݵĽṹ��

tiff_num = length(InputFiles.tiff); % ��ṹ�峤��

data_struct = struct('frames', cell(1, tiff_num));

% ѭ������ÿ�� TIFF �ļ�
for fileIndex = 1:tiff_num
    fileFullName = fullfile(InputFiles.tiff(fileIndex).folder, InputFiles.tiff(fileIndex).name);
    
    % ��ȡ TIFF �ļ�
    tiffInfo = imfinfo(fileFullName);
    tureFrameNum = length(tiffInfo);
    
    % ������Ч�Ŀ�ʼ֡�ͽ���֡���������������ȷ��
    [StartFrame, EndFrame] = setStartEndFrames(StartFrame, EndFrame, tureFrameNum);
    
    % ��ʼ����ǰ�ļ���ͼ��֡
    framesData = cell(1, EndFrame - StartFrame + 1);
    
    % ѭ������֡ͼ��
    for frameIndex = 1:(EndFrame - StartFrame + 1)
        % ��ȡ��ǰ֡
        frame = imread(fileFullName, frameIndex);
        
        % �洢����ǰ�ļ���ͼ��֡������
        framesData{frameIndex} = setMatrixDataFormat(frame, OutputDataFormat);
    end
    
    % ����ǰ�ļ���ͼ��֡���ݴ洢���ṹ����
    data_struct(fileIndex).frames = framesData;
    
end

end


%% ��������ת����ת��Ϊ�����������ʽ
function InputMatrix = setMatrixDataFormat(InputMatrix, OutputDataFormat)
% Convert to the given output foramt��ת��Ϊ�����������ʽ��

switch(OutputDataFormat)
    case 'uint8'
        cInt8 = 255; % coefficient to convert to 8 bit integer����ϵ��ת��Ϊ8λ������
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


%% �趨��Ч�Ŀ�ʼ֡�ͽ���֡
function [StartFrame, EndFrame] = setStartEndFrames(StartFrame, EndFrame, Frames)
% Check validity and set start and end frames�������Ч�ԣ����ÿ�ʼ֡�ͽ���֡��

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