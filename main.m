
clear all; close all; clc;
addpath(genpath('D:\WOOYONG\MATLAB\PhDYeom\Ground'));

% 파라미터 설정
sf = 1000;                % 샘플링 주파수
ch_n = 62;                % 채널 수
wnd_size = [-1 4];        % 윈도우 크기 (초 단위)
baseline = [-1 0];        % 기준선
f_scale = 1;              % 주파수 스케일
freq_band = [0.1 100];    % 분석할 주파수 대역
normal = 1;               % 기준선 정규화 여부
fullscreen = get(0, 'ScreenSize'); % 전체 화면 크기
position = EEG_62ch_layout_Brain_Products; % EEG 채널 위치

% 데이터 로드
load('sess01_subj36_EEG_MI.mat'); % 학습용 데이터
EEG = EEG_MI_train.x;             % EEG 데이터 단축 이름

% Common Average Referencing (CAR)
EEG = EEG - repmat(mean(EEG, 1), ch_n, 1);

% 이벤트 추출
events{1} = EEG_MI_train.t(find(EEG_MI_train.y_dec == 1)); % 오른손
events{2} = EEG_MI_train.t(find(EEG_MI_train.y_dec == 2)); % 왼손

% EEG 에포킹
for i = 1:length(events)
    for tr = 1:size(events{i}, 2)
        e_EEG{i}(:,:,tr) = EEG(:, round(events{i}(tr) + wnd_size(1) * sf): ...
                                   round(events{i}(tr) + wnd_size(2) * sf));
    end
end

% 실행 시간 저장용 초기화
FFT_times = [];
DFT_times = [];
STFT_times = [];

% FFT, DFT, STFT 구현
for i = 1:length(events)
    for tr = 1:size(e_EEG{i}, 3)
        % FFT
        tic;
        FFT_results{i}(:,:,tr) = fft(e_EEG{i}(:,:,tr), [], 2);
        FFT_times = [FFT_times, toc];

        % DFT
        tic;
        N = size(e_EEG{i}, 2);
        dft_matrix = exp(-2i * pi * (0:N-1).' * (0:N-1) / N);
        DFT_results{i}(:,:,tr) = e_EEG{i}(:,:,tr) * dft_matrix';
        DFT_times = [DFT_times, toc];

        % STFT
        tic;
        window_length = 256;
        for ch = 1:ch_n
            [S, F, T] = spectrogram(e_EEG{i}(ch,:,tr), window_length, window_length/2, window_length, sf);
            STFT_results{i}(ch,:,:,tr) = S;
        end
        STFT_times = [STFT_times, toc];
    end
end

% 시각화할 채널 선택
channel_to_display = 1;

% 결과 시각화
for i = 1:length(events)
    figure('Position', fullscreen);
    for tr = 1:size(FFT_results{i}, 3)
        subplot(5, 10, tr);
        plot(abs(FFT_results{i}(channel_to_display,:,tr)));
        title([num2str(tr)]);
    end

    figure('Position', fullscreen);
    for tr = 1:size(DFT_results{i}, 3)
        subplot(5, 10, tr);
        plot(abs(DFT_results{i}(channel_to_display,:,tr)));
        title([num2str(tr)]);
    end

    figure('Position', fullscreen);
    for tr = 1:size(STFT_results{i}, 4)
        subplot(5, 10, tr);
        imagesc(T, F, abs(STFT_results{i}(channel_to_display,:,:,tr)));
        axis tight;
        colormap jet;
        title([num2str(tr)]);
    end
end
