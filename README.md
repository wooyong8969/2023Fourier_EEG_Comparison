# 뇌파 분석에서 세 가지 푸리에 변환 방식의 차이 분석

> 기간: 2023년  
> 형태: 고등학교 심화수학I 수행평가  
> 주제: 푸리에 변환의 이론과 EEG 데이터 적용 비교 분석

## 개요

본 프로젝트는 『수학으로 배우는 파동의 법칙』을 읽은 후,  
파동과 관련된 수학적 도구인 푸리에 변환에 대한 심화 이해를 목표로 수행되었습니다.  
특히 이산 푸리에 변환(DFT), 고속 푸리에 변환(FFT), 단시간 푸리에 변환(STFT)  
각각의 원리와 장단점을 비교하고, 실제 뇌파 데이터에 적용하여 효율성과 특성을 분석하고자 하였습니다.

※ 본 탐구는 학교 수학 수행평가의 일환으로 진행되었으며, 수학적 이론과 데이터 분석의 융합을 목표로 하였습니다.

## 사용 데이터

- KaraOne Dataset (University of Toronto)  
  - 참가자가 특정 음운(iy, uw, piy, tiy, diy 등)을 발음하거나 상상할 때의 EEG 및 음성 데이터가 포함되어 있습니다.
  - [KaraOne 데이터셋 링크](https://www.cs.toronto.edu/~complingweb/data/karaOne/karaOne.html)

- BCI Competition IV Dataset 2a  
  - http://dx.doi.org/10.5524/100542

## 이론적 배경

### 1. 푸리에 급수

- 복잡한 주기적 파동은 sin, cos 함수의 합으로 표현할 수 있으며, 이를 푸리에 급수라고 합니다.
- 푸리에 계수는 삼각함수의 직교성을 이용하여 계산할 수 있습니다.

### 2. 푸리에 변환

- **DFT (Discrete Fourier Transform)**: 이산 신호를 주파수 영역으로 변환
- **FFT (Fast Fourier Transform)**: DFT의 계산을 효율적으로 수행하는 알고리즘
- **STFT (Short-Time Fourier Transform)**: 시간에 따른 주파수 변화를 분석할 수 있는 기법

### 3. 오일러 공식

- 복소수를 이용해 sin과 cos 파동을 하나의 지수 함수로 표현할 수 있습니다.  
- 푸리에 급수와 변환의 복소 표현 도출에 활용되었습니다.

## 탐구 방법

1. KaraOne 음성 및 뇌파 데이터 수집
2. EEG 데이터에 대해 DFT, FFT, STFT 각각 적용
3. FFT는 MATLAB 내장 함수를 이용, DFT와 STFT는 직접 수식에 따라 구현
4. 각각의 변환 결과를 시각화하여 비교
5. 계산 시간과 변환 결과 품질을 평가

## 결과

- DFT, FFT, STFT 모두 뇌파 데이터의 주요 주파수 성분을 유사하게 포착하는 것으로 나타났습니다.
- 계산 시간 면에서는 FFT가 가장 우수한 성능을 보였습니다.
- STFT는 시간-주파수 해상도 분석에 유리하였으나, 신호를 나누는 윈도우 크기 설정에 따라 성능이 크게 변동할 수 있음을 확인하였습니다.

## 한계점 및 향후 연구

- 세 가지 변환 방식 모두 매우 유사한 결과를 보여, 장단점 분석이 다소 어려웠습니다.
- 향후 연구에서는 STFT 윈도우 크기 및 중첩 조정, 웨이브릿 변환 도입 등을 통해 시간-주파수 해상도 비교를 강화할 계획입니다.
- 노이즈 제거, 신호 분리 능력 평가 등 추가적인 성능 지표를 도입할 필요가 있습니다.

## 탐구 소감

본 탐구를 통해 이론적 수학 지식(삼각함수, 푸리에 급수, 오일러 공식)과  
신호 처리 실습(푸리에 변환)을 연계하는 경험을 쌓을 수 있었습니다.  
특히, 이론에서 실제 EEG 분석까지 연결하는 과정을 통해  
수학의 실제 응용 가능성과 한계를 직접 체감할 수 있었다는 점에서 의미가 컸습니다.

