# RPT 피처 추출 상세 명세서 (Mathematical Formulation)

> **관련 스크립트 (Related Scripts)**:
> 1. **피처 추출**: `RPT_Feature_Extraction_Advanced.m` (메인 로직)
> 2. **상관 분석**: `RPT_Correlation_Heatmap.m`
> 3. **시각화**: `RPT_Pipeline_Visualization.m`

## 1. 데이터 전처리 (Data Preprocessing)

### 1.1 리샘플링 (균일 전압 그리드)
- **목적**: 모든 채널과 사이클의 데이터 비교 및 모델 입력을 위한 통일된 기준 마련
- **소스 데이터**: 원본 시계열 측정 데이터 (`.csv` 파일)
- **방법**: `0.001V` 간격으로 선형 보간 (Linear Interpolation)
- **필터링 여부**: **없음** (원본 데이터를 그대로 보간 후, 각 단계별로 별도 스무딩 적용)

---

## 2. 피처 추출 (Feature Extraction) - 총 14개

모든 피처는 **C-rate별 충/방전 곡선**에서 추출됩니다. (OCV 데이터 사용 안 함)

| No. | Feature Name | Data Source (데이터 소스) | Voltage Window (전압 윈도우) | Preprocessing (전처리) | Extraction Formula (추출 수식) |
|---|---|---|---|---|---|
| 1–5 | **Charge ΔQ** (Seg1~5) | **C-rate Charge Profile**<br>`Q_crate_chg(V)` | `3.70 V` ~ `3.95 V` | 0.001V Resampling | `F_dQ = |Q(V_end) - Q(V_start)|` |
| 6–10 | **Discharge ΔQ** (Seg1~5) | **C-rate Discharge Profile**<br>`Q_crate_dch(V)` | `3.75 V` ~ `3.88 V` | 0.001V Resampling | `F_dQ = |Q(V_end) - Q(V_start)|` |
| 11 | **Chg PkH** (Height) | **C-rate Charge Profile**<br>`Q_crate_chg(V)` | `3.70 V` ~ `3.95 V` | **Moving Average**<br>(Window 21) | `max( dQ/dV_filt )` |
| 12 | **Chg PkA** (Area) | **C-rate Charge Profile**<br>`Q_crate_chg(V)` | `3.70 V` ~ `3.95 V` | **Moving Average**<br>(Window 21) | `trapz( dQ/dV_filt )` |
| 13 | **Dch PkH** (Height) | **C-rate Discharge Profile**<br>`Q_crate_dch(V)` | `3.75 V` ~ `3.88 V` | **Moving Average**<br>(Window 21) | `max( dQ/dV_filt )` |
| 14 | **Dch PkA** (Area) | **C-rate Discharge Profile**<br>`Q_crate_dch(V)` | `3.75 V` ~ `3.88 V` | **Moving Average**<br>(Window 21) | `trapz( dQ/dV_filt )` |

---

### 2.3 데이터 로드 및 처리 시점

| 단계 | 상태 | 설명 |
|---|---|---|
| **데이터 로드** | **0.001V Resampled** | `RPT_VQ_grid.mat`에서 로드되는 데이터는 **SG 필터가 적용되지 않은 상태**입니다. |
| **피처 추출 시** | **Dynamic Filtering** | 코드 실행 중(`extract_features_half`)에 `gradient` 계산 후 **실시간으로 Moving Average**를 적용합니다. |

---

## 3. 라벨 추출 (Label Generation) - 총 3개

라벨은 배터리의 물리적/화학적 상태를 나타내며, **OCV 데이터** 또는 **정적 용량 데이터**에서 추출됩니다.

| Label Name | Data Source (데이터 소스) | ROI Window (분석 구간) | Preprocessing (전처리) | Extraction Method (추출 방법) |
|---|---|---|---|---|
| **SOH** | **Static Capacity Test**<br>(`Capacity_Data_Static.mat`) | `3.0 V` ~ `4.2 V`<br>(전체 구간) | N/A | `max(Q_discharge)`<br>(최대 방전 용량) |
| **LLI** | **OCV Charge Profile**<br>`Q_ocv_chg(V)` | `3.40 V` ~ `4.00 V` | **Moving Average**<br>(Window 21, 0.021V) | `V_peak(Classic) - V_peak(Fresh)`<br>(피크 전압 이동) |
| **LAM** | **OCV Charge Profile**<br>`Q_ocv_chg(V)` | `3.40 V` ~ `4.00 V` | **Moving Average**<br>(Window 21, 0.021V) | `H_peak(Classic) - H_peak(Fresh)`<br>(피크 높이 변화) |

## 4. Outlier 분석 및 Peak Shift 검증 (Phase 6)

상관계수 저하의 원인인 비단조성(Non-monotonic) 샘플을 시각적으로 검증합니다.

| 분석 대상 | 특이 사항 | 원인 가설 |
|---|---|---|
| **Ch16** | Seg3(0.1C)에서 400cyc 이후 용량 급증 | **Peak Shift**: 물리적 피크가 고정 전압 세그먼트를 통과함 |
| **Ch10** | Seg1(2C, 3C) 초반 용량 치솟음 | 초기 활성화 및 전압 윈도우 경계 거동 |

### 검증 방법
- **Ch14(Normal)**와 **Ch16(Outlier)**의 $V$-$Q$ 곡선을 겹쳐서 비교.
- 동일 전압 구간 내에서 사이클 진행에 따른 곡선 형태(기울기, 피크 위치) 시프트 확인.
- 시작 용량을 `0`으로 정규화하여 델타 Q($dQ$)의 전압별 변화 시각화.

> [!IMPORTANT]
> **결론**: 이러한 물리적 **피크 시프트(Peak Shift)**가 고정된 전압 세그먼트 경계를 넘나들기 때문에 상관계수가 낮게 측정될 수 있습니다. 이는 데이터 오류가 아닌 NMC 배터리의 물리적 특성입니다.

## 5. 모델 학습 및 검증 전략 (Phase 7): PCA + MLR & RF (LOCO-CV)

> **스크립트**: `RPT_Modeling.m` (신규 생성)

### 5.0 데이터 구성

| 항목 | 내용 |
|---|---|
| **입력 ($X$)** | `X_Normalized` (14개 피처: 충전dQ×5, 방전dQ×5, Chg PkH/A×2, Dch PkH/A×2) |
| **라벨 ($Y$)** | `Y_Labels` (SOH, LLI, LAM) — LLI/LAM은 **OCV Charge** 기준 |
| **그룹** | `CellID` (Ch09~Ch16, 총 8개 그룹) |
| **C-rate** | 5개 C-rate 데이터를 **모두 합쳐서** 하나의 모델로 학습 |

---

### 5.1 [1단계] 전체 데이터 정합성 확인 (Sanity Check)

전체 8개 채널 데이터를 모두 사용하여 피처가 라벨을 설명할 수 있는지 빠르게 확인한다.

1. **PCA**: 14개 피처 → 누적 설명력 **95%** 기준으로 PC 개수(k) 결정
2. **MLR**: k개 PC로 전체 데이터 학습 → R² 확인 (Train = Test, 과적합 허용)
3. **출력**: Scree Plot, 라벨별 R², Predicted vs Actual 산점도

---

### 5.2 [2단계] LOCO-CV (Leave-One-Cell-Out Cross-Validation)

```
for i = 1:8  (각 셀을 한 번씩 Test로 지정)
  ├── Train = 나머지 7개 셀
  ├── Test  = i번째 셀
  │
  ├── [스케일링] Train의 μ, σ 로 정규화 → Test에 동일 적용
  ├── [PCA] Train으로 PC 축 학습 → Test 투영 (Data Leakage 방지)
  │
  ├── [MLR] Train 학습 → Test 예측 → RMSE/MAE/R² 저장
  │
  └── [RF]  Train 내 HP 그리드 서치 (NumLearningCycles × MinLeafSize)
            → 최적 HP로 Train 재학습 → Test 예측 → RMSE/MAE/R² 저장
end
```

#### HP 튜닝 그리드 (RF)

| 파라미터 | 탐색 범위 |
|---|---|
| `NumLearningCycles` (나무 개수) | [50, 100, 200, 500] |
| `MinLeafSize` (최소 잎사귀 크기) | [1, 3, 5, 10] |

---

### 5.3 최종 출력물

| 출력물 | 설명 |
|---|---|
| **성적표** | 8회 루프 평균 RMSE, MAE, R² (SOH, LLI, LAM 각각) |
| **추정 그래프** | 가장 오차가 **컸던 셀** vs **잘 맞춘 셀** 비교 |
| **Scree Plot** | PCA 누적 설명력 그래프 |
| **Feature Importance** | RF의 피처 중요도 분석 |

---

### 5.4 스크립트 섹션 구성

```
RPT_Modeling.m
├── Section 1: 데이터 로드 & 전처리
├── Section 2: PCA (Scree Plot, 95% 기준 PC 선택)
├── Section 3: [1단계] 전체 데이터 MLR Sanity Check
├── Section 4: [2단계] LOCO-CV 루프 (MLR + RF)
├── Section 5: 성능 요약 테이블 출력
└── Section 6: 추정 그래프 & Feature Importance 시각화
```
