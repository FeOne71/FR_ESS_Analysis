---
description: Lab_RPT_Analysis Pipeline 및 Workflow 상세 플랜
---

# Lab RPT Analysis 파이프라인 워크플로우 (RPT_model_workflow)

## 1. 개요
* **목적**: Lab 환경에서 확보한 RPT(Reference Performance Test) 데이터 기반, 배터리 노화 지표(SOH, LLI, LAM) 추정을 위한 머신러닝 모델(Random Forest, MLR) 구축
* **Input Data**:
  * 0.001V 단위로 일정하게 선형 보간(Linear Interpolation)된 V-Q(Voltage-Capacity) 데이터 (`RPT_VQ_grid.mat`)
  * 각 채널/사이클별 Static Capacity 평가 결과 (`Capacity_Data_Static.mat`)
* **핵심 방법론**: 
  * 기준 구간(Voltage Window) 내에서 시간을 균등 분할(Time-Balanced)하여 개별 채널의 편차를 줄이는 **Global Master Ruler** 방식 적용.
  * 이동평균(Moving Average, Window: RPT 20, DC 20)을 통한 노이즈 제거 후 `dQ` 및 `dQ/dV` Feature 추출.
  * PCA (Principal Component Analysis) 기반 차원 축소와 Random Forest 및 MLR(Multiple Linear Regression)을 활용한 5-Fold Cross-Validation 모델링.

---

## 2. Pre-processing 및 Ruler 생성
**Target Script:** `RPT_Feature_Extraction_Advanced.m`

### 2.1. 이동평균 적용 및 0.001V Uniform Grid 보간
원본 Raw 데이터의 비균일한 Sample 간격을 해결하고 노이즈를 억제하기 위해, 전처리 단계에서 먼저 **이동평균(Moving Average, Window: RPT 20, DC 20)** 을 적용하여 데이터를 스무딩(Smoothing)한 후, 0.001V 단위의 Uniform Voltage Grid로 선형 보간($V_{grid}$)된 데이터를 사용한다.

### 2.2. Voltage Window 설정
분석 대상이 되는 C-rate 데이터(예: 0.1C, 0.5C 등)의 충/방전 구간을 다음과 같이 제한한다.
* **Charge**: 3.70 V ~ 3.95 V
* **Discharge**: 3.75 V ~ 3.88 V
*(해당 Window를 벗어나는 데이터는 분석에서 제외하여 노이즈 개입 최소화)*

### 2.3. Global Master Ruler 생성 (Time-Balanced)
충/방전 곡선의 특징을 일관되게 추출하기 위해 단순 Voltage 등분할이 아닌 시간 균등 분할 방식을 채택한다.
1. `cyc0`의 Static 방전 데이터를 기준으로 8개 채널의 Time Grid를 추출한다.
2. 8개 채널의 V-t 곡선을 평균화하여 **대표 평균 V-t 곡선**을 생성한다.
3. 평균 시간 축을 기준으로 5개 Segments로 **균등 분할(Equally divided by Time)** 한다.
4. 균등 분할된 시간 경계에 대응하는 Voltage 경계값(`V_bounds`)을 산출하여 전 채널, 전 사이클에 공통으로 적용되는 **Global Master Ruler**를 만든다.

---

## 3. Feature 추출 및 Standardization
**Target Script:** `RPT_Feature_Extraction_Advanced.m`

### 3.1. Feature 추출 (총 14개)
각 Cycle, 각 Channel, 각 C-rate 곡선마다 Global Master Ruler 전압 경계를 적용하여 구간 데이터를 추출한다.
* **Segment Capacity ($\Delta Q$, 10 Features):**
  * Charge 5개 Segments에 대한 용량 차이 ($\Delta Q$)
  * Discharge 5개 Segments에 대한 용량 차이 ($\Delta Q$)
* **Differential Capacity (dQ/dV, 4 Features):**
  * `gradient` 함수를 통해 이산 미분(Discrete Differentiation)을 수행한 뒤, **이동평균(Moving Average)** (Window: RPT 20, DC 20)을 적용해 Smoothing.
  * Charge dQ/dV 곡선의 **Peak Height**, **Peak Area**
  * Discharge dQ/dV 곡선의 **Peak Height**, **Peak Area**

### 3.2. Standardization (Z-score Normalization)
서로 다른 scale을 가진 Feature들을 모델에 입력하기 위해 표준 정규화(Z-score Normalization)를 수행한다.
* C-rate label 그룹별로 개별적인 Mean($\mu$)과 Standard Deviation($\sigma$) 계산
* 수식: $X_{norm} = \frac{X - \mu}{\sigma}$

---

## 4. Label 생성
추출된 Feature 행렬($X_{norm}$)과 매핑될 Target Variable인 SOH, LLI, LAM Label을 구한다.
* **SOH (Capacity)**: 해당 사이클 Static 평가의 Maximum Discharge Capacity.
* **LLI (Loss of Lithium Inventory)** & **LAM (Loss of Active Material)**: 
  * Fresh(cyc0) 상태의 OCV 곡선과 현재 OCV 곡선 간의 ICA(Incremental Capacity Analysis) Peak-tracking을 활용하거나 물리적 노화 방정식을 적용해 퍼센티지(%) 단위 지표로 추정.

---

## 5. Diagnostic Modeling
**Target Script:** `RPT_Modeling.m`

### 5.1. PCA (Principal Component Analysis)
14 Dimension의 Feature 간 다중공선성(Multicollinearity)을 방지하고 Noise를 억제하기 위해 Feature 행렬($X_{norm}$)에 대해 PCA를 수행한다.
* 95%의 Cumulative Variance를 설명하는 최소 수의 Principal Components (PC) 선택 (보통 5~7개).
* Scree Plot을 그려 개별 PC의 분산 기여도 검증.

### 5.2. 5-Fold Cross-Validation 설계
Data Leakage를 방지하고 범용적 검증을 위해 전체 Cell/Cycle 데이터를 무작위로 섞어 5-Fold CV 분할을 수행한다.
* 매 Fold마다 Train Set 기준으로 PCA Coefficients를 fit하고 이를 Test Set에 투영시켜 차원 축소 수행.

### 5.3. MLR (Multiple Linear Regression)
* 선택된 PC들을 독립 변수(Predictors)로 사용하여 SOH, LLI, LAM을 예측하는 단순 선형 모델 (Baseline Model).
* 1단계 전체 데이터 대상 Sanity Check 및 최종 CV 비교군으로 활용.

### 5.4. Random Forest 
* Bagging 방식의 Ensemble Tree Algorithm (`fitrensemble` 활용).
* **Hyperparameter Tuning (Grid Search)**:
  * OOB (Out-of-Bag) RMSE를 최소화하는 파라미터 조합 탐색.
  * Search Space: `NumTrees` = [50, 100, 200, 500], `MinLeafSize` = [1, 3, 5, 10]
* 각 Fold 내에서 튜닝된 최적 Parameter로 최종 RF 모델을 Train 한 뒤 Test Set 평가.

---

## 6. Evaluation & Export
* **Performance Metrics**: RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), $R^2$ (R-squared).
* **Visualization Outputs**:
  1. MLR Sanity Check Scatter Plot (Overall Data)
  2. 5-Fold RF Predicted vs Actual Scatter Plot
  3. Feature Importance (Decision Tree 분기 기여도 기반)
  4. Error Distribution Histogram & Residual Plot
  5. MLR vs RF 성능 한눈에 보기 바 차트 (Fold별 RMSE 비교 및 $R^2$ 개선 요약)
  6. Cell 단위 Cycle 경과 SOH Trace Plot.
* **Workspace Export**: 시각화를 위해 저장된 Performance Metrics, K-Fold Prediction 결과, PCA Matrices 등을 `.mat` 혹은 `.xlsx` (Excel) 포맷으로 내보낸다.

위 Pipeline 및 Workflow는 수정사항 발생 시 지속적으로 업데이트되며, Lab Data를 기반으로 안정성이 확보되면 동일한 Feature Extraction & Modeling 방법론을 Field Data (`Field_Process` Directory)에 Deploy(횡개전)한다.

---

## 7. App Deployment UI 기획 (app1.mlapp 기반)
현재 구성된 `app1.mlapp`의 초기 UI 구조를 바탕으로 앱(App) 자동화 파이프라인의 레이아웃과 데이터 흐름을 정의 및 고도화한다.

### 7.1. 현재 구현된 UI 레이아웃 분석
* **좌측 패널 (실행 버튼)**: `DataLoadButton` -> `FeaturesLabelsButton` -> `ModelButton` -> `DeploymentButton` 순으로 파이프라인의 단계를 순차적으로 수행하는 흐름 타기.
* **중앙 패널 (옵션 선택 Tree)**: 
  * 데이터 선택 트리: RPT (Static Capacity, OCV, DCIR 등) 및 Drive Cycle (DC 1~8) 데이터 타입 지정.
  * 모델 선택 트리: 머신러닝/딥러닝 알고리즘 (Random Forest, SVM, LSTM, Transfer Learning) 종류 지정.
* **우측 패널 (상태 모니터링)**: `StatusLamp` (정상 동작 여부 피드백) 및 `Gauge` (RMSE 등 예측 성능 지표).

### 7.2. UI/UX 가독성 및 사용성 향상 제안 (추후 반영 권장사항)
사용자가 "데이터만 툭 던져주면 알아서 분석/시각화/모델링을 해주는" 목표를 달성하기 위해 UI를 다음과 같이 개선하는 것을 추천한다.

1. **파일 업로드 기능 강화 (사용자화)**: RPT 및 Drive Cycle 데이터의 csv(또는 mat) 파일 위치를 사용자가 직접 클릭해서 불러올 수 있도록 `EditField(텍스트 창)`와 `Select Folder/File (Browse)` 버튼 연결 추가 (`uigetfile` / `uigetdir` 함수 활용).
2. **논리적 구역(Panel) 분리**: 각 버튼과 트리를 기능별로 묶어(e.g., [Step 1: Data], [Step 2: Model]) `UIPanel` 안에 배치하면 시선의 흐름(좌에서 우, 혹은 위에서 아래)이 더 자연스럽고 그룹화되어 가독성이 높아진다.
3. **게이지(Gauge) UI의 현대적 대시보드화**: 투박한 선형 게이지(Linear Gauge) 대신, 성과 지표(RMSE, $R^2$)는 아주 커다란 **숫자 라벨(`uilabel`) 텍스트**와 함께 심플한 진행바(Progress bar) 혹은 반원형 게이지(Semicircular) 하나만 최소화하여 배치하면 훨씬 세련되고 전문적인 대시보드 룩을 완성할 수 있다.
4. **내장 그래프 영역 (UIAxes) 확보**: 우측의 큰 빈 공간을 활용해, 추출된 `dQ/dV` 곡선이나 학습된 SOH 예측 트렌드 그래프를 즉각적으로 확인할 수 있는 시각화(Plot) 영역을 앱 내부에 매립한다.
5. **실시간 로그 창 (Log Console) 추가**: 램프(Lamp) 통과/실패만으로는 내부에서 어떤 파일을 처리 중이고 어떤 에러가 났는지 알 수 없다. 하단에 `TextArea` (텍스트 영역)를 눕혀 배치하여 "XX 모델 학습 중..." 같은 실시간 메시지를 띄워준다.
6. **진행 상태바 (Progress Dialog)**: 데이터 로딩이나 모델 학습 동작 시, 실행 버튼 작동 후 앱 전체 화면에 반투명한 로딩 스피너(`uiprogressdlg`)가 나타나 현재 '열일하는 중' 임을 알게 해야 한다.

---

## 8. App Backend: Feature & Label 추출 상세

**Target Scripts:** `App_DataLoader.m`, `App_FeatureExtractor.m`, `process_single_rpt_file.m`

### 8.1. Data Loading & Pre-processing (`App_DataLoader.m`)
앱의 [Data Load] 단계에서 사용자가 선택한 폴더를 스캔하고, [Features & Labels] 단계에서 실제 전처리를 수행한다.

| Step | Description | Function |
|---|---|---|
| 1 | CSV 파일 스캔 (`*_RPT_*.csv`) | `dir()` in callback |
| 2 | 파일별 파싱 (textscan) + Moving Average (Window=20) | `process_single_rpt_file()` |
| 3 | 0.001V Uniform Grid 선형 보간 | `interp1()` inside helper |

**출력 구조체** `App_VQ_grid`:
```
App_VQ_grid
  └── cyc0 / cyc200 / ... / cyc1000
        └── Ch09 / Ch10 / ... / Ch16
              ├── Static        (V_grid, Q, t)
              ├── OCV_charge    (V_grid, Q, t)
              └── OCV_discharge (V_grid, Q, t)
```

### 8.2. Feature Extraction (`App_FeatureExtractor.m`)

#### Features (14개) — OCV Charge/Discharge 데이터에서 추출

| # | Feature | Source | Method |
|---|---|---|---|
| 1~5 | `Chg_dQ(1-5)` | `OCV_charge` | MasterRuler `V_bounds_chg` 경계에서 $\|\Delta Q\|$ |
| 6~10 | `Dch_dQ(1-5)` | `OCV_discharge` | MasterRuler `V_bounds_dch` 경계에서 $\|\Delta Q\|$ |
| 11 | `Chg_PkH` | `OCV_charge` | dQ/dV Peak Height [3.70~3.95V] |
| 12 | `Chg_PkA` | `OCV_charge` | dQ/dV Peak Area [3.70~3.95V] |
| 13 | `Dch_PkH` | `OCV_discharge` | dQ/dV Peak Height [3.75~3.88V] |
| 14 | `Dch_PkA` | `OCV_discharge` | dQ/dV Peak Area [3.75~3.88V] |

#### Labels (3개) — Static + OCV 데이터에서 추출

| Label | Source | Method |
|---|---|---|
| **SOH (Ah)** | `Static.Q` | `max(abs(Q))` — Static Capacity 최대 용량 |
| **LLI (%)** | `OCV_charge` vs `Fresh_OCV_Charge` | ICA Peak 위치 이동: $\text{LLI} = \frac{\|V_{peak,aged} - V_{peak,fresh}\| \times H_{peak,fresh}}{Q_{rated}} \times 100$ |
| **LAM (%)** | `OCV_charge` vs `Fresh_OCV_Charge` | ICA Peak 면적 비율: $\text{LAM} = (1 - \frac{A_{peak,aged}}{A_{peak,fresh}}) \times 100$ [3.4~4.0V] |

<details>
<summary><b>💡 상세 설명: Master Ruler와 Fresh_OCV_Charge의 역할 분리 및 LLI/LAM 추출 과정</b></summary>

이 파이프라인에서 `MasterRulers.mat` 파일은 사실 **두 가지 완전히 다른 목적**을 위해 사용됩니다.

1. **X 피처 추출용 룰러 (Static Discharge 기반):** `V_bounds_chg`, `V_bounds_dch` 
   - 5개의 C-rate(0.1C~3C) V-Q 곡선을 동일한 5개 구간(Segment)으로 자르기 위한 **절대적인 전압 기준선**입니다. 배터리가 낼 수 있는 최대 전압 구간을 잡기 위해 가장 느린 **Static Capacity 방전 곡선**을 평균 내어 생성했습니다.
2. **Y 라벨 계산용 기준점 (OCV Charge 기반):** `Fresh_OCV_Charge`
   - SOH(용량)는 현재 사이클의 최대 Q값만 보면 되지만, **LLI(리튬 손실)와 LAM(활물질 손실)**은 현재 값이 절대적으로 얼마인지가 아니라 **"최초 상태(0 사이클) 대비 얼마나 열화되었는가?"**를 비교해야만 알 수 있습니다.
   - 하지만 400, 800 사이클 등의 데이터를 개별적으로 로드하여 분석할 때 메모리에는 0 사이클 데이터가 남아있지 않으므로, 아예 0 사이클의 OCV 충전 곡선(`OCV_charge`) 통째로 `MasterRulers.mat`에 `Fresh_OCV_Charge`라는 이름으로 영구 보관해두고 두고두고 "기준점(Anchor)"으로 꺼내 쓰는 것입니다.

**ICA (Incremental Capacity Analysis) 기반 LLI, LAM 추출 알고리즘**
*   **`get_main_peak(data, minV, maxV)`**:
    - OCV 곡선을 미분하여 $dQ/dV$를 계산합니다. 미분은 노이즈를 극도로 증폭시키므로, 반드시 **크기 21의 이동평균(Moving Average)** 필터를 적용하여 곡선을 부드럽게 깎아냅니다.
    - 이후 주요 화학 반응이 일어나는 `[3.40V ~ 4.00V]` 구간 내에서 가장 우뚝 솟은 메인 피크의 전압 위치($V_{peak}$), 높이($H_{peak}$), 그리고 면적($A_{peak}$)을 찾아 반환합니다.
*   **`analyze_ica_aging(current_ocv, fresh_ocv, Q_rated)`**:
    - `Fresh` 상태의 피크와 노화된 `Aged` 상태의 피크를 각각 `get_main_peak()`로 찾아 비교합니다.
    - **LLI (%)**: $dQ/dV$ 피크가 0 사이클 대비 오른쪽(더 높은 전압)으로 얼마나 밀려났는지(**Peak Shift**)를 측정합니다. 리튬 이온이 고갈되고 내부 저항이 증가할수록 피크가 지연되어 나타나는 현상을 정량화한 것입니다. 피크 이동량($\Delta V$)에 초기 피크 높이($H_{fresh}$)를 곱하고 정격 용량($Q_{rated}$)으로 나누어 퍼센트로 환산합니다.
    - **LAM (%)**: $dQ/dV$ 피크의 면적이 0 사이클 대비 얼마나 쪼그라들었는지(**Peak Shrinkage**)를 측정합니다. 면적이 줄어든다는 것은 해당 전압 대역에서 리튬을 수용할 수 있는 활물질(Active Material) 구조 자체가 무너져(Loss) 사라졌음을 의미합니다. (1 - 면적 비율)을 통해 손실률을 퍼센트로 환산합니다.
</details>

### 8.3. 최종 출력 (FeatureTable)

```
FeatureTable (table)
  ├── CellID          : string    ('Ch09', 'Ch10', ...)
  ├── Cycle           : double    (0, 200, 400, ...)
  ├── X_Features      : [N×14]    (Raw Features)
  ├── Y_Labels        : [N×3]     (SOH, LLI, LAM)
  └── X_Normalized    : [N×14]    (Z-score Normalized Features)
```

---

## 9. App Visualization Plan

현재 앱에는 각 단계(Data Load, Features & Labels, Model)마다 **Visualization 체크박스**가 있다.
체크된 상태에서 해당 단계를 실행하면, 분석 결과를 **외부 figure 창**으로 팝업하여 시각화한다.

**참조 스크립트 목록:**
| Script | 역할 |
|---|---|
| `RPT_Pipeline_Visualization.m` | 전처리 과정 시각화 (Raw vs Interp, V-Q overlay, Voltage Window, Master Ruler) |
| `RPT_Label_Visualization.m` | Label 추이 시각화 (SOH/LLI/LAM trend, OCV dQ/dV Peak Evolution) |
| `RPT_Result_Visualization_NoPCA.m` | 모델 결과 시각화 (R², Pred vs Actual, Error, Feature Importance, Residual) |

---

### 9.1. Data Load → Visualization (`VisualizationCheckBox`)
**대상 데이터**: `App_VQ_grid` (전처리 완료된 V-Q 구조체)

Data Load 완료 후 불러온 데이터를 확인하기 위한 시각화.

| # | Figure | 참조 스크립트 | 내용 |
|---|---|---|---|
| 1 | **V-Q Overlay (Charge/Discharge)** | `RPT_Pipeline_Visualization.m` Phase 0-2 | 전 채널 × 전 사이클 충/방전 V-Q 곡선 (2-subplot) |
| 2 | **Voltage Window 표시** | `RPT_Pipeline_Visualization.m` Phase 0-3 | 대표 채널 V-Q + 분석 윈도우(Chg 3.70~3.95V, Dch 3.75~3.88V) 음영 |
| 3 | **SOH Trend (Static Capacity)** | `RPT_Pipeline_Visualization.m` Phase 0-4 | 사이클별 Static 용량 변화 (Ah + %, 8채널) |

---

### 9.2. Features & Labels → Visualization (`VisualizationCheckBox_2`)
**대상 데이터**: `FeatureTable` (14 Features + 3 Labels)

추출된 피처와 라벨을 확인하기 위한 시각화. Tree_3 체크 항목에 따라 선택적 표시.

**Tree_3 노드별 선택 시각화:**

| Tree_3 Node | Figure | 참조 스크립트 | 내용 |
|---|---|---|---|
| **Equilibrium** | OCV dQ/dV Peak Evolution | `RPT_Label_Visualization.m` §3 | 채널별 8-subplot, 사이클별 dQ/dV overlay (jet colormap) |
| **Label → Available Capacity** | SOH Trend | `RPT_Label_Visualization.m` §2 | X=Cycle, Y=SOH(Ah), 8채널 라인 |
| **Label → LLI** | LLI Trend | `RPT_Label_Visualization.m` §2 | X=Cycle, Y=LLI(%), 8채널 라인 |
| **Label → LAMp** | LAM Trend | `RPT_Label_Visualization.m` §2 | X=Cycle, Y=LAM(%), 8채널 라인 |

**기본 시각화 (항상 표시):**

| # | Figure | 참조 스크립트 | 내용 |
|---|---|---|---|
| 1 | **Feature Heatmap** | `RPT_Pipeline_Visualization.m` Phase 2-3 | 14개 Feature × N samples 히트맵 (정규화 후) |
| 2 | **Label Summary (3-subplot)** | `RPT_Label_Visualization.m` §2 | SOH / LLI / LAM vs Cycle (전채널 한 번에) |

---

### 9.3. Model → Visualization (`VisualizationCheckBox_3`)
**트리거**: `ModelButton` 완료 후, 체크 시 팝업.

| # | Figure | 참조 스크립트 | 내용 |
|---|---|---|---|
| 1 | **R² Bar Chart** | `RPT_Result_Visualization_NoPCA.m` §2 | SOH/LLI/LAM 별 R² 점수 바 차트 |
| 2 | **Predicted vs Actual** | `RPT_Result_Visualization_NoPCA.m` §3 | 3-subplot 산점도 + 45° 기준선 + R² 표시 |
| 3 | **Error Distribution** | `RPT_Result_Visualization_NoPCA.m` §4 | 3-subplot 히스토그램 + RMSE 값 |
| 4 | **Feature Importance** | `RPT_Result_Visualization_NoPCA.m` §5 | 3-subplot 수평 바 차트 (14개 Feature 기여도 순위) |
| 5 | **Convergence Plot** | `RPT_Result_Visualization_NoPCA.m` §5.5 | NumTrees × MinLeafSize별 OOB RMSE 수렴 곡선 |
| 6 | **Residual Analysis** | `RPT_Result_Visualization_NoPCA.m` §6 | Residual vs Actual 산점도 (편향 진단) |

---

### 9.4. 구현 방식 요약

```matlab
% FeaturesLabelsButton 콜백 끝에 추가:
if app.VisualizationCheckBox_2.Value
    App_Visualizer_Features(app.ProcessData_VQ, app.FeatureTable, MasterRulers);
end

% ModelButton 콜백 끝에 추가:
if app.VisualizationCheckBox_3.Value
    App_Visualizer_Model(modelResults);
end
```

* 시각화 함수를 `App_Visualizer_Features.m`, `App_Visualizer_Model.m`으로 별도 스크립트 작성
* 각 figure는 `saveas(fig, path, 'fig')` + `exportgraphics(fig, path, 'pdf')`로 `App/Results/` 폴더에 자동 저장

