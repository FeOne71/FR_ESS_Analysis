---
description: Field Data Processing and Analysis Pipeline
---

### 1. Data Ingestion & Preprocessing
- Load year-indexed raw data (e.g., `Raw_20210603.mat`).
- Identify Charge/Discharge segments using a C-rate threshold (e.g., > 0.05C / Np).
- 

### 2. Feature Extraction (14-Feature Matrix)
- 너의 지식을 활용해서 저항성분, 에너지효율, 용량 등 raw.mat의 변수를 사용하던 새롭게 추출하던 feature과 label 후보를 생성하고 연도별로 어떻게 진행되는지 확인할것. 즉, 피쳐와 라벨을 추출하는작업을 할것
- 

### 3. State Estimation (SOH, Energy Efficiency, Power (LLI, LAM if possible))
- Random forest, Transfer Learning, LSTM model

### 4. Visualization & Documentation
- Plot SOH, LLI, and LAM trends over time.
- Save summary figures in `.fig` and `.pdf` formats.
- Log estimation results for further analysis.