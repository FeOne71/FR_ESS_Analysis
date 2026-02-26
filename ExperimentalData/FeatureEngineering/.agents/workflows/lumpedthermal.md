---
description: Degradation model using lumped thermal model
---

1. FeatureEngineering 폴더의 하위폴더 Reference의 3가지 논문을 참고하여 열화모델 개발
2. 랩데이터 구조는 Lab_RPT_Anylsis 폴더의 하위폴더 App 폴더안의 스크립트 참조. 랩데이터는 8개 셀 실험 데이터임.
3. 필드데이터 구조는 Reference폴더의 하위폴더 Event_Detection.m 참조.
4. 논문과 같이 랩데이터를 사용해서 배터리 셀 모델 (electrical, thermal)을 만들고 필드데이터에 시연해보는 것이 목적.
5. 모든 작업은 FeatureEngineering 폴더에 하위폴더 "Lumped_Model_0224" 폴더를 만들고 해당 폴더안에서 작업할것
6. Matlab R2025b 버전으로 진행. 시각화는 반드시 fig로만 저장. png 사용 금지.
7. implentation plan을 Lumped_Model_0224 폴더에 생성 및 저장. 코드 수정 및 방향 수정이 되었을 경우에 업데이트 할것.
8. 플랜은 항상 한국어로 (용어제외)
9. 수식도 함께 작성할것. latex 수식이 깨지지는 일이 자주 발생하므로 깨지지 않게 할것.

*** Lab Data ***
1. RPT는 Static capacity (0.5C) - OCV Chg(0.05C) - OCV Dchg(0.05C) - SOC-DCIR - C-rate (0.1C, 0.5C, 1C, 2C, 3C) 로 구성됨
2. RPT 후처리 스크립트: D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\RPT_Postprocessing.m
3. RPT 후처리 데이터 (V,Q,C-rate등등): D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat
4. RPT 후처리 데이터 (ocv 데이터): D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat
5. RPT 후처리 데이터 (DCIR):
D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Archive\RPT_DCIR_Final.m


*** Field Data ***
1. 필드데이터 csv to mat 스크립트가 2개 있는데 뭐가 다른지 기억이 안남.
D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New_rack2mat.m
D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\raw2mat_ver05.m
2. 필드데이터에 HVAC 팬이 4개 있는데 뭐를 사용해야할지 모르겠음. 그리고 HVAC Setting Temp도 4개가 존재.
3. 필드데이터 mat 저장위치: D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New, D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old
2개가 존재. Old, New 폴더안에는 연도별 폴더 YYYY가 있고, 연도별 폴더 안에는 월별폴더가 존재. YYYYMM. 월별 폴더 안에 일자별 mat 파일 존재.