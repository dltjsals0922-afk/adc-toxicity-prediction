# Multi-task Machine Learning Pipeline for Predicting Cytotoxic Potency and Selectivity of ADC Payload-like Compounds from ChEMBL Assay Data

[![Python](https://img.shields.io/badge/Python-3.10-blue)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2023-green)](https://www.rdkit.org/)
[![XGBoost](https://img.shields.io/badge/XGBoost-3.x-orange)](https://xgboost.readthedocs.io/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0-red)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 프로젝트 개요

ChEMBL의 세포독성(cytotoxicity) assay 데이터를 기반으로  
ADC(Antibody-Drug Conjugate) 페이로드 후보 화합물의 독성을 예측하는  
Multi-task Machine Learning 파이프라인.

### 연구 질문
> "분자 구조 피처(Morgan Fingerprint + RDKit Descriptors)만으로  
>  ADC 페이로드 후보의 세포독성을 예측할 수 있는가?"

---

## 배경

ADC는 항체(Antibody), 링커(Linker), 페이로드(Payload)로 구성된 항암제다.  
페이로드는 강력한 세포독성 물질(MMAE, DM1, SN-38 등)로,  
후보 스크리닝 단계에서 독성 예측 모델이 있으면 비용과 시간을 절감할 수 있다.

---

## 데이터

| 항목 | 내용 |
|------|------|
| 출처 | ChEMBL (EMBL-EBI) |
| 수집 방법 | ChEMBL REST API |
| Assay 조건 | assay_type=F, standard_type=IC50, cytotoxicity 키워드 |
| 원본 데이터 | 2,000개 |
| 정제 후 | 1,758개 |
| 고유 화합물 | 1,111개 |
| 고유 Assay | 250개 |

### 라벨링 기준

| Task | 유형 | 기준 |
|------|------|------|
| Task 1 | pIC50 회귀 | pIC50 = -log10(IC50[M]) |
| Task 2 | 독성 분류 | IC50 < 100nM → 고독성(1) |
| Task 3 | 등급 분류 | <50nM(high) / 50~500nM(mid) / >500nM(low) |

---

## 성능

### Task 2 (독성 이진 분류) AUC-ROC

| 모델 | AUC-ROC | 비고 |
|------|---------|------|
| Logistic Regression | 0.8963 | Morgan FP + Descriptors |
| Random Forest | 0.9298 | Morgan FP + Descriptors |
| XGBoost | 0.9437 | Morgan FP + Descriptors (2058 피처) |
| GCN (v3.0) | 0.9280 | 원자 피처 6개만 |
| GAT (v3.0) | 0.9416 | 원자 피처 6개 + 어텐션 |

### Task 1 / Task 3

| Task | 지표 | 값 |
|------|------|-----|
| Task 1 (pIC50 회귀) | R² | 0.747 |
| Task 3 (등급 분류) | AUC-ROC (macro) | 0.883 |

### GNN 핵심 인사이트

피처 엔지니어링 없이 분자 그래프 구조만으로 XGBoost(0.9437)에 근접한 GAT(0.9416) 달성.
GAT 어텐션 메커니즘이 GCN보다 우수 (0.9280 → 0.9416).
향후 노드 피처 확장 시 XGBoost 초과 가능.

---

## 주요 발견 (SHAP 분석)

### 1. LogP × TPSA 조합이 핵심 예측 인자

단독 피처보다 조합이 2.5배 강한 분리력을 보임:

| 조건 | 고독성 비율 |
|------|------------|
| LogP < 4 AND TPSA < 100 | 9.6% (가장 안전) |
| LogP > 4 AND TPSA > 100 | 40.0% (가장 위험) |
| 전체 평균 | 27.8% |

### 2. TPSA는 MolWt의 proxy (r = 0.922)

TPSA와 MolWt의 상관계수가 0.922. TPSA 단독 인과 해석 주의.

### 3. TPSA 100~140 구간이 고독성 피크 (39.2%)

ADC 페이로드 최적 분자량 범위(MolWt ~487)와 일치.
endocytosis 메커니즘으로 LogP와 무관하게 독성 가능.

---

## v2.0 NCI-60 구조-활성 관계 (SAR) 분석

### 데이터
- NCI-60 (CellMiner): 58,953개 화합물 × 60개 세포주 pGI50
- 선택적 고독성 화합물: 141개 (pGI50 > 6 AND 선택성 지수 > 2)

### 핵심 발견: Scaffold → 암종 선택성

| Scaffold | 주요 암종 | 선택성 | 메커니즘 |
|----------|----------|--------|---------|
| Anilinopyrimidine | 흑색종 (ME) | 80% | MEK/BRAF 억제 |
| Antifolate | 전립선암 (PR) | 55% | DHFR 억제 |
| Alkaloid | 백혈병 (LE) | 67% | 리보솜/위상이성질화효소 억제 |
| Macrolactone | 유방암 (BR) | 40% | 미세소관 억제 |

### ADC 설계 시사점

유방암 타깃 ADC → Macrolactone 계열 페이로드 우선 검토  
흑색종 타깃 ADC → Anilinopyrimidine 계열  
전립선암 타깃 ADC → Antifolate 계열  
백혈병 타깃 ADC → Alkaloid 계열

---

## 한계 및 향후 과제

### 한계
- 다양한 세포주/assay 조건의 데이터가 혼재 → assay 간 노이즈 존재
- 상관관계 ≠ 인과관계 (데이터셋 선택 편향 가능성)
- GI50 ≠ IC50 → 외부 검증 시 assay 이질성 존재
- 정상세포 데이터 없음 → 진정한 selectivity 측정 불가

### 향후 과제
- v4.0: 정상세포 IC50 데이터 통합 → 진정한 selectivity 계산
- v5.0: 노드 피처 확장 → GAT 성능 향상

---

## 폴더 구조

```
adc-toxicity-prediction/
├── data/
│   ├── raw/
│   └── processed/
├── notebooks/
│   ├── 01_data_collection.ipynb
│   ├── 02_eda.ipynb
│   ├── 03_feature_engineering.ipynb
│   ├── 04_modeling.ipynb
│   ├── 05_shap_deep.ipynb
│   ├── 06_nci60_selectivity.ipynb
│   └── 07_gnn.ipynb
├── src/
│   └── fetch_chembl.py
├── environment.yml
├── requirements.txt
└── README.md
```

---

## 환경 설정

```bash
conda env create -f environment.yml
conda activate adc-env

# PyTorch + PyG (GNN용)
pip install torch==2.0.1 --index-url https://download.pytorch.org/whl/cpu
pip install torch-geometric
pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.0.1+cpu.html
```

---

## 기술 스택

| 분류 | 라이브러리 |
|------|-----------|
| 화학정보학 | RDKit, chembl-webresource-client |
| 머신러닝 | scikit-learn, XGBoost |
| 딥러닝 | PyTorch, PyTorch Geometric (GCN, GAT) |
| 해석 | SHAP |
| 데이터 | pandas, numpy |
| 시각화 | matplotlib, seaborn |

---

## 참고문헌

- Petersson, E.J. et al. (2021). ADC payload potency and linker design. *J. Med. Chem.*
- Lundberg, S.M. & Lee, S.I. (2017). A unified approach to interpreting model predictions. *NeurIPS*
- Veličković, P. et al. (2018). Graph Attention Networks. *ICLR*

---

## 라이선스

MIT License