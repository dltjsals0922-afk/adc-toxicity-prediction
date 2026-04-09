# Multi-task Machine Learning Pipeline for Predicting Cytotoxic Potency and Selectivity of ADC Payload-like Compounds from ChEMBL Assay Data

[![Python](https://img.shields.io/badge/Python-3.10-blue)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2023-green)](https://www.rdkit.org/)
[![XGBoost](https://img.shields.io/badge/XGBoost-3.x-orange)](https://xgboost.readthedocs.io/)
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

```
항체 (Antibody)  —  링커 (Linker)  —  페이로드 (Payload)
 표적 인식                            세포독성 물질
```

---

## 데이터

| 항목 | 내용 |
|------|------|
| 출처 | ChEMBL (EMBL-EBI) |
| 수집 방법 | ChEMBL REST API (`chembl-webresource-client`) |
| Assay 조건 | `assay_type=F`, `standard_type=IC50`, `cytotoxicity` 키워드 |
| 원본 데이터 | 2,000개 |
| 정제 후 | 1,758개 |
| 고유 화합물 | 1,111개 |
| 고유 Assay | 250개 |

### 라벨링 기준

| Task | 유형 | 기준 |
|------|------|------|
| Task 1 | pIC50 회귀 | `pIC50 = -log10(IC50[M])` |
| Task 2 | 독성 분류 | IC50 < 100nM → 고독성(1) |
| Task 3 | 등급 분류 | <50nM(high) / 50~500nM(mid) / >500nM(low) |

> Task 3 기준 근거: ADC 페이로드 문헌 기준(Petersson et al., 2021) + 데이터 분포 중앙값(990nM) 고려

---

## 모델 구조

```
SMILES
  ↓
RDKit → Morgan Fingerprint (2048bit, radius=2)
      + RDKit Descriptors (10개: LogP, MolWt, HBD, HBA, TPSA 등)
  ↓
Feature matrix (2058개 피처)
  ↓
┌─────────────────────────────────────┐
│  Baseline 1: Logistic Regression    │
│  Baseline 2: Random Forest          │
│  Proposed:   XGBoost (Multi-output) │
│    ├─ Task 1: pIC50 회귀            │
│    ├─ Task 2: 독성 이진 분류        │
│    └─ Task 3: 독성 등급 분류        │
└─────────────────────────────────────┘
  ↓
SHAP 해석
```

---

## 성능

### Task 2 (독성 이진 분류) AUC-ROC

| 모델 | AUC-ROC |
|------|---------|
| Baseline 1: Logistic Regression | 0.8963 |
| Baseline 2: Random Forest | 0.9298 |
| **Proposed: XGBoost** | **0.9437** |

### Task 1 / Task 3

| Task | 지표 | 값 |
|------|------|-----|
| Task 1 (pIC50 회귀) | R² | 0.747 |
| Task 3 (등급 분류) | AUC-ROC (macro) | 0.883 |

---

## 주요 발견 (SHAP 분석)

### 1. LogP × TPSA 조합이 핵심 예측 인자

단독 피처보다 조합이 **2.5배 강한 분리력**을 보임:

| 조건 | 고독성 비율 |
|------|------------|
| LogP < 4 AND TPSA < 100 | **9.6%** (가장 안전) |
| LogP > 4 AND TPSA > 100 | **40.0%** (가장 위험) |
| 전체 평균 | 27.8% |

### 2. TPSA는 MolWt의 proxy (r = 0.922)

TPSA와 MolWt의 상관계수가 0.922로 거의 동일한 정보를 담고 있음.  
→ TPSA 단독 인과 해석 주의.

### 3. TPSA 100~140 구간이 고독성 피크 (39.2%)

ADC 페이로드 최적 분자량 범위(MolWt ~487)와 일치.  
→ endocytosis 메커니즘으로 LogP와 무관하게 독성 가능.

### 4. Morgan FP 서브구조가 예측 주도

Top SHAP 피처 대부분이 Morgan FP 비트.  
→ 분자 구조 패턴이 단순 물리화학적 특성보다 중요.

---

## 한계 및 향후 과제

### 한계
- 다양한 세포주/assay 조건의 데이터가 혼재 → assay 간 노이즈 존재
- 상관관계 ≠ 인과관계 (데이터셋 선택 편향 가능성)
- 역설 케이스(LogP<1 고독성) n=6으로 통계적 신뢰도 제한

### 향후 과제
```
v2.0: NCI-60 데이터 추가 → 외부 검증 세트 구성
v3.0: PyTorch Geometric GNN → 분자 그래프 구조 학습
      "GNN이 복잡한 LogP × TPSA 상호작용을 학습할 것으로 기대"
```

---

## 폴더 구조

```
adc-toxicity-prediction/
├── data/
│   ├── raw/                        # ChEMBL 원본 데이터
│   └── processed/                  # 정제 + 피처 포함 데이터
├── notebooks/
│   ├── 01_data_collection.ipynb    # ChEMBL API 데이터 수집
│   ├── 02_eda.ipynb                # 탐색적 데이터 분석
│   ├── 03_feature_engineering.ipynb # RDKit 피처 추출
│   ├── 04_modeling.ipynb           # 모델 학습 + 성능 평가
│   └── 05_shap_deep.ipynb          # SHAP 심층 분석
├── src/
│   └── fetch_chembl.py             # ChEMBL 데이터 수집 모듈
├── environment.yml                 # conda 환경 설정
├── requirements.txt                # pip 의존성
└── README.md
```

---

## 환경 설정

```bash
# conda 환경 생성 (권장)
conda env create -f environment.yml
conda activate adc-env

# 또는 pip
pip install -r requirements.txt
```

---

## 실행 순서

```bash
# 1. 데이터 수집
jupyter nbconvert --to notebook --execute notebooks/01_data_collection.ipynb

# 2. 순서대로 실행
# notebooks/ 폴더의 01 → 02 → 03 → 04 → 05 순서
```

---

## 기술 스택

| 분류 | 라이브러리 |
|------|-----------|
| 화학정보학 | RDKit, chembl-webresource-client |
| 머신러닝 | scikit-learn, XGBoost |
| 해석 | SHAP |
| 데이터 | pandas, numpy |
| 시각화 | matplotlib, seaborn |

---

## 참고문헌

- Bray, M.A. et al. (2016). Cell Painting, a high-content image-based assay for morphological profiling using multiplexed fluorescent dyes. *Nature Protocols*
- Petersson, E.J. et al. (2021). ADC payload potency and linker design. *J. Med. Chem.*
- Lundberg, S.M. & Lee, S.I. (2017). A unified approach to interpreting model predictions. *NeurIPS*

---

## 라이선스

MIT License
