"""
fetch_chembl.py
ChEMBL에서 ADC 페이로드 후보 화합물의 세포독성 데이터 수집

수집 전략:
  A) cytotoxicity assay 키워드 검색 (넓은 그물)
  B) 알려진 ADC 페이로드 계열 직접 검색 (MMAE, DM1, SN-38 등)
  두 결과를 합쳐서 중복 제거 후 저장
"""

import pandas as pd
from tqdm import tqdm
from chembl_webresource_client.new_client import new_client

# ── ChEMBL 클라이언트 ──────────────────────────────────────────────
activity = new_client.activity
assay    = new_client.assay
molecule = new_client.molecule


# ──────────────────────────────────────────────────────────────────
# 전략 A: cytotoxicity assay 검색
# ──────────────────────────────────────────────────────────────────
def fetch_by_assay_keyword(keyword: str = "cytotoxicity",
                           standard_type: str = "IC50",
                           max_assays: int = 50) -> pd.DataFrame:
    """
    assay description에 keyword가 포함된 assay를 찾고
    해당 assay의 IC50 활성 데이터를 수집한다.

    Parameters
    ----------
    keyword      : assay 검색 키워드
    standard_type: 활성 타입 (IC50 고정)
    max_assays   : 검색할 assay 최대 수 (API 부하 조절)

    Returns
    -------
    DataFrame: molecule_chembl_id, smiles, IC50(nM), assay_description
    """
    print(f"[A] assay 검색 중: '{keyword}'")

    assays = list(
        assay.filter(
            description__icontains=keyword,
            assay_type="B"          # Binding/세포 기반 assay
        ).only(["assay_chembl_id", "description"])[:max_assays]
    )
    print(f"    → {len(assays)}개 assay 발견")

    records = []
    for a in tqdm(assays, desc="    assay 순회"):
        acts = list(
            activity.filter(
                assay_chembl_id=a["assay_chembl_id"],
                standard_type=standard_type,
                standard_units="nM",
            ).only([
                "molecule_chembl_id",
                "canonical_smiles",
                "standard_value",
                "standard_units",
                "assay_chembl_id",
            ])
        )
        for act in acts:
            if act.get("standard_value") and act.get("canonical_smiles"):
                records.append({
                    "molecule_chembl_id": act["molecule_chembl_id"],
                    "smiles":             act["canonical_smiles"],
                    "ic50_nM":            float(act["standard_value"]),
                    "assay_id":           act["assay_chembl_id"],
                    "assay_description":  a["description"],
                    "source":             "assay_keyword",
                })

    df = pd.DataFrame(records)
    print(f"    → 총 {len(df)}개 활성 데이터 수집\n")
    return df


# ──────────────────────────────────────────────────────────────────
# 전략 B: 알려진 ADC 페이로드 계열 직접 검색
# ──────────────────────────────────────────────────────────────────
ADC_PAYLOAD_KEYWORDS = [
    "monomethyl auristatin",   # MMAE / MMAF
    "auristatin",
    "maytansine",              # DM1 / DM4
    "emtansine",
    "calicheamicin",
    "SN-38",
    "dxd",                     # Dxd (Enhertu 페이로드)
    "pyrrolobenzodiazepine",   # PBD
    "duocarmycin",
    "camptothecin",
]

def fetch_by_payload_keywords(standard_type: str = "IC50") -> pd.DataFrame:
    """
    ADC 페이로드로 알려진 화합물 계열을 이름 기반으로 검색하고
    해당 화합물들의 IC50 데이터를 수집한다.

    Returns
    -------
    DataFrame: molecule_chembl_id, smiles, IC50(nM), pref_name
    """
    print("[B] ADC 페이로드 계열 검색")
    records = []

    for kw in tqdm(ADC_PAYLOAD_KEYWORDS, desc="    페이로드 키워드"):
        mols = list(
            molecule.filter(
                pref_name__icontains=kw
            ).only(["molecule_chembl_id", "pref_name",
                    "molecule_structures"])
        )

        for mol in mols:
            chembl_id = mol.get("molecule_chembl_id")
            pref_name = mol.get("pref_name", "")
            smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles")

            if not smiles:
                continue

            acts = list(
                activity.filter(
                    molecule_chembl_id=chembl_id,
                    standard_type=standard_type,
                    standard_units="nM",
                ).only(["standard_value", "assay_chembl_id"])
            )

            for act in acts:
                if act.get("standard_value"):
                    records.append({
                        "molecule_chembl_id": chembl_id,
                        "smiles":             smiles,
                        "ic50_nM":            float(act["standard_value"]),
                        "assay_id":           act.get("assay_chembl_id", ""),
                        "assay_description":  pref_name,
                        "source":             "payload_keyword",
                    })

    df = pd.DataFrame(records)
    print(f"    → 총 {len(df)}개 활성 데이터 수집\n")
    return df


# ──────────────────────────────────────────────────────────────────
# 데이터 정제 + 저장
# ──────────────────────────────────────────────────────────────────
def clean_and_save(df: pd.DataFrame, save_path: str) -> pd.DataFrame:
    """
    두 전략의 데이터를 합치고 정제한 뒤 CSV로 저장한다.

    정제 기준:
      - IC50 값이 0 이하이거나 10,000,000 nM 초과인 것 제거 (명백한 이상값)
      - SMILES 결측 제거
      - 같은 화합물 + assay 중복 제거
      - 이진 라벨 추가: ic50_nM < 100 → label=1 (고독성)
    """
    print("[정제] 데이터 클리닝 시작")
    print(f"    원본 행 수: {len(df)}")

    # 결측 및 이상값 제거
    df = df.dropna(subset=["smiles", "ic50_nM"])
    df = df[(df["ic50_nM"] > 0) & (df["ic50_nM"] <= 1_000_000)]

    # 중복 제거 (같은 화합물 + 같은 assay)
    df = df.drop_duplicates(subset=["molecule_chembl_id", "assay_id"])

    # 이진 라벨링 (기준: 100 nM)
    df["label"] = (df["ic50_nM"] < 100).astype(int)

    print(f"    정제 후 행 수: {len(df)}")
    print(f"    고독성(1): {df['label'].sum()}개 / 저독성(0): {(df['label']==0).sum()}개")

    df.to_csv(save_path, index=False)
    print(f"    저장 완료: {save_path}\n")
    return df


# ──────────────────────────────────────────────────────────────────
# 메인 실행
# ──────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import os

    RAW_PATH = "data/raw/chembl_raw.csv"
    os.makedirs("data/raw", exist_ok=True)

    # 두 전략으로 수집
    df_a = fetch_by_assay_keyword(keyword="cytotoxicity", max_assays=30)
    df_b = fetch_by_payload_keywords()

    # 합치기
    df_all = pd.concat([df_a, df_b], ignore_index=True)

    # 정제 + 저장
    df_clean = clean_and_save(df_all, RAW_PATH)

    print("=" * 50)
    print(df_clean.head())
    print(f"\n컬럼: {list(df_clean.columns)}")
