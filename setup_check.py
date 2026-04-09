import sys

REQUIRED = [
    ("numpy",                     "numpy"),
    ("pandas",                    "pandas"),
    ("rdkit",                     "rdkit"),
    ("chembl_webresource_client", "chembl-webresource-client"),
    ("sklearn",                   "scikit-learn"),
    ("xgboost",                   "xgboost"),
    ("imblearn",                  "imbalanced-learn"),
    ("shap",                      "shap"),
    ("matplotlib",                "matplotlib"),
    ("seaborn",                   "seaborn"),
]

print(f"Python {sys.version}\n")
print(f"{'패키지':<32} {'상태':<10} {'버전'}")
print("-" * 60)

all_ok = True
for import_name, pip_name in REQUIRED:
    try:
        mod = __import__(import_name)
        version = getattr(mod, "__version__", "n/a")
        print(f"  {pip_name:<30} {'OK':<10} {version}")
    except ImportError:
        print(f"  {pip_name:<30} {'없음':<10}")
        all_ok = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    fp  = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    mw  = Descriptors.MolWt(mol)
    print(f"\nRDKit 기능 테스트: OK (아스피린 MW={mw:.2f}, FP={len(fp)}bits)")
except Exception as e:
    print(f"\nRDKit 기능 테스트: 실패 ({e})")
    all_ok = False

try:
    from chembl_webresource_client.new_client import new_client
    result = list(new_client.molecule.filter(chembl_id="CHEMBL25").only(["chembl_id"]))
    print("ChEMBL API 연결:   OK")
except Exception as e:
    print(f"ChEMBL API 연결:   실패 ({e})")
    all_ok = False

print("\n" + "=" * 60)
print("환경 세팅 완료! 2일차로 넘어가자." if all_ok else "위 실패 항목 확인 필요.")
print("=" * 60)
