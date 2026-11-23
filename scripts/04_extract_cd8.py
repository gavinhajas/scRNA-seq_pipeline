import sys
from extract_common import run

DATASETS = ["GSE149614", "GSE156625", "GSE140228"]

if __name__ == "__main__":
    targets = [sys.argv[1]] if len(sys.argv) > 1 else DATASETS
    for ds in targets:
        run(ds, "cd8")