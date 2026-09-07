# Visium HD count-accumulation validation

Supporting evidence for [HEST #143](https://github.com/mahmoodlab/HEST/issues/143). This folder lives on a fork-only evidence branch; the upstream fix consists of the pooling change and its regression test, not this dataset-specific validation bundle.

Upstream pull request: [Fix count loss when pooling Visium HD bins #145](https://github.com/mahmoodlab/HEST/pull/145).

- [Runnable validation script](validate_real_data.py)
- [Machine-readable results and input checksums](results.json)
- [Fix and regression test](https://github.com/psgundla/HEST/commit/bba8421b74c43a490a89a1ecbb3c7fd777406fa0)
- [Original upstream code](https://github.com/mahmoodlab/HEST/commit/3ddb5eaf5bd2a8133e0c0e8015816489a3d99dc3)

## What was checked

The [10x Genomics human colorectal-cancer Visium HD dataset](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc), processed with Space Ranger 3.0.0, supplies 137,051 filtered 16 µm bins, 18,085 genes, and 266,691,885 UMIs. All bins and genes are used, without normalization or cropping. Data are provided by 10x Genomics under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); raw files are not redistributed here.

The script imports the patched function from the checked-out HEST package when run with `--import-hest`. The baseline is the exact function body from the upstream commit, loaded with AST extraction to avoid changing the checkout. An independent integer sparse membership-matrix product supplies expected counts, using HEST's existing spatial-grid contract.

Checks include barcode alignment, unchanged input counts, gene order, every output coordinate, every pooled bin/gene count, per-gene totals, total UMIs, and invariance to chunk size and reversed input order. The raw HDF5 count total is also checked against the [10x metrics CSV](https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_metrics_summary.csv). Its bin-count × mean-UMI product differs from the integer total by about 0.0000015 UMI due to decimal rounding; this check uses less than half a UMI tolerance. Matrix-entry comparisons use zero tolerance.

| Implementation | Chunk rows | Order | Output UMIs | Lost UMIs | Wrong bin/gene entries |
| --- | ---: | --- | ---: | ---: | ---: |
| Original | 1,024 | Source | 213,383,291 | 53,308,594 | 13,675,539 |
| Fixed | 1,024 | Source | 266,691,885 | 0 | 0 |
| Original | 4,096 | Source | 121,371,935 | 145,319,950 | 19,403,587 |
| Fixed | 4,096 | Source | 266,691,885 | 0 | 0 |
| Fixed | 1,024 | Reversed | 266,691,885 | 0 | 0 |

All three fixed runs produce 2,281 nonempty 128 µm bins and the same output-matrix SHA256: `c0f0e5857e20dbf902a6018a6c10f5f41f5a6b4e24359febca044dc6d86226e0`.

## Reproduce

Tested on macOS with Python 3.11.16. Library versions and timings are recorded in `results.json`; the tested TRIDENT dependency resolved to `f3eb7f301ce34f875306b545e6cfefc5d3335a5c`. Follow HEST's main README for system dependencies. These are opt-in checks, not automatic CI downloads.

Clone this evidence branch with full history so the original commit is available:

```sh
git clone --branch validation/visiumhd-143 https://github.com/psgundla/HEST.git HEST-143
cd HEST-143
python3.11 -m venv .venv
.venv/bin/python -m pip install -e .

# Exercise the actual package and all 72 regression combinations.
PYTHONPATH=src:tests .venv/bin/python -m unittest hest_tests.TestVisiumHDPooling -v
```

Download the original archive (15,886,623,172 bytes, approximately 15.9 GB). Only three members are extracted; no WSI image or FASTQ is needed.

```sh
HEST143_DATA="$PWD/data/visiumhd-143"
HEST143_URL="https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer"
mkdir -p "$HEST143_DATA"
curl -fL --retry 2 -o "$HEST143_DATA/binned_outputs.tar.gz" \
  "$HEST143_URL/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz"
curl -fL --retry 2 -o "$HEST143_DATA/metrics_summary.csv" \
  "$HEST143_URL/Visium_HD_Human_Colon_Cancer_metrics_summary.csv"
tar -xzf "$HEST143_DATA/binned_outputs.tar.gz" -C "$HEST143_DATA" \
  binned_outputs/square_016um/filtered_feature_bc_matrix.h5 \
  binned_outputs/square_016um/spatial/tissue_positions.parquet \
  binned_outputs/square_016um/spatial/scalefactors_json.json

PYTHONPATH=src .venv/bin/python validation/visiumhd-143/validate_real_data.py \
  . "$HEST143_DATA/binned_outputs/square_016um" "$HEST143_DATA/rerun.json" \
  --metrics "$HEST143_DATA/metrics_summary.csv" --import-hest
```

Archive SHA256: `a65a7fe9564cba2fb1344c4f586368cfc2cab2d5f8c934657f6994f0eeacb568`. The script records individual file hashes in the output JSON; compare them with the committed results. A successful run exits zero and ends with `PASS: baseline loss reproduced; fixed counts match integer oracle exactly`.

## Limits and tradeoffs

- One real tissue sample was checked, not all 25 slides reported in the issue.
- Normal HEST imports and the pooling function were exercised, not the complete WSI reader, image processing, or broader HEST data suite. Linux CI remains separate.
- The real-data runs use chunks of 1,024 and 4,096 to limit RAM use. The default 50,000-row chunk is covered by the small regression test, not the full real-data matrix.
- This fix retains dense destination storage and densifies each sparse input chunk. It fixes accumulation, not memory scaling or performance. Timings are single-run observations, not a benchmark.
- The reference preserves the current coordinate convention; this is not an independent biological validation of image registration. Source gene-symbol duplicates are retained in their original column order.
- Maximum expected bin/gene count is 62,502, below float32's exact-integer ceiling. This experiment does not establish exactness for counts above that ceiling.
