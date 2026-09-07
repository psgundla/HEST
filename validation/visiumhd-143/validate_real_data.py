"""Validate exact HEST pooling functions on unmodified public 10x count data.

By default this isolates both pooling functions via AST. With --import-hest,
the fixed function is imported from the installed HEST checkout. Neither mode
tests whole-slide image loading or the complete VisiumHDReader pipeline.
See README.md in this directory for data provenance and reproduction commands.
This opt-in check requires the full public dataset; it is not part of normal CI.
"""
from __future__ import annotations

import argparse
import ast
import gc
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import resource
import sys
import subprocess
import time

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

BASELINE = "3ddb5eaf5bd2a8133e0c0e8015816489a3d99dc3"
SOURCE = "src/hest/readers.py"


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load_pool(source, filename):
    node = next(node for node in ast.parse(source).body
                if isinstance(node, ast.FunctionDef) and node.name == "pool_bins_visiumhd")
    future = ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0)
    module = ast.fix_missing_locations(ast.Module(body=[future, node], type_ignores=[]))
    namespace = {"np": np, "pd": pd, "math": math}
    exec(compile(module, filename, "exec"), namespace)
    return namespace[node.name]


def matrix_digest(matrix):
    digest = hashlib.sha256()
    for array in (matrix.data, matrix.indices, matrix.indptr):
        digest.update(memoryview(np.ascontiguousarray(array)))
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo", type=Path)
    parser.add_argument("data", type=Path, help="Unmodified square_016um directory")
    parser.add_argument("report", type=Path)
    parser.add_argument("--metrics", type=Path, help="Original 10x metrics_summary.csv")
    parser.add_argument("--import-hest", action="store_true", help="Test the installed HEST function instead of AST-isolating the patch")
    args = parser.parse_args()
    files = [args.data / "filtered_feature_bc_matrix.h5",
             args.data / "spatial/tissue_positions.parquet",
             args.data / "spatial/scalefactors_json.json"]
    fixed_source = (args.repo / SOURCE).read_text()
    original_source = subprocess.check_output(
        ["git", "-C", str(args.repo), "show", f"{BASELINE}:{SOURCE}"], text=True)
    pools = {"original": load_pool(original_source, f"{BASELINE}:{SOURCE}"),
             "fixed": load_pool(fixed_source, str(args.repo / SOURCE))}
    if args.import_hest:
        import inspect
        from hest.readers import pool_bins_visiumhd
        assert Path(inspect.getfile(pool_bins_visiumhd)).resolve() == (args.repo / SOURCE).resolve()
        pools["fixed"] = pool_bins_visiumhd
    data = sc.read_10x_h5(files[0])
    positions = pd.read_parquet(files[1]).set_index("barcode", verify_integrity=True)
    assert data.obs_names.is_unique
    assert data.obs_names.isin(positions.index).all(), "Missing barcode coordinates"
    reader_aligned = pd.merge(data.obs, positions, how="inner", left_index=True, right_index=True)
    assert reader_aligned.index.equals(data.obs_names), "HEST reader merge changes barcode order"
    pd.testing.assert_frame_equal(positions.loc[data.obs_names], reader_aligned, check_names=False)
    data.obs = reader_aligned
    assert data.obs["in_tissue"].eq(1).all()
    assert np.isfinite(data.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]]).all().all()
    pixel_size = float(json.loads(files[2].read_text())["microns_per_pixel"])
    assert np.isfinite(pixel_size) and pixel_size > 0
    assert sparse.isspmatrix_csr(data.X)
    assert np.all(data.X.data >= 0) and np.all(data.X.data == np.floor(data.X.data))
    input_digest = matrix_digest(data.X)
    input_total = int(data.X.sum(dtype=np.int64))
    with h5py.File(files[0], "r") as handle:
        assert set(handle["matrix/features/feature_type"].asstr()[:]) == {"Gene Expression"}
        assert input_total == int(handle["matrix/data"][:].sum(dtype=np.int64))
    input_genes = np.asarray(data.X.sum(axis=0, dtype=np.int64)).ravel()
    vendor = None
    if args.metrics:
        metrics = pd.read_csv(args.metrics).iloc[0]
        vendor_bins = int(metrics["Number of Bins Under Tissue 16 µm"])
        vendor_total = vendor_bins * float(metrics["Mean UMIs Under Tissue per Bin 16 µm"])
        assert vendor_bins == data.n_obs
        # CSV means are rounded floating-point values; compare at whole-UMI resolution.
        assert abs(vendor_total - input_total) < 0.5
        vendor = {"file": args.metrics.name, "sha256": sha256(args.metrics),
                  "bins": vendor_bins, "implied_total_umis": vendor_total,
                  "rounding_difference_umis": vendor_total - input_total}

    # Independent integer sparse membership product. Keep HEST's existing spatial
    # grid contract; do not use its accumulation code as the expected answer.
    xy = data.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    origin = xy.min(axis=0)
    source_pixels, target_pixels = 16 / pixel_size, 128 / pixel_size
    grid_shape = np.ceil((xy.max(axis=0) - origin + source_pixels) / target_pixels).astype(int)
    grid_xy = np.floor((xy - origin + source_pixels // 2) / target_pixels).astype(int)
    destinations = grid_xy[:, 1] * grid_shape[0] + grid_xy[:, 0]
    membership = sparse.csr_matrix(
        (np.ones(data.n_obs, dtype=np.int64), (destinations, np.arange(data.n_obs))),
        shape=(int(np.prod(grid_shape)), data.n_obs))
    oracle = membership @ data.X.astype(np.int64)
    kept = np.flatnonzero(np.asarray(oracle.sum(axis=1)).ravel())
    oracle = oracle[kept].tocsr()
    expected_xy = origin + np.column_stack((kept % grid_shape[0], kept // grid_shape[0])) * target_pixels + target_pixels // 2
    assert int(oracle.sum()) == input_total
    assert oracle.data.max() < 2**24, "Counts exceed exact float32 integer range"

    report = {
        "scope": ("Full filtered 16um matrix, all genes; fixed function imported through HEST, original AST-isolated; no WSI pipeline"
                  if args.import_hest else "Full filtered 16um matrix, all genes; exact-function isolation, not full HEST integration"),
        "default_chunk_limit": "50000-row real-data run omitted for RAM safety; covered by synthetic regressions",
        "baseline_commit": BASELINE,
        "python_version": sys.version.split()[0],
        "platform": sys.platform,
        "source_sha256": {"original": hashlib.sha256(original_source.encode()).hexdigest(),
                          "fixed": hashlib.sha256(fixed_source.encode()).hexdigest()},
        "files": {path.relative_to(args.data).as_posix(): {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in files},
        "versions": {name: importlib.metadata.version(name) for name in
                     ("numpy", "pandas", "scipy", "anndata", "scanpy", "h5py", "pyarrow")},
        "input": {"bins": data.n_obs, "genes": data.n_vars, "nonzero_entries": data.X.nnz,
                  "umis": input_total, "microns_per_pixel": pixel_size,
                  "all_barcodes_matched": True, "matrix_sha256": input_digest},
        "oracle": {"nonempty_128um_bins": len(kept), "umis": int(oracle.sum()),
                   "max_bin_gene_count": int(oracle.data.max())},
        "vendor_metrics": vendor,
        "runs": [],
    }
    print(json.dumps({"input": report["input"], "oracle": report["oracle"]}), flush=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    fixed_digest = None
    for version, chunk, order in [
        ("original", 1024, "source"), ("fixed", 1024, "source"),
        ("original", 4096, "source"), ("fixed", 4096, "source"),
        ("fixed", 1024, "reversed"),
    ]:
        case = data if order == "source" else data[::-1].copy()
        print(f"Running {version}, chunk={chunk}, order={order}", flush=True)
        started = time.perf_counter()
        result = pools[version](case, pixel_size=pixel_size, chunk_len=chunk)
        elapsed = time.perf_counter() - started
        total = int(result.X.sum(dtype=np.float64))
        assert result.var_names.equals(data.var_names)
        # Align to oracle by physical coordinates, not HEST's compacted labels.
        output_keys = pd.MultiIndex.from_arrays(result.obsm["spatial"].T)
        oracle_keys = pd.MultiIndex.from_arrays(expected_xy.T)
        assert output_keys.is_unique
        locations = oracle_keys.get_indexer(output_keys)
        assert np.all(locations >= 0)
        mismatched = 0
        max_error = 0
        for start in range(0, result.n_obs, 128):
            actual = result.X[start:start + 128]
            expected = oracle[locations[start:start + 128]].toarray()
            mismatched += int(np.count_nonzero(actual != expected))
            max_error = max(max_error, int(np.abs(actual - expected).max()))
        missing = np.setdiff1d(np.arange(len(kept)), locations)
        mismatched += int(oracle[missing].nnz)
        if len(missing):
            max_error = max(max_error, int(oracle[missing].max()))
        changed_genes = int(np.count_nonzero(result.X.sum(axis=0, dtype=np.float64) != input_genes))
        digest = hashlib.sha256(memoryview(np.ascontiguousarray(result.X))).hexdigest()
        row = {"version": version, "chunk_len": chunk, "row_order": order,
               "output_bins": result.n_obs, "output_umis": total,
               "lost_umis": input_total - total, "loss_percent": 100 * (input_total - total) / input_total,
               "missing_output_bins": len(missing), "mismatched_bin_gene_entries": mismatched,
               "max_absolute_count_error": max_error, "genes_with_changed_totals": changed_genes,
               "seconds": round(elapsed, 3), "output_matrix_sha256": digest,
               "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)}
        report["runs"].append(row)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(row), flush=True)
        if version == "fixed":
            assert total == input_total and mismatched == 0 and changed_genes == 0 and len(missing) == 0
            np.testing.assert_array_equal(result.obsm["spatial"], expected_xy)
            if fixed_digest is None:
                fixed_digest = digest
            assert digest == fixed_digest, "Chunk size or row order changed output"
        else:
            assert total < input_total and mismatched > 0, "Expected baseline bug not reproduced"
        assert matrix_digest(data.X) == input_digest, "Original input counts mutated"
        del result, case
        gc.collect()
    report["status"] = "PASS: baseline loss reproduced; fixed counts match integer oracle exactly"
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(report["status"], flush=True)


if __name__ == "__main__":
    main()
