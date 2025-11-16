#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Themis profiling pipeline (paired/single).

双端：
  ganon classify/report → 物种/菌株表 → 生成 color_mapping.in →
  ggcat build → fastp 合并 R1/R2 → ggcat query →
  阈值过滤 + length-corrected abundance →
  与 ganon 结果自适应混合（需 predict_spy.ID.abundance 为两列：ID\tabundance）→ 输出 species_abundance.txt

单端：
  ganon classify/report → 物种/菌株表 → 生成 color_mapping.in →
  ggcat build → ggcat query（单端）→
  阈值过滤 + length-corrected abundance（不与 ganon 混合）→ 输出 species_abundance.txt
"""

import subprocess
from pathlib import Path
import math
import shutil
import argparse
import importlib

from .utils import ensure_dir
from . import ganon_wrapper, ggcat_wrapper

# 延迟导入 themis_scripts 子模块（避免构建/测试阶段提前触发 pandas 等重依赖）
def _ts(mod: str):
    return importlib.import_module(f"themis_scripts.{mod}")

# =========================
# 对外主入口
# =========================

def run(
    reads,
    single,
    db_prefix,
    out_prefix,
    report_type="abundance",
    file_info=None,
    threads=8,
    k=31,
):
    """
    Themis profile 主入口

    参数
    ----
    reads:
      - paired: [R1, R2]
      - single: [R]
    single: 是否单端（True=单端；False=双端）
    db_prefix: ganon 数据库前缀
    out_prefix: 输出目录前缀，例如 results/sample1
    file_info: Ref 信息表（RefDB_xxx_genomes_info.txt），用于物种-菌株→FASTA 路径映射
    threads: 线程数（传递给 ganon/ggcat/fastp）
    k: ggcat 的 k-mer 长度
    """
    out_prefix = Path(out_prefix).absolute()
    ensure_dir(out_prefix)

    if single:
        return run_single(
            reads=reads,
            db_prefix=db_prefix,
            out_prefix=out_prefix,
            report_type=report_type,
            ref_info=file_info,
            threads=threads,
            k=k,
        )
    else:
        return run_paired(
            reads=reads,
            db_prefix=db_prefix,
            out_prefix=out_prefix,
            report_type=report_type,
            ref_info=file_info,
            threads=threads,
            k=k,
        )

# =========================
# 双端流程
# =========================

def run_paired(reads, db_prefix, out_prefix,
               report_type, ref_info, threads, k):
    if len(reads) != 2:
        raise SystemExit("[Themis] Paired-end mode requires two -r reads (R1 and R2).")
    if not ref_info:
        raise SystemExit("[Themis] --ref-info is required.")

    r1, r2 = map(str, reads)

    work_ganon = out_prefix / "query_r1"
    work_db = out_prefix / "database_filter_ccDBG"
    work_q = out_prefix / "res_themis"

    ensure_dir(work_ganon)
    ensure_dir(work_db)
    ensure_dir(work_q)

    # 1) ganon classify (paired)
    ganon_out_prefix = work_ganon / "results"
    ganon_wrapper.run_classify_paired(
        db_prefix=db_prefix,
        read1=r1,
        read2=r2,
        out_prefix=str(ganon_out_prefix),
        threads=threads,
        report_type=report_type,
    )

    # 2) ganon report + species_abundance + strain_abundance (+ predict_spy.ID.abundance 两列)
    ganon_tre = work_ganon / "tax_profile.tre"
    ganon_species = work_ganon / "species_abundance.txt"
    ganon_strain = work_ganon / "strain_abundance.txt"
    ganon_predict_spy = work_ganon / "predict_spy.ID.abundance"

    ganon_wrapper.run_report_and_postprocess(
        db_prefix=db_prefix,
        classify_prefix=str(ganon_out_prefix),
        tre_out=str(ganon_tre),
        species_out=str(ganon_species),
        strain_out=str(ganon_strain),
        predict_spy_out=str(ganon_predict_spy),  # 必须是两列：ID\tabundance
    )

    # 3) 基于 species 行数选择 topk 或普通筛选，生成 tmp_id_table.tsv
    line_count = _count_non_empty_lines(str(ganon_species))
    tmp_id_table = work_q / "tmp_id_table.tsv"

    if line_count > 1000:
        top_tsv = work_db / "ganon_species_strain_top.tsv"
        _ts("make_ganon_pred_symlinks_topk_singleton_filter").run(
            ref_info=str(ref_info),
            ganon_species=str(ganon_species),
            ganon_strain=str(ganon_strain),
            top_k=3,
            singleton_min_abund=7e-7,
            out_tsv=str(top_tsv),
        )
        shutil.copyfile(top_tsv, tmp_id_table)
    elif line_count >= 1:
        sel_tsv = work_db / "ganon_species_strain_selected.tsv"
        _ts("make_ganon_pred_symlinks").run(
            ref_info=str(ref_info),
            ganon_species=str(ganon_species),
            ganon_strain=str(ganon_strain),
            out_tsv=str(sel_tsv),
        )
        shutil.copyfile(sel_tsv, tmp_id_table)
    else:
        raise SystemExit("[Themis] Too few ganon species hits in paired-end mode.")

    # 4) 生成 color_mapping.in：new_strain_taxid \t fasta_path
    color_map = work_db / "color_mapping.in"
    _make_color_mapping_from_mapping_tsv(str(tmp_id_table), str(color_map))

    # 5) ggcat build
    ggcat_db = work_db / "ggcatDB.fasta.lz4"
    ggcat_wrapper.run_build(
        k=k,
        threads=threads,
        color_mapping=str(color_map),
        output=str(ggcat_db),
    )

    # 6) fastp 合并双端 reads → ggcat 用的单 fastq
    ggcat_reads = work_q / "reads_for_ggcat.fastq"
    _prepare_ggcat_reads_paired(r1, r2, str(ggcat_reads), threads=threads)

    # 7) ggcat query
    ggcat_prefix = work_q / "query_ggcatDB"
    ggcat_wrapper.run_query(
        db=str(ggcat_db),
        reads=[str(ggcat_reads)],
        k=k,
        threads=threads,
        out_prefix=str(ggcat_prefix),
        single=False,
    )

    # 8) 阈值 + length-corrected abundance + 与 ganon 自适应混合
    final_abundance = out_prefix / "species_abundance.txt"
    _run_threshold_and_mix_for_paired(
        ggcat_prefix=str(ggcat_prefix),
        tmp_id_table=str(tmp_id_table),
        ganon_predict_spy=str(ganon_predict_spy),
        tre_file=str(ganon_tre),
        reads_path=str(ggcat_reads),
        output=str(final_abundance),
    )

    print(f"[Themis] Paired-end profiling done → {final_abundance}")
    return str(final_abundance)

# =========================
# 单端流程（不与 ganon 混合）
# =========================

def run_single(reads, db_prefix, out_prefix,
               report_type, ref_info, threads, k):
    if len(reads) != 1:
        raise SystemExit("[Themis] --single mode requires exactly one -r reads file.")
    if not ref_info:
        raise SystemExit("[Themis] --ref-info is required.")

    read = str(reads[0])

    work_ganon = out_prefix / "query_r1"
    work_db = out_prefix / "database_filter_ccDBG"
    work_q = out_prefix / "res_themis"

    ensure_dir(work_ganon)
    ensure_dir(work_db)
    ensure_dir(work_q)

    # 1) ganon classify (single)
    ganon_out_prefix = work_ganon / "results"
    ganon_wrapper.run_classify_single(
        db_prefix=db_prefix,
        read=read,
        out_prefix=str(ganon_out_prefix),
        threads=threads,
        report_type=report_type,
    )

    # 2) ganon report + species/strain abundance
    ganon_tre = work_ganon / "tax_profile.tre"
    ganon_species = work_ganon / "species_abundance.txt"
    ganon_strain = work_ganon / "strain_abundance.txt"

    ganon_wrapper.run_report_and_postprocess(
        db_prefix=db_prefix,
        classify_prefix=str(ganon_out_prefix),
        tre_out=str(ganon_tre),
        species_out=str(ganon_species),
        strain_out=str(ganon_strain),
    )

    # 3) species/strain + ref_info → tmp_id_table
    tmp_tsv = work_db / "ganon_species_strain_selected.tsv"
    tmp_id_table = work_q / "tmp_id_table.tsv"

    _ts("make_ganon_pred_symlinks").run(
        ref_info=str(ref_info),
        ganon_species=str(ganon_species),
        ganon_strain=str(ganon_strain),
        out_tsv=str(tmp_tsv),
    )
    shutil.copyfile(tmp_tsv, tmp_id_table)

    # 4) color_mapping.in
    color_map = work_db / "color_mapping.in"
    _make_color_mapping_from_mapping_tsv(str(tmp_id_table), str(color_map))

    # 5) ggcat build
    ggcat_db = work_db / "ggcatDB.fasta.lz4"
    ggcat_wrapper.run_build(
        k=k,
        threads=threads,
        color_mapping=str(color_map),
        output=str(ggcat_db),
    )

    # 6) ggcat query（单端直接用原始 read）
    ggcat_prefix = work_q / "query_ggcatDB"
    ggcat_wrapper.run_query(
        db=str(ggcat_db),
        reads=[read],
        k=k,
        threads=threads,
        out_prefix=str(ggcat_prefix),
        single=True,
    )

    # 7) 阈值 + length-corrected abundance（不与 ganon 混合）
    final_abundance = out_prefix / "species_abundance.txt"
    _run_threshold_and_lc_for_single(
        ggcat_prefix=str(ggcat_prefix),
        tmp_id_table=str(tmp_id_table),
        reads_path=read,
        output=str(final_abundance),
    )

    print(f"[Themis] Single-end profiling done → {final_abundance}")
    return str(final_abundance)

# =========================
# 工具函数 & 阈值/混合
# =========================

def _count_non_empty_lines(path):
    c = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                c += 1
    return c

def _prepare_ggcat_reads_paired(r1, r2, out_fastq, threads=8):
    """
    使用 fastp 将双端 reads 合并为一个 fastq，供 ggcat 使用。
    fastp 通过 conda 依赖提供，需在 PATH 中。
    """
    cmd = [
        "fastp",
        "-i", r1,
        "-I", r2,
        "--thread", str(threads),
        "--stdout",
    ]
    with open(out_fastq, "w") as out_f:
        try:
            subprocess.run(cmd, check=True, stdout=out_f)
        except FileNotFoundError:
            raise SystemExit(
                "[Themis][error] fastp not found. "
                "Please install via conda (fastp) or use the Themis conda package."
            )
        except subprocess.CalledProcessError as e:
            raise SystemExit(
                f"[Themis][error] fastp failed with exit code {e.returncode}"
            )

def _make_color_mapping_from_mapping_tsv(mapping_tsv, out_path):
    """
    输入:  TSV（含列：species_taxid, strain_taxid, new_strain_taxid, id）
    输出:  color_mapping.in （两列：new_strain_taxid\tfasta_path）
    """
    with open(mapping_tsv, "r", encoding="utf-8") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        header = next(fin, None)
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            new_strain_taxid = parts[2]
            fasta_path = parts[3]
            if new_strain_taxid and fasta_path:
                fout.write(f"{new_strain_taxid}\t{fasta_path}\n")

def _run_threshold_and_mix_for_paired(
    ggcat_prefix,
    tmp_id_table,
    ganon_predict_spy,
    tre_file,
    reads_path,
    output,
):
    """
    双端：threshold + length-corrected DBG abundance + DBG/Ganon 自适应混合。
    需要 ganon_predict_spy 为两列：ID \t abundance
    """
    species_counts = f"{ggcat_prefix}.species_counts.tsv"

    # 1) remains & species
    remains = 0.0
    species = 0
    with open(species_counts, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            remains += v
            species += 1

    # 2) total_reads
    with open(reads_path, "r", encoding="utf-8") as rf:
        total_lines = sum(1 for _ in rf)
    total_reads = max(total_lines // 4, 1)

    if species == 0:
        raise SystemExit("[Themis] No species in ggcat species_counts.tsv (paired).")

    ratio = remains / total_reads
    expr = (remains / 1_000_000.0) * (ratio ** 2) * (remains / species)
    threshold = math.sqrt(expr) if expr > 0 else 0.0

    # 3) 过滤 > threshold
    filtered_counts = f"{ggcat_prefix}.species_counts_more{threshold:.6g}.tsv"
    with open(species_counts, "r", encoding="utf-8", errors="ignore") as fin, \
         open(filtered_counts, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            if v > threshold:
                fout.write(line)

    # 4) length_corrected_abundance（DBG）
    abund_more = f"{ggcat_prefix}.species_abundance_more{threshold:.6g}.tsv"
    _ts("length_corrected_abundance").run(
        counts_file=filtered_counts,
        mapping_file=tmp_id_table,
        output_file=abund_more,
    )

    # 5) DBG 两列：id \t abundance
    dbg_tmp = f"{abund_more}.tmp"
    with open(abund_more, "r", encoding="utf-8", errors="ignore") as fin, \
         open(dbg_tmp, "w", encoding="utf-8") as fout:
        header = fin.readline()  # 跳过表头
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            fout.write(f"{parts[0]}\t{parts[5]}\n")

    # 6) 初次混合（w=0）：看 Ganon 结构
    raw_mix = f"{ggcat_prefix}.raw_filtered_abundance_prediction.tsv"
    _ts("mix_predictions").run(
        dbg_file=dbg_tmp,
        ganon_file=ganon_predict_spy,
        weight=0.0,
        output=raw_mix,
    )

    # 7) 计算 top-2 差值
    def top2_diff(path):
        vals = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                try:
                    v = float(parts[1])
                except ValueError:
                    continue
                vals.append(v)
        if not vals:
            return 0.0
        vals.sort(reverse=True)
        max1 = vals[0]
        max2 = vals[1] if len(vals) > 1 else 0.0
        return max1 - max2

    diffDBG = top2_diff(dbg_tmp)
    diffGANON = top2_diff(raw_mix)
    i = diffDBG - diffGANON

    # 8) TRE 权重（unclassified 百分比 / 5；裁剪到 [0,1]）
    def tre_weight_calc(tre_file):
        w = 0.0
        found = False
        try:
            with open(tre_file, "r", encoding="utf-8") as f:
                for line in f:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if parts and parts[0] == "unclassified" and len(parts) >= 9:
                        try:
                            val = float(parts[8])
                        except ValueError:
                            val = 0.0
                        w = val / 100.0 / 5.0
                        found = True
                        break
        except FileNotFoundError:
            return 0.0
        if not found:
            w = 0.0
        return max(0.0, min(1.0, w))

    tre_w = tre_weight_calc(tre_file)

    # 9) 最终权重规则
    if i > 0.2:
        weight = tre_w
    elif 0.02 < i <= 0.2:
        weight = 1.0
    elif -0.2 < i <= 0.02:
        weight = tre_w
    else:
        weight = 1.0

    # 10) 最终混合
    final_tmp = f"{ggcat_prefix}.mix_abundance_prediction.tsv"
    _ts("mix_predictions").run(
        dbg_file=dbg_tmp,
        ganon_file=ganon_predict_spy,
        weight=weight,
        output=final_tmp,
    )

    with open(final_tmp, 'r', encoding='utf-8') as fin, \
        open(output, 'w', encoding='utf-8') as fout:
        # 写入表头
        fout.write("speciesID\tabundance\n")
        # 复制内容
        shutil.copyfileobj(fin, fout)

def _run_threshold_and_lc_for_single(
    ggcat_prefix,
    tmp_id_table,
    reads_path,
    output,
):
    """
    单端流程：
      1) 从 ggcat 的 species_counts.tsv 计算阈值并过滤；
      2) 用映射表（tmp_id_table：含 species_taxid 与多个 strain fasta 路径）
         调用 length_corrected_abundance.run() 进行长度校正并归一化；
      3) 输出最终的 species_abundance.txt。
      注：单端不进行 DBG/Ganon 混合。
    """
    species_counts = f"{ggcat_prefix}.species_counts.tsv"

    remains = 0.0
    species = 0
    with open(species_counts, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            remains += v
            species += 1

    with open(reads_path, "r", encoding="utf-8") as rf:
        total_lines = sum(1 for _ in rf)
    total_reads = max(total_lines // 4, 1)

    if species == 0:
        raise SystemExit("[Themis] No species in ggcat species_counts.tsv (single-end).")

    ratio = remains / total_reads
    expr = (remains / 1_000_000.0) * (ratio ** 2) * (remains / species)
    threshold = math.sqrt(expr) if expr > 0 else 0.0

    filtered_counts = f"{species_counts}.more{threshold:.6g}.tsv"
    with open(species_counts, "r", encoding="utf-8") as fin, \
         open(filtered_counts, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            if v > threshold:
                fout.write(line)

    abund_more = f"{species_counts}.abundance_more{threshold:.6g}.tsv"
    _ts("length_corrected_abundance").run(
        counts_file=filtered_counts,
        mapping_file=tmp_id_table,
        output_file=abund_more,
    )

    pairs = []
    with open(abund_more, "r", encoding="utf-8", errors="ignore") as fin:
        _ = fin.readline()  # 跳过表头
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                sid = parts[0].strip()
                try:
                    rel = float(parts[5])  # rel_abundance
                except ValueError:
                    continue
                if sid:
                    pairs.append((sid, rel))

    pairs.sort(key=lambda x: x[1], reverse=True)  # 如不想排序可删除

    with open(output, "w", encoding="utf-8") as fout:
        fout.write("speciesID\tabundance\n")
        for sid, rel in pairs:
            fout.write(f"{sid}\t{rel}\n")

# =========================
# CLI 封装
# =========================

def cli():
    p = argparse.ArgumentParser("themis-profile")
    p.add_argument(
        "-r", "--reads", action="append", required=True,
        help="Input reads. Paired-end: -r R1 -r R2. Single-end: -r R.",
    )
    p.add_argument(
        "--single", action="store_true",
        help="Single-end mode.",
    )
    p.add_argument(
        "--db-prefix", required=True,
        help="Ganon DB prefix built by themis-build-custom.",
    )
    p.add_argument(
        "--ref-info", required=True,
        help="Reference info TSV (RefDB_xxx_genomes_info.txt).",
    )
    p.add_argument(
        "--out", required=True,
        help="Output directory.",
    )
    p.add_argument(
        "--threads", type=int, default=8,
        help="Number of threads.",
    )
    p.add_argument(
        "-k", "--kmer", type=int, default=31,
        help="k-mer size for ggcat.",
    )
    args = p.parse_args()

    run(
        reads=args.reads,
        single=args.single,
        db_prefix=args.db_prefix,
        out_prefix=args.out,
        report_type="abundance",
        file_info=args.ref_info,
        threads=args.threads,
        k=args.kmer,
    )
