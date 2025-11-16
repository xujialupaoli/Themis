use crate::structs::query_colored_counters::{ColorsRange, QueryColoredCountersSerializer};
use crate::ColoredQueryOutputFormat;
use colors::colors_manager::{ColorMapReader, ColorsManager, ColorsMergeManager};
use config::{
    get_compression_level_info, get_memory_mode, ColorIndexType, SwapPriority,
    DEFAULT_PREFETCH_AMOUNT, KEEP_FILES, QUERIES_COUNT_MIN_BATCH,
};
use flate2::Compression;
use ggcat_logging::UnrecoverableErrorLogging;
use hashes::HashFunctionFactory;
use nightly_quirks::prelude::*;
use parallel_processor::buckets::readers::compressed_binary_reader::CompressedBinaryReader;
use parallel_processor::buckets::readers::BucketReader;
use parallel_processor::buckets::writers::compressed_binary_writer::CompressedBinaryWriter;
use parallel_processor::buckets::{LockFreeBucket, SingleBucket};
use parallel_processor::memory_fs::RemoveFileMode;
use parallel_processor::phase_times_monitor::PHASES_TIMES_MONITOR;
use parking_lot::{Condvar, Mutex};
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::DerefMut;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

// =============================
// Optional JSONL output writer
// =============================

enum QueryOutputFileWriter {
    Plain(File),
    LZ4Compressed(lz4::Encoder<File>),
    GzipCompressed(flate2::write::GzEncoder<File>),
}

impl Write for QueryOutputFileWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            QueryOutputFileWriter::Plain(w) => w.write(buf),
            QueryOutputFileWriter::LZ4Compressed(w) => w.write(buf),
            QueryOutputFileWriter::GzipCompressed(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            QueryOutputFileWriter::Plain(w) => w.flush(),
            QueryOutputFileWriter::LZ4Compressed(w) => w.flush(),
            QueryOutputFileWriter::GzipCompressed(w) => w.flush(),
        }
    }
}

// =============================
// Species aggregation helpers
// =============================

#[derive(Clone)]
enum Decision {
    Species(String),
    Discard(&'static str, String),
}

#[inline]
fn species_of(name: &str) -> &str {
    name.splitn(2, '_').next().unwrap_or(name)
}

#[inline]
fn winners_of_max(matches: &[(String, f64)]) -> (f64, Vec<&str>) {
    if matches.is_empty() {
        return (0.0, vec![]);
    }
    let mut best = 0.0f64;
    for &(_, v) in matches {
        if v > best {
            best = v;
        }
    }
    let wins = matches
        .iter()
        .filter(|(_, v)| *v == best)
        .map(|(k, _)| k.as_str())
        .collect();
    (best, wins)
}

#[inline]
fn accept_single(_idx: usize, matches: &[(String, f64)], t: f64) -> Decision {
    let (rmax, winners) = winners_of_max(matches);

    if rmax <= t {
        return Decision::Discard(
            "below_threshold",
            format!("Rmax={:.6} T={:.6}", rmax, t),
        );
    }

    if winners.is_empty() {
        return Decision::Discard("no_match", String::new());
    }

    let mut sp_set: Vec<&str> = winners.iter().map(|w| species_of(w)).collect();
    sp_set.sort_unstable();
    sp_set.dedup();

    if sp_set.len() > 1 {
        return Decision::Discard(
            "multi_species",
            format!("winners_species={}", sp_set.join(",")),
        );
    }

    Decision::Species(sp_set[0].to_string())
}

// ---------- 阈值采样逻辑 ----------

// FIRST_N 是“最多采多少个 Rmax 样本用来估计阈值”，
// 如果总 query 数 < FIRST_N，就“有多少用多少”。
const FIRST_N: usize = 50_000;
const PCT: f64 = 0.30; // p30
const THR_MIN: f64 = 1.0 / 3.0;
const THR_MAX: f64 = 2.0 / 3.0;

struct RmaxSampler {
    vals: Vec<f64>,
}

impl RmaxSampler {
    fn new() -> Self {
        Self {
            vals: Vec::with_capacity(FIRST_N),
        }
    }

    #[inline]
    fn try_push(&mut self, r: f64) {
        if self.vals.len() < FIRST_N {
            self.vals.push(r);
        }
    }

    /// 根据当前采样计算阈值:
    /// - 有样本: p30 裁剪到 [THR_MIN, THR_MAX]
    /// - 无样本: 退化为 THR_MIN
    fn compute_threshold(&mut self) -> f64 {
        if self.vals.is_empty() {
            return THR_MIN;
        }

        self.vals
            .sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let n = self.vals.len();
        let pos = PCT * (n as f64 - 1.0);
        let i = pos.floor() as usize;
        let j = pos.ceil() as usize;

        let p = if i == j {
            self.vals[i]
        } else {
            let frac = pos - (i as f64);
            self.vals[i] * (1.0 - frac) + self.vals[j] * frac
        };

        p.clamp(THR_MIN, THR_MAX)
    }
}

// ---------- 全局聚合状态 ----------

struct GlobalAgg {
    is_paired: bool,

    // species -> (reads_count, indices)
    species: HashMap<String, (usize, Vec<usize>)>,

    // idx -> (reason, detail)
    discard: HashMap<usize, (&'static str, String)>,

    // reason -> count
    reason_stats: HashMap<&'static str, usize>,

    // 仅 paired 模式用：等待 mate 的决策
    pair_buffer: HashMap<usize, Decision>,

    // 所有 query 的 (idx, matches)，等算出 T 后再统一决策
    pending: Vec<(usize, Vec<(String, f64)>)>,

    sampler: RmaxSampler,
}

impl GlobalAgg {
    fn new(is_paired: bool) -> Self {
        Self {
            is_paired,
            species: HashMap::new(),
            discard: HashMap::new(),
            reason_stats: HashMap::new(),
            pair_buffer: HashMap::new(),
            pending: Vec::new(),
            sampler: RmaxSampler::new(),
        }
    }

    #[inline]
    fn inc_reason(&mut self, r: &'static str, n: usize) {
        *self.reason_stats.entry(r).or_insert(0) += n;
    }

    // single / HiFi: 每条 read 独立计数
    fn apply_decision_single(&mut self, idx: usize, dec: Decision) {
        match dec {
            Decision::Species(sp) => {
                let entry = self.species.entry(sp).or_insert((0, Vec::new()));
                entry.0 += 1;
                entry.1.push(idx);
            }
            Decision::Discard(r, d) => {
                self.discard.insert(idx, (r, d));
                self.inc_reason(r, 1);
            }
        }
    }

    // paired: 要求两个 mate 同物种才 +2
    fn merge_pair(&mut self, idx: usize, dec: Decision) {
        let mate = if idx % 2 == 0 { idx + 1 } else { idx - 1 };

        if let Some(other) = self.pair_buffer.remove(&mate) {
            match (dec, other) {
                (Decision::Species(a), Decision::Species(b)) => {
                    if a == b {
                        let entry = self.species.entry(a).or_insert((0, Vec::new()));
                        entry.0 += 2;

                        let (i1, i2) = if idx % 2 == 0 { (idx, mate) } else { (mate, idx) };
                        entry.1.push(i1);
                        entry.1.push(i2);
                    } else {
                        let msg = format!("{} vs {}", a, b);
                        self.discard.insert(idx, ("pair_conflict", msg.clone()));
                        self.discard.insert(mate, ("pair_conflict", msg));
                        self.inc_reason("pair_conflict", 2);
                    }
                }
                (Decision::Species(_), Decision::Discard(r, d))
                | (Decision::Discard(r, d), Decision::Species(_)) => {
                    self.discard
                        .insert(idx, ("pair_missing", format!("mate={}", mate)));
                    self.discard.insert(mate, (r, d));
                    self.inc_reason("pair_missing", 1);
                    self.inc_reason(r, 1);
                }
                (Decision::Discard(r1, d1), Decision::Discard(r2, d2)) => {
                    self.discard.insert(idx, (r1, d1));
                    self.discard.insert(mate, (r2, d2));
                    self.inc_reason(r1, 1);
                    self.inc_reason(r2, 1);
                }
            }
        } else {
            self.pair_buffer.insert(idx, dec);
        }
    }

    fn apply_decision(&mut self, idx: usize, dec: Decision) {
        if self.is_paired {
            self.merge_pair(idx, dec);
        } else {
            self.apply_decision_single(idx, dec);
        }
    }

    // paired: 标记孤立 mate；single: 清空即可
    fn finalize_orphans(&mut self) {
        if !self.is_paired {
            self.pair_buffer.clear();
            return;
        }

        let keys: Vec<usize> = self.pair_buffer.keys().cloned().collect();
        for i in keys {
            self.discard
                .insert(i, ("pair_missing", "file_end_or_missing_mate".to_string()));
            self.inc_reason("pair_missing", 1);
        }
        self.pair_buffer.clear();
    }
}

// =============================
// 主入口
// =============================

pub fn colored_query_output<MH: HashFunctionFactory, CX: ColorsManager>(
    colormap: &<CX::ColorsMergeManagerType as ColorsMergeManager>::GlobalColorsTableReader,
    mut colored_query_buckets: Vec<SingleBucket>,
    output_file: PathBuf,
    temp_dir: PathBuf,
    query_kmers_count: &[u64],
    colored_query_output_format: ColoredQueryOutputFormat,
    emit_jsonl: bool,
    is_paired: bool,
) -> anyhow::Result<()> {
    PHASES_TIMES_MONITOR
        .write()
        .start_phase("phase: colored query output".to_string());

    let buckets_count = colored_query_buckets.len();

    let max_bucket_queries_count = (((query_kmers_count.len() + 1) as u64)
        .nq_div_ceil(QUERIES_COUNT_MIN_BATCH)
        * QUERIES_COUNT_MIN_BATCH) as usize;

    static OPS_COUNT: AtomicUsize = AtomicUsize::new(0);
    static COL_COUNT: AtomicUsize = AtomicUsize::new(0);

    colored_query_buckets.reverse();
    let buckets_channel = Mutex::new(colored_query_buckets);

    // ===== 可选 JSONL 输出 =====
    let (maybe_query_output, output_sync_condvar, output_path_final) = if emit_jsonl {
        let output_file = if output_file.extension().is_none() {
            output_file.with_extension("jsonl")
        } else {
            output_file.clone()
        };

        let query_output_file = File::create(&output_file).log_unrecoverable_error_with_data(
            "Cannot create output file",
            output_file.display(),
        )?;

        let query_output = Mutex::new((
            BufWriter::new(match output_file.extension().and_then(|e| e.to_str()) {
                Some("lz4") => QueryOutputFileWriter::LZ4Compressed(
                    lz4::EncoderBuilder::new()
                        .level(4)
                        .build(query_output_file)
                        .unwrap(),
                ),
                Some("gz") => QueryOutputFileWriter::GzipCompressed(
                    flate2::GzBuilder::new().write(query_output_file, Compression::default()),
                ),
                _ => QueryOutputFileWriter::Plain(query_output_file),
            }),
            0usize,
        ));

        (Some(query_output), Some(Condvar::new()), Some(output_file))
    } else {
        (None, None, None)
    };

    // ===== 全局聚合 =====
    let global_agg = Mutex::new(GlobalAgg::new(is_paired));

    (0..rayon::current_num_threads())
        .into_par_iter()
        .for_each(|_| {
            #[derive(Copy, Clone)]
            struct QueryColorListItem {
                color: ColorIndexType,
                count: u64,
                next_index: usize,
            }

            let mut queries_colors_list_pool = vec![];
            let mut queries_results = vec![(0u32, 0usize); max_bucket_queries_count];
            let mut temp_colors_list = vec![];
            let mut epoch = 0u32;

            let mut local_rmax: Vec<f64> = Vec::new();
            let mut local_pending: Vec<(usize, Vec<(String, f64)>)> = Vec::new();

            while let Some(input) = {
                let mut lock = buckets_channel.lock();
                lock.pop()
            } {
                epoch = epoch.wrapping_add(1);
                queries_colors_list_pool.clear();
                local_rmax.clear();
                local_pending.clear();

                let start_query_index =
                    input.index as usize * max_bucket_queries_count / buckets_count;

                CompressedBinaryReader::new(
                    &input.path,
                    RemoveFileMode::Remove {
                        remove_fs: !KEEP_FILES.load(Ordering::Relaxed),
                    },
                    DEFAULT_PREFETCH_AMOUNT,
                )
                .decode_all_bucket_items::<QueryColoredCountersSerializer, _>(
                    (Vec::new(), Vec::new()),
                    &mut (),
                    |counters, _| {
                        for query in counters.queries {
                            let (entry_epoch, colors_map_index) = &mut queries_results
                                [query.query_index as usize - start_query_index - 1];

                            if *entry_epoch != epoch {
                                *entry_epoch = epoch;
                                *colors_map_index = usize::MAX;
                            }

                            assert_eq!(counters.colors.len() % 2, 0);

                            for range in counters.colors.chunks(2) {
                                let ColorsRange::Range(range) = ColorsRange::from_slice(range);

                                OPS_COUNT.fetch_add(1, Ordering::Relaxed);
                                COL_COUNT.fetch_add(range.len(), Ordering::Relaxed);

                                for color in range {
                                    queries_colors_list_pool.push(QueryColorListItem {
                                        color,
                                        count: query.count,
                                        next_index: *colors_map_index,
                                    });
                                    *colors_map_index = queries_colors_list_pool.len() - 1;
                                }
                            }
                        }
                    },
                );

                let bucket_index = input.index;

                // JSONL 中间流（按 bucket 写）
                let mut maybe_compressed_stream = if emit_jsonl {
                    Some(CompressedBinaryWriter::new(
                        &temp_dir.join("query-data"),
                        &(
                            get_memory_mode(SwapPriority::ColoredQueryBuckets),
                            CompressedBinaryWriter::CHECKPOINT_SIZE_UNLIMITED,
                            get_compression_level_info(),
                        ),
                        bucket_index as usize,
                        &(),
                    ))
                } else {
                    None
                };

                let mut jsonline_buffer = vec![];

                // 遍历当前 bucket 的 queries
                for (query, mut query_colors_list_index) in
                    queries_results.iter().enumerate().filter_map(|(i, r)| {
                        if r.0 != epoch {
                            None
                        } else {
                            Some((i + start_query_index, r.1))
                        }
                    })
                {
                    // 1) flatten (color, count)
                    temp_colors_list.clear();
                    while query_colors_list_index != usize::MAX {
                        let el = &queries_colors_list_pool[query_colors_list_index];
                        temp_colors_list.push((el.color, el.count));
                        query_colors_list_index = el.next_index;
                    }
                    temp_colors_list.sort_unstable_by_key(|r| r.0);

                    // 2) 聚合同色 => fraction，构建 matches
                    let mut matches: Vec<(String, f64)> = Vec::new();
                    for grp in temp_colors_list.nq_group_by(|a, b| a.0 == b.0) {
                        let color_index = grp[0].0;
                        let color_presence =
                            grp.iter().map(|x| x.1).sum::<u64>();
                        let denom =
                            query_kmers_count[query as usize].max(1) as f64;
                        let frac = (color_presence as f64) / denom;

                        let key = match colored_query_output_format {
                            ColoredQueryOutputFormat::JsonLinesWithNumbers => {
                                color_index.to_string()
                            }
                            ColoredQueryOutputFormat::JsonLinesWithNames => {
                                colormap.get_color_name(color_index, true).to_string()
                            }
                        };

                        matches.push((key, frac));
                    }

                    // 3) 收集 Rmax + pending（后面统一按阈值判）
                    let mut rmax = 0.0f64;
                    for &(_, v) in &matches {
                        if v > rmax {
                            rmax = v;
                        }
                    }
                    local_rmax.push(rmax);
                    local_pending.push((query, matches.clone()));

                    // 4) 可选：写 JSONL 行
                    if let Some(ref mut stream) = maybe_compressed_stream {
                        if emit_jsonl {
                            jsonline_buffer.clear();
                            write!(
                                jsonline_buffer,
                                "{{\"query_index\":{}, \"matches\":{{",
                                query
                            )
                            .unwrap();
                            for (i, (k, v)) in matches.iter().enumerate() {
                                if i != 0 {
                                    write!(jsonline_buffer, ",").unwrap();
                                }
                                write!(jsonline_buffer, "\"{}\": {:.2}", k, v).unwrap();
                            }
                            writeln!(jsonline_buffer, "}}}}").unwrap();
                            stream.write_data(&jsonline_buffer[..]);
                        }
                    }
                }

                // 5) 合并本 bucket 的采样和 pending 到全局
                {
                    let mut g = global_agg.lock();
                    for r in &local_rmax {
                        g.sampler.try_push(*r);
                    }
                    g.pending
                        .extend(local_pending.drain(..));
                }

                // 6) 如果有 JSONL，则把本 bucket 流合并进最终文件
                if let (Some(mut stream), Some(ref condvar), Some(ref query_output)) =
                    (maybe_compressed_stream, output_sync_condvar.as_ref(), maybe_query_output.as_ref())
                {
                    let stream_path = stream.get_path();
                    stream.finalize();

                    let mut decompress_stream = CompressedBinaryReader::new(
                        stream_path,
                        RemoveFileMode::Remove { remove_fs: true },
                        DEFAULT_PREFETCH_AMOUNT,
                    );

                    let mut guard = query_output.lock();
                    while guard.1 != bucket_index {
                        condvar.wait(&mut guard);
                    }

                    let (queries_file, query_write_index) = guard.deref_mut();
                    std::io::copy(
                        &mut decompress_stream.get_single_stream(),
                        queries_file,
                    )
                    .unwrap();
                    *query_write_index += 1;
                    condvar.notify_all();
                }
            }
        });

    // ===== 收尾：根据采样计算 T，处理 pending，处理孤立 mate，输出结果 =====
    let mut g = global_agg.lock();

    let t = g.sampler.compute_threshold();

    // 用最终阈值处理所有 pending reads
    let pending = std::mem::take(&mut g.pending);
    for (idx, matches) in pending {
        let dec = accept_single(idx, &matches, t);
        g.apply_decision(idx, dec);
    }

    // paired 模式下处理没等到 mate 的 read
    g.finalize_orphans();

    // 基名：emit_jsonl=true 用 JSONL 前缀，否则用传入的 output_file
    let base = if let Some(p) = output_path_final.as_ref() {
        p.with_extension("")
    } else {
        output_file.with_extension("")
    };

    let species_out = base.with_extension("species_counts.tsv");
    let discard_out = base.with_file_name(format!(
        "{}_discard_log.tsv",
        base.file_name().unwrap().to_string_lossy()
    ));
    let summary_out = base.with_file_name(format!(
        "{}_summary.txt",
        base.file_name().unwrap().to_string_lossy()
    ));

    // 1) species_counts.tsv
    {
        let mut v: Vec<_> = g
            .species
            .iter()
            .map(|(sp, (cnt, idxs))| (sp.clone(), *cnt, idxs.clone()))
            .collect();
        v.sort_by_key(|t| Reverse(t.1));

        let mut w = File::create(&species_out)?;
        for (sp, cnt, idxs) in v {
            writeln!(
                w,
                "{}\t{}\t{}",
                sp,
                cnt,
                idxs
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
        }
    }

    // 2) discard_log.tsv
    {
        let mut v: Vec<_> = g
            .discard
            .iter()
            .map(|(i, (r, d))| (*i, *r, d.clone()))
            .collect();
        v.sort_by_key(|t| t.0);

        let mut w = File::create(&discard_out)?;
        writeln!(w, "index\treason\tdetail")?;
        for (i, r, d) in v {
            writeln!(w, "{}\t{}\t{}", i, r, d)?;
        }
    }

    // 3) summary.txt
    {
        let kept: usize = g.species.values().map(|(n, _)| *n).sum();
        let mut w = File::create(&summary_out)?;
        writeln!(
            w,
            "threshold(T)={:.6} clipped_to=[{:.6},{:.6}]",
            t, THR_MIN, THR_MAX
        )?;
        writeln!(w, "discard_reasons:")?;
        let mut rs: Vec<_> = g.reason_stats.iter().map(|(k, v)| (*k, *v)).collect();
        rs.sort_by(|a, b| a.0.cmp(&b.0));
        for (k, v) in rs {
            writeln!(w, "  {}: {}", k, v)?;
        }
        writeln!(w, "kept_species_reads={}", kept)?;
    }

    ggcat_logging::info!(
        "Operations count: {} vs real {}",
        OPS_COUNT.load(Ordering::Relaxed),
        COL_COUNT.load(Ordering::Relaxed)
    );

    Ok(())
}
