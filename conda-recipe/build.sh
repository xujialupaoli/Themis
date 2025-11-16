#!/usr/bin/env bash
set -euo pipefail
echo "[Themis] unified build started"

export CMAKE_BUILD_PARALLEL_LEVEL=${CMAKE_BUILD_PARALLEL_LEVEL:-${CPU_COUNT:-1}}
export RUST_BACKTRACE=1

# 避免和 conda 注入 flags 冲突，只补充必要的
export CPPFLAGS="${CPPFLAGS:-} -I${PREFIX}/include"
export LDFLAGS="${LDFLAGS:-} -L${PREFIX}/lib"

# -----------------------------------------------------------------------------
# 1) Build & install ganon C++  (产生 ELF: ganon-build / ganon-classify)
# -----------------------------------------------------------------------------
echo "[Themis] building ganon (C++) ..."
pushd "${SRC_DIR}/thirdparty/ganon_mod"
rm -rf build_cpp
cmake -S . -B build_cpp -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCONDA=ON \
  -DVERBOSE_CONFIG=ON
cmake --build build_cpp --parallel "${CMAKE_BUILD_PARALLEL_LEVEL}"
cmake --install build_cpp
popd

# 确认 ELF 存在
file "${PREFIX}/bin/ganon-classify"
file "${PREFIX}/bin/ganon-build"

# -----------------------------------------------------------------------------
# 2) Install ganon Python 包（提供 `ganon` 顶层 CLI 以及 `ganon-report`）
# -----------------------------------------------------------------------------
echo "[Themis] installing ganon (Python) ..."
pushd "${SRC_DIR}/thirdparty/ganon_mod"
# 用宿主 env 的 setuptools/pip，避免隔离构建拉多余依赖
"${PYTHON}" -m pip install . --no-deps --no-build-isolation -vv
popd

# 加一个明确的 shim，确保 ganon-report 可执行（有的发行版没单独装脚本）
cat > "${PREFIX}/bin/ganon-report" <<'SH'
#!/usr/bin/env bash
exec "${CONDA_PREFIX}/bin/python" -m ganon.report "$@"
SH
chmod +x "${PREFIX}/bin/ganon-report"

# -----------------------------------------------------------------------------
# 3) Build & install ggcat (Rust)
# -----------------------------------------------------------------------------
echo "[Themis] building ggcat (Rust) ..."
pushd "${SRC_DIR}/thirdparty/ggcat_mod"
# crates/cmdline/ 是 ggcat 主 CLI crate
cargo install --locked --root "${PREFIX}" --path crates/cmdline/
# 生成第三方依赖清单（不失败）
cargo-bundle-licenses --format yaml --output "${SRC_DIR}/THIRDPARTY.yml" || true
popd

# -----------------------------------------------------------------------------
# 4) Install Python package 'themis'（只暴露一个顶层 CLI: themis）
# -----------------------------------------------------------------------------
echo "[Themis] installing Themis (Python) ..."
pushd "${SRC_DIR}"
"${PYTHON}" -m pip install . --no-deps --no-build-isolation -vv
popd

# sanity check：import themis
"${PREFIX}/bin/python" - <<'PY'
import importlib, sys
m = importlib.import_module("themis")
print("[sanity] themis module:", m.__file__)
PY

echo "[Themis] unified build finished"
