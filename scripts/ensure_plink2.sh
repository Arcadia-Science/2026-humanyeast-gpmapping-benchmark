#!/bin/bash
# ensure_plink2.sh — Make sure plink2 is available before running any plink2 commands.
#
# Resolution order:
#   1. plink2 already on PATH                  → nothing to do
#   2. ./bin/plink2 exists                     → add ./bin/ to PATH and done
#   3. conda install -c bioconda plink2        → done if successful
#   4. Direct download from cog-genomics.org   → installed to ./bin/plink2
#
# Direct-download version:
#   Set PLINK2_VERSION below (or via environment variable) to match the
#   filenames listed on https://www.cog-genomics.org/plink/2.0/
#   Example filenames from that page:
#     plink2_linux_avx2_20260311.zip
#     plink2_mac_arm64_20260311.zip
#   → PLINK2_VERSION=20260311
#
# Usage:
#   bash scripts/ensure_plink2.sh              # run standalone
#   source scripts/ensure_plink2.sh            # run and inherit PATH change
#   Called automatically by scripts/plink.sh

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────

PLINK2_VERSION="${PLINK2_VERSION:-20260311}"    # date tag from the download filename
PLINK2_BIN_DIR="${PLINK2_BIN_DIR:-${PWD}/bin}" # where to install if downloading

# ── Already on PATH? ──────────────────────────────────────────────────────────

if command -v plink2 &>/dev/null; then
    echo "[ensure_plink2] plink2 found: $(command -v plink2)"
    exit 0
fi

# ── Local bin/ install? ───────────────────────────────────────────────────────

if [[ -x "${PLINK2_BIN_DIR}/plink2" ]]; then
    echo "[ensure_plink2] Found local binary: ${PLINK2_BIN_DIR}/plink2"
    export PATH="${PLINK2_BIN_DIR}:${PATH}"
    exit 0
fi

echo "[ensure_plink2] plink2 not found. Attempting installation..."

# ── Try conda ─────────────────────────────────────────────────────────────────

if command -v conda &>/dev/null; then
    echo "[ensure_plink2] Trying: conda install -y -c bioconda -c conda-forge plink2"
    if conda install -y -c bioconda -c conda-forge plink2; then
        if command -v plink2 &>/dev/null; then
            echo "[ensure_plink2] plink2 installed via conda: $(command -v plink2)"
            exit 0
        fi
        echo "[ensure_plink2] conda ran but plink2 still not on PATH — falling back."
    else
        echo "[ensure_plink2] conda install failed — falling back to direct download."
    fi
fi

# ── Direct download from cog-genomics ─────────────────────────────────────────

OS=$(uname -s)
ARCH=$(uname -m)

case "${OS}" in
    Linux)
        case "${ARCH}" in
            x86_64)
                if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
                    PLATFORM="linux_avx2"
                else
                    PLATFORM="linux_x86_64"
                fi
                ;;
            aarch64|arm64)
                PLATFORM="linux_arm64"
                ;;
            *)
                echo "[ensure_plink2] ERROR: Unsupported Linux architecture: ${ARCH}" >&2
                echo "  Download manually from https://www.cog-genomics.org/plink/2.0/" >&2
                exit 1
                ;;
        esac
        ;;
    Darwin)
        case "${ARCH}" in
            arm64)  PLATFORM="mac_arm64"  ;;
            x86_64) PLATFORM="mac_x86_64" ;;
            *)
                echo "[ensure_plink2] ERROR: Unsupported macOS architecture: ${ARCH}" >&2
                exit 1
                ;;
        esac
        ;;
    *)
        echo "[ensure_plink2] ERROR: Unsupported OS: ${OS}" >&2
        exit 1
        ;;
esac

PLINK2_RELEASE_TAG="${PLINK2_RELEASE_TAG:-v2.0.0-a.6.33}"
PLINK2_URL="https://github.com/chrchang/plink-ng/releases/download/${PLINK2_RELEASE_TAG}/plink2_${PLATFORM}.zip"
echo "[ensure_plink2] Downloading: ${PLINK2_URL}"
echo "[ensure_plink2] If this URL fails, check https://www.cog-genomics.org/plink/2.0/"
echo "               and update PLINK2_VERSION at the top of this script."

mkdir -p "${PLINK2_BIN_DIR}"
PLINK2_TMPDIR=$(mktemp -d)
trap 'rm -rf "${PLINK2_TMPDIR}"' EXIT

if command -v wget &>/dev/null; then
    wget -q --show-progress -O "${PLINK2_TMPDIR}/plink2.zip" "${PLINK2_URL}"
elif command -v curl &>/dev/null; then
    curl -fSL -o "${PLINK2_TMPDIR}/plink2.zip" "${PLINK2_URL}"
else
    echo "[ensure_plink2] ERROR: Neither wget nor curl is available." >&2
    exit 1
fi

unzip -q "${PLINK2_TMPDIR}/plink2.zip" plink2 -d "${PLINK2_BIN_DIR}"
chmod +x "${PLINK2_BIN_DIR}/plink2"
export PATH="${PLINK2_BIN_DIR}:${PATH}"

echo "[ensure_plink2] Installed to: ${PLINK2_BIN_DIR}/plink2"
plink2 --version
