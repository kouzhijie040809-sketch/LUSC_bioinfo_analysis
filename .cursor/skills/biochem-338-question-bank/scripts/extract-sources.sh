#!/usr/bin/env bash
# Extract text from sdu/生化资料/338（生物化学）资料 into _extracted/
# Usage: .cursor/skills/biochem-338-question-bank/scripts/extract-sources.sh
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "$ROOT" ]]; then
  ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
fi
SRC="$ROOT/sdu/生化资料/338（生物化学）资料"
OUT="$SRC/_extracted"
mkdir -p "$OUT"

if [[ ! -d "$SRC" ]]; then
  echo "missing $SRC" >&2
  exit 1
fi

count=0
while IFS= read -r -d '' f; do
  rel="${f#"$SRC"/}"
  base="$(basename "$f")"
  [[ "$base" == "README.md" || "$base" == ".gitkeep" ]] && continue
  dest="$OUT/${rel}.txt"
  mkdir -p "$(dirname "$dest")"
  ext="${base##*.}"
  ext_lc="$(printf '%s' "$ext" | tr 'A-Z' 'a-z')"
  case "$ext_lc" in
    pdf)
      if command -v pdftotext >/dev/null 2>&1; then
        pdftotext -layout "$f" "$dest"
      elif command -v python3 >/dev/null 2>&1; then
        python3 - "$f" "$dest" <<'PY'
import sys
src, dest = sys.argv[1], sys.argv[2]
try:
    from pypdf import PdfReader
    reader = PdfReader(src)
    parts = []
    for i, page in enumerate(reader.pages, 1):
        parts.append(f"\n\n===== PAGE {i} =====\n")
        parts.append(page.extract_text() or "")
    open(dest, "w", encoding="utf-8").write("".join(parts))
except Exception as e:
    open(dest, "w", encoding="utf-8").write(f"[extract failed] {e}\n")
    sys.exit(0)
PY
      else
        echo "[extract skipped: no pdftotext/python] $rel" > "$dest"
      fi
      ;;
    txt|md|markdown|csv)
      cp "$f" "$dest"
      ;;
    *)
      echo "[binary or unsupported; open with agent Read tool] $rel" > "$dest"
      ;;
  esac
  echo "extracted: $rel"
  count=$((count + 1))
done < <(find "$SRC" -type f ! -path '*/_extracted/*' -print0)

echo "done: $count file(s) -> $OUT"
