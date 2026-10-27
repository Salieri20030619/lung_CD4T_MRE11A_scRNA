#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
RUN="$ROOT/02_analysis/GSE174725/suppC2_20260911_164417_l173TC"
CONDA=/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda
export RUN

exec 9>"$ROOT/tmp/GSE174725_suppC2.lock"
if ! flock -n 9; then echo SUPPC2_ALREADY_RUNNING; exit 0; fi
trap 'rc=$?; echo "GSE174725_SUPPC2_FAIL EXIT_CODE=$rc"; exit "$rc"' ERR

test -s "$RUN/input/counts.mtx"
test -s "$RUN/input/cells.tsv.gz"
cp -n "$RUN/score.R" "$RUN/score.R.before_boolfix"

sed -i '/^ take<-m\$sample==ss /c\ take<-m$sample==ss \& if(scope=="Core") {flag<-tolower(as.character(m$core));stopifnot(all(flag %in% c("true","false")));flag=="true"} else rep(TRUE,nrow(m))' "$RUN/score.R"

"$CONDA" run --no-capture-output -n R_4.5.2 Rscript -e 'invisible(parse(file=file.path(Sys.getenv("RUN"),"score.R")));m<-read.delim(gzfile(file.path(Sys.getenv("RUN"),"input/cells.tsv.gz")));flag<-tolower(as.character(m$core));stopifnot(all(flag %in% c("true","false")),sum(flag=="true")==2958);cat("R_PARSE_PASS\nBOOLEAN_CHECK_PASS\n")'
"$CONDA" run --no-capture-output -n batch_bench python -m py_compile "$RUN/summarize.py"
echo SYNTAX_CHECK=PASS

"$CONDA" run --no-capture-output -n R_4.5.2 Rscript "$RUN/score.R" "$RUN"
"$CONDA" run --no-capture-output -n batch_bench python -u "$RUN/summarize.py" "$RUN"

printf 'GSE174725_SUPPC2_PASS\n' > "$RUN/STATUS.txt"
PACK="$ROOT/GSE174725_suppC2_return.tar.gz"
tar -C "$RUN" -czf "$PACK.part" \
  tables figures input/cells.tsv.gz input/genes.tsv \
  input/expected_samples.tsv input/provenance.txt \
  export.py score.R summarize.py captions.txt methods.txt STATUS.txt
mv "$PACK.part" "$PACK"
echo "RETURN_PACKAGE=$PACK"
echo GSE174725_SUPPC2_PASS
