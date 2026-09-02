#!/usr/bin/env bash
# Build and install regInspect into jamovi desktop and/or a running jamovi
# Docker container.
#
#   bash install.sh              both targets, whichever are available
#   bash install.sh desktop
#   bash install.sh docker [container]     (default container: jamovi)
#
# Unlike jamovi-corr/jamovi-jmvplus's install.sh, the desktop target here
# does NOT set a custom R_LIBS_USER ($HOME/R/.Rlib-arm etc.) -- that path
# doesn't exist on an R installation managed by rig. Instead it just uses
# whichever R `Rscript` resolves to (respecting ~/.Rprofile, which appends
# jamovi.app's bundled module library -- see the user's ~/.Rprofile and
# jamovi-linreg/task_plan.md's "Errors log" for why). Pick the active R
# version with `rig default <version>` before running this if needed; it
# must match the R version jamovi.app itself bundles, or jmvcore segfaults
# on load.
set -euo pipefail

TARGET="${1:-both}"
CONTAINER="${2:-jamovi}"

HERE="$(cd "$(dirname "$0")/../regInspect" && pwd)"
MODULE=regInspect
VERSION="$(awk -F': *' '$1 == "Version" { print $2; exit }' "$HERE/DESCRIPTION")"
ARTIFACT="$HERE/${MODULE}_${VERSION}.jmo"

# ── desktop ──────────────────────────────────────────────────────────────────
install_desktop() {
  local APP APP_R LOG
  APP=/Applications/jamovi.app
  APP_R="$APP/Contents/Frameworks/R.framework/Versions/Current/Resources/bin/R"
  [ -x "$APP_R" ] || { echo "error: no R inside $APP" >&2; return 1; }

  echo ">> desktop: building regInspect with $(Rscript -e 'cat(R.version.string)')"
  cd "$HERE"

  LOG="$(mktemp)"
  Rscript -e 'jmvtools::install()' 2>&1 | tee "$LOG" | grep -vE '^\s*$' || true

  # jmvtools::install() can report errors on stdout while exiting successfully.
  # It can also claim installation succeeded after a SingletonLock failure.
  [ -f "$ARTIFACT" ] || {
    echo "error: jmvtools did not produce $ARTIFACT" >&2
    rm -f "$LOG"; return 1
  }
  if grep -q 'SingletonLock' "$LOG"; then
    echo
    echo "!! jamovi.app could not be driven (SingletonLock denied)."
    echo "!! The .jmo was still built. Install it by hand:"
    echo "!!   jamovi -> Modules -> jamovi library -> Sideload -> $ARTIFACT"
    rm -f "$LOG"
    return 0
  fi
  if ! grep -q 'Module installed successfully' "$LOG"; then
    echo "error: jmvtools::install() did not install the module (see above)" >&2
    rm -f "$LOG"; return 1
  fi
  rm -f "$LOG"

  local MODDIR="$HOME/Library/Application Support/jamovi/modules/$MODULE"
  if [ -d "$MODDIR" ]; then
    echo ">> desktop: installed at $MODDIR"
  else
    echo "!! desktop: install reported success but $MODDIR does not exist."
    echo "!! Install $ARTIFACT by hand (Modules -> jamovi library -> Sideload)."
    return 1
  fi
}

# ── docker ───────────────────────────────────────────────────────────────────
install_docker() {
  if ! docker ps --format '{{.Names}}' 2>/dev/null | grep -qx "$CONTAINER"; then
    echo "!! docker: container '$CONTAINER' is not running — skipping"
    return 0
  fi

  if ! docker exec "$CONTAINER" sh -c 'command -v jmc >/dev/null 2>&1'; then
    echo "!! docker: jmc is not in the container." >&2
    echo "!! Install the jamovi compiler in the image before using this target." >&2
    return 1
  fi

  echo ">> docker: copying source into $CONTAINER"
  # --no-mac-metadata/--no-xattrs: AppleDouble ._ files otherwise land in the
  # container and jmc tries to compile them.
  tar --no-mac-metadata --no-xattrs -C "$HERE" -cf - DESCRIPTION NAMESPACE R jamovi inst \
    | docker exec -i "$CONTAINER" sh -c \
        'rm -rf /tmp/regInspect-src && mkdir -p /tmp/regInspect-src && tar -C /tmp/regInspect-src -xf -'

  echo ">> docker: jmc --install"
  docker exec -i "$CONTAINER" bash -s <<'INCONTAINER'
set -euo pipefail
source /usr/lib/jamovi/bin/env.conf 2>/dev/null || true
RHOME="${R_HOME:-$(R RHOME 2>/dev/null || true)}"
[ -n "$RHOME" ] || { echo "   error: no R in the container" >&2; exit 1; }
RLIBS=/usr/lib/jamovi/modules/base/R

# --skip-deps assumes car/emmeans/multcomp/gridExtra already resolve from
# RLIBS. jamovi bundles ggplot2/dplyr/tidyr/tibble for its own modules in
# base/R, but car/emmeans/multcomp/gridExtra come from jmv's OWN module
# library -- if the smoke test below fails on "there is no package called
# 'car'" (or emmeans/multcomp/gridExtra), add
# /usr/lib/jamovi/modules/jmv/R to RLIBS too, or install jmv first.
jmc --install /tmp/regInspect-src \
    --to /usr/lib/jamovi/modules \
    --rhome "$RHOME" \
    --rlibs "$RLIBS" \
    --patch-version --skip-deps

[ -f /usr/lib/jamovi/modules/regInspect/jamovi.yaml ] || {
  echo "   error: jmc did not install regInspect" >&2; exit 1; }
INCONTAINER

  echo ">> docker: restarting $CONTAINER to load the module"
  docker restart "$CONTAINER" >/dev/null
  echo ">> docker: testing model coefficients against a plain lm()"
  docker exec -i "$CONTAINER" bash -s <<'INCONTAINER'
set -euo pipefail
Rscript --vanilla -e '
    .libPaths(c(
        "/usr/lib/jamovi/modules/jmv/R",
        "/usr/lib/jamovi/modules/base/R",
        "/usr/lib/jamovi/modules/regInspect/R",
        .libPaths()
    ))
    library(regInspect)

    data(mtcars)
    results <- regInspect::regInspect(
        data = mtcars, dep = mpg, covs = c("wt", "hp"),
        blocks = list(list("wt", "hp"))
    )
    tbl <- results$models[[1]]$coef$asDF
    oracle <- coef(lm(mpg ~ wt + hp, data = mtcars))
    stopifnot(isTRUE(all.equal(tbl$est, unname(oracle), tolerance = 1e-6)))
    cat(sprintf("   coefficients smoke test passed: wt = %.3f, hp = %.4f\n",
                tbl$est[2], tbl$est[3]))
'
INCONTAINER
  echo ">> docker: testing partial-residual points on the emm plot"
  docker exec -i "$CONTAINER" bash -s <<'INCONTAINER'
set -euo pipefail
Rscript --vanilla -e '
    .libPaths(c(
        "/usr/lib/jamovi/modules/jmv/R",
        "/usr/lib/jamovi/modules/base/R",
        "/usr/lib/jamovi/modules/regInspect/R",
        .libPaths()
    ))

    data(mtcars)
    options <- regInspect:::regInspectOptions$new(
        dep = "mpg", covs = "wt", blocks = list(list("wt")),
        emMeans = list(list("wt")), emmPoints = "partial"
    )
    analysis <- regInspect:::regInspectClass$new(options = options, data = mtcars)
    analysis$run()

    priv <- analysis$.__enclos_env__$private
    pts <- priv$.emmPointData(1, c("wt", NA, NA), c(jmvcore::toB64("wt"), NA, NA),
                               list(x = jmvcore::toB64("wt")))
    stopifnot(isTRUE(all.equal(sort(pts$.y), sort(mtcars$mpg), tolerance = 1e-6)))
    cat("   partial-residual smoke test passed (single predictor == raw y)\n")
'
INCONTAINER
  echo ">> docker: installed regInspect; open it in jamovi to verify the panel and plots"
}

case "$TARGET" in
  desktop) install_desktop ;;
  docker)  install_docker ;;
  both)    install_desktop || true; echo; install_docker || true ;;
  *)       echo "usage: install.sh [desktop|docker|both] [container]" >&2; exit 1 ;;
esac
