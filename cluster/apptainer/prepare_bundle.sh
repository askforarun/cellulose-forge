#!/bin/bash
set -euo pipefail
repo="$(git rev-parse --show-toplevel)"
here="${repo}/cluster/apptainer"
ff_source="${CELLULOSE_FORGE_CHARMM36_DIR:-${HOME}/.local/share/cellulose-forge/forcefields/charmm36-feb2026_cgenff-5.0.ff}"
git -C "${repo}" diff --quiet
git -C "${repo}" diff --cached --quiet
test -d "${ff_source}"
git -C "${repo}" bundle create "${here}/cellulose-forge.gitbundle" HEAD
rm -rf "${here}/charmm36-feb2026_cgenff-5.0.ff"
cp -a "${ff_source}" "${here}/charmm36-feb2026_cgenff-5.0.ff"
echo "Prepared private bundle at ${here}"
echo "Do not publish it: it contains the locally installed CHARMM36 distribution."

