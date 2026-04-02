#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "============================================================"
echo "RNA Hypothesis Workbench"
echo "============================================================"
echo

if command -v ss >/dev/null 2>&1; then
  if ss -ltn 2>/dev/null | grep -q ':3838 '; then
    echo "The workbench already appears to be running on http://127.0.0.1:3838"
    echo
    echo "The launcher will open that existing app instead of starting a second copy."
    echo
    if command -v xdg-open >/dev/null 2>&1; then
      xdg-open "http://127.0.0.1:3838" >/dev/null 2>&1 || true
    fi
    read -r -p "Press Enter to exit..."
    exit 0
  fi
elif command -v netstat >/dev/null 2>&1; then
  if netstat -ltn 2>/dev/null | grep -q ':3838 '; then
    echo "The workbench already appears to be running on http://127.0.0.1:3838"
    echo
    echo "The launcher will open that existing app instead of starting a second copy."
    echo
    if command -v xdg-open >/dev/null 2>&1; then
      xdg-open "http://127.0.0.1:3838" >/dev/null 2>&1 || true
    fi
    read -r -p "Press Enter to exit..."
    exit 0
  fi
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found in this environment."
  echo
  echo "Please open:"
  echo "  ${ROOT_DIR}/START_HERE.txt"
  echo
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "${ROOT_DIR}/START_HERE.txt" >/dev/null 2>&1 || true
  fi
  read -r -p "Press Enter to exit..."
  exit 1
fi

echo "Starting the app..."
echo
echo "This terminal will stay open while the app runs."
echo "If the browser does not open, visit http://127.0.0.1:3838"
echo

cd "${ROOT_DIR}"
Rscript deseq2_workbench.R

echo
read -r -p "The app has stopped. Press Enter to exit..."
