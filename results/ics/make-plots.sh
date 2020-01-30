#!/bin/bash

if [ -d ./venv ]; then
  source ./venv/bin/activate
else
  python3 -m venv venv
  source ./venv/bin/activate
  pip install --upgrade pip
  pip install -r requirements.txt
fi

./tidy.py all.csv -o tidy.csv

./analysis.py baseline tidy.csv -o baseline.csv
./plot.py baseline baseline.csv
./plot.py baseline_bench baseline.csv

./analysis.py expert tidy.csv -o expert.csv
./plot.py expert expert.csv

./analysis.py speeds tidy.csv -o speeds.csv
sed '/slow/d' speeds.csv > speeds-fixed.csv
./plot.py distribution speeds-fixed.csv

./analysis.py marshall tidy.csv -o marshall.csv
./plot.py marshall marshall.csv

out_dir="$1"
shift

if [ -d "$out_dir" ]; then
  cp *.pdf "$out_dir"
fi
