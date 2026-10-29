#!/bin/bash
# After the single-satellite atlas runs: compact their releases, build the snpp/n20/n21 comparison layers, redeploy.
AX=/c/Users/jelle/Desktop/PCB_Analysis_1/INPUT_DATA/atlas_extension
MAP=/c/Users/jelle/Desktop/PCB_Analysis_1/PCB_exploratory
until grep -q "all single-satellite queues done" $AX/logs/run_multi_all.log; do sleep 120; done
cd $AX
for v in snpp n20 n21; do
  for y in $(ls release_$v | grep -E '^[0-9]{4}$'); do
    [ -d release_$v/$y/compact ] || python -u scripts/compact_release.py $y --workers 1 --release release_$v
  done
done
cd $MAP
python -u build_compare_map.py --layers snpp n20 n21
cp compare/*.json.gz .deploy_tmp/compare/ && cp compare/index.json .deploy_tmp/compare/index.json
cd .deploy_tmp && vercel deploy --prod --yes --archive=tgz
echo "FINISHED $(date)"
