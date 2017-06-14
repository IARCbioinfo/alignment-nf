#!/bin/bash
cd ~/alignment-nf/
git add dag.png
git add dag.html
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/iarcbioinfo/alignment-nf/trigger/5360932c-62ae-4969-9581-312f4089660e/ 


