#!/bin/bash
cd ~/project/
commitID=`git log -n 1 --pretty="%h" -- environment.yml`
sed -i '/^# environment.yml/d' Singularity && echo -e "\n# environment.yml commit ID: $commitID\n" >> Singularity
git config --global user.email "alcalan@fellows.iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add dag.png
git add dag.html
git status
git commit -m "circle CI deployment [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/iarcbioinfo/alignment-nf/trigger/5360932c-62ae-4969-9581-312f4089660e/ 


