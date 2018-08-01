#!/usr/bin/env bash

git tag -a etomica_modules_$(date +%Y%m%d%H%M%S) -m "Etomica modules build tag"
echo "Tag created, push with git push --follow-tags"
