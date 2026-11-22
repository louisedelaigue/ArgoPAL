#!/bin/sh
# Update the coverage report and badge in .misc/ (run from the repository
# root, with internet access: the GDAC tests are included).
set -e
python -m pytest -q -m "network or not network" --cov=argopal --cov-report=xml:.misc/coverage.xml
python -m coverage report > .misc/coverage.txt
genbadge coverage -l -i .misc/coverage.xml -o .misc/coverage.svg
rm .misc/coverage.xml
cat .misc/coverage.txt
