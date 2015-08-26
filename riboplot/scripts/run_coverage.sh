#!/bin/bash
# Run coverage and geneated HTML reports
echo "INFO: If this fails with command not found, please install the coverage python package and re-run"
coverage run --source tests/test_riboplot.py --source riboplot setup.py test
coverage html
