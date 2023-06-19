#!/usr/bin/bash

'''
Run ntSynt demo to test installation
'''

set -eux -o pipefail
echo "Running ntSynt Demo..."
echo "This demo assumes that you have the 'ntSynt.py' script in your PATH"

echo "Running ntSynt with 2 input genomes"
ntSynt.py celegans-chrII-III.fa celegans-chrII-III.A.fa --prefix celegans-A-ntSynt

echo "Running ntSynt with 3 input genomes"
ntSynt.py celegans-chrII-III.fa celegans-chrII-III.A.fa celegans-chrII-III.B.fa -k20 --prefix celegans-A-B-ntSynt

echo "DONE! Compare your output files with those in 'expected_result' to ensure that your installation is generating the expected results"
