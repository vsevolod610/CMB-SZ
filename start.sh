rm data/N/result*
cd core
python check.py
python start_analyze.py
python estimate.py > ../res/T0.txt
cd ..
