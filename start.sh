rm data/N/result*
cd core
python start_analyze.py
python estimate.py > ../res/T0.txt
cd ..
