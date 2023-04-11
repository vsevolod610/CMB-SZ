rm data/set/result*
cd core
python check.py
python create_data.py
python start_analyze.py
python estimate.py > ../data/set/res/T0.txt
cd ..
