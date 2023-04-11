# CMB-SZ: python-only

### Package dependeces
You need install `python` (`python3`) and python-packages below:
- `numpy` 
- `scipy`
- `matplotlib`
- `chainconsumer`
- `emcee`
- `astropy`
- `tdqm` (optionally)

(You can install it with `pip`)

### Run code
##### To analyze one cluster
1. In `./data`: put `SZ_data.txt` and `prior.dat`;
    set paramenters in `startSZ.sss`
2. Run `python show/analyze_one.py`
3. Look for the result in puctures `./data`

##### To analyze multiple clusters
1. Put `szdata1.txt`, `szdata2.txt`, ... in `./data/N/szdata`  
    and `prior1.dat`, `prior2.dat`, ... in `./data/N/priors`  
    set paramenters in `./data/startSZ.sss`
2. Run `python ./core/start_analyze.py`
    and look for the result in `./data/N/pic_*`
3. To obtain common $T_0$ estimation run `python ./core/estimate.py`

or just run `bash start.sh` from this directory!

##### To generate artificial data
1. Run `python ./core/create_data.py`
