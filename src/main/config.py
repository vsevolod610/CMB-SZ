import paths # manage path imports by paths.py
from core_sz.model import *

N = 10
case = "100"

models = {'000': model_000, 
          '100': model_100, 
          '112': model_112, 
          '212': model_212 }


if case[0] == '1': # case == '1..'
    path_startSZ = paths.src + '/main/input/startSZ_1xx.sss'
    params_names = [r"$T_0$, K", r"$T_e$, K", r"$\beta$", r"$\Tau$"]

    params_include = lambda T0, Te, beta, Tau, z=0: [T0, Te, beta, Tau]
    def modelT(T0, z): return T0 + (0 * z)


if case[0] == '2': # case == '2..'
    path_startSZ = paths.src + '/main/input/startSZ_2xx.sss'
    params_names = [r"$T_z$, K", r"$T_e$, K", r"$\beta$", r"\Tau"]

    params_include = lambda T0, Te, beta, Tau, z: [T0 * (1 + z), Te, beta, Tau]
    def modelT(T0, z): return T0 * (1 + z)


if case[0] == '0': # case == '0..'
    path_startSZ = paths.src + '/main/input/startSZ_0xx.sss'
    params_names = [r"$T_0$, K", r"$T_e$, K"]

    params_include = lambda T0, Te, beta, Tau, z=0: [T0, Te]
    def modelT(T0, z): return T0 + (0 * z)


config = {
        "N": N, 
        "case": case,
        "model": models[case], 
        'path_startSZ': path_startSZ,
        "params_names": params_names,
        "params_include": params_include,   # create
        "modelT": modelT                    # estimate
        }


if __name__ == "__main__":
    print(f"{N = }")
    print(f"{case = }")
    print(f"startSZ = {path_startSZ}:\n")
    with open(path_startSZ, 'r') as file:
        print(*file)


