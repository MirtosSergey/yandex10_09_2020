import PreSolver as PreS
import PostSolver as PosS
import MathModel as MM

if __name__ == "__main__":
    config = PreS.read_config()
    PreS.rewrite_res()
    res = MM.Calc_one_out_path()
    DF = PosS.write_res(res)
    #PosS.write_graf(DF, 1)