import argparse
import yaml
#import multiprocessing as mp
import modules.class_fourteen_matricesS
import modules.class_fourteen_matricesP
import modules.class_fourteen_matrices_spinS
import modules.class_fourteen_matrices_spinP
import modules.current_volt
import generate_matrices.system_def
import generate_matrices.system_def_CP2K


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--build", help="Building system", nargs="+")
    parser.add_argument("-bCP2K", "--buildCP2K", help="Building system with CP2K matrices", nargs="+")
    parser.add_argument("-cs", "--config_sequential", help="Config file sequential calculation", nargs="+")
    parser.add_argument("-cp", "--config_parallel", help="Config file parallel calculation", nargs="+")
    parser.add_argument("-scs", "--spincfg_sequential", help="Config file with spin sequential", nargs="+")
    parser.add_argument("-scp", "--spincfg_parallel", help="Config file with spin parallel", nargs="+")
    parser.add_argument("-ivc", "--iv_characteristics", help="Calculate I-V characteristics", nargs="+")
    parser.add_argument("-p", "--printing", help="Printing of files", nargs="+")
    args = parser.parse_args()
    
    
    if args.build:
        for i in args.build:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            generate_matrices.system_def.build_and_dump(cfg)

    elif args.buildCP2K:
        for i in args.buildCP2K:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            generate_matrices.system_def_CP2K.build_and_dump(cfg)

    elif args.config_sequential:
        for i in args.config_sequential:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            system = modules.class_fourteen_matricesS.fourteen_matrices(cfg)
            system.NEGF()
    
    elif args.config_parallel:
        for i in args.config_parallel:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            system = modules.class_fourteen_matricesP.fourteen_matrices(cfg)
            system.NEGF()
    
    elif args.spincfg_sequential:
        for i in args.spincfg_sequential:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            system = modules.class_fourteen_matrices_spinS.fourteen_matrices_spin(cfg)
            system.NEGF()

    elif args.spincfg_parallel:
        for i in args.spincfg_parallel:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            system = modules.class_fourteen_matrices_spinP.fourteen_matrices_spin(cfg)
            #pool = mp.Pool(processes=mp.cpu_count())
            #pool = mp.Pool(2)
            #pool.map(system.NEGF, ["alpha","beta"])
            system.NEGF()
    
    elif args.iv_characteristics:
        for i in args.iv_characteristics:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            system = modules.current_volt.iv_characteristics(cfg)
            system.iv_curves()

    elif args.printing:
        for i in args.printing:
            with open(i,"r") as yamlfile:
                cfg = yaml.load(yamlfile)
            system = modules.class_fourteen_matrices.fourteen_matrices(cfg)
            system.plot()
