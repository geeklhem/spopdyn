def libpSSA_config(sbml_file, outpath, species,
                   dt=0.01, tend=1000, n=3,
                   spatial=True, gridpoints=10, K=5000):

    if spatial:
        volume = "homogeneous2d" 
    else:
        volume = "single"

    conf = {
        "output_path": outpath,
        "mode": "trial",
        "methods":  0,
        "species": ",".join(species),
        "gnuplot": "no",
        "verbose": "no",
        "dt": dt,
        "tstart": 0,
        "tend": tend,
        "ntrajectories": n,
        "sbml_file": sbml_file,
        "boundary": "periodic",
        "volumetype": volume,
        "gridpoints": gridpoints,
        "initialpop":"multiply",
        "omega":K,
    }

    s = "\n".join(["{} = {}".format(k,v) for k,v in conf.items()])
    return s 


