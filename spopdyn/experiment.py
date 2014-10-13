import os
import subprocess
import logging
import glob
import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

import spopdyn.sbml_writer
import spopdyn.libpssa_config
import spopdyn.extract
import spopdyn.display

PSSALIB = os.path.dirname(__file__)+"/../spopdyn/libpssa/pssa_cli/pssa"

logger = logging.getLogger("spopdyn")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(levelname)s:%(name)s: %(message)s'))
logger.addHandler(ch)

logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))
logging.addLevelName( logging.INFO, "\033[1;42m%s\033[1;0m" % logging.getLevelName(logging.INFO))
logging.addLevelName( logging.DEBUG, "\033[1;43m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
logging.getLogger().addFilter(logging.Filter("vcontact"))


    
def growth_rate(temperature,habitat,muH,sH,muT,sT,rho):
    l = len(temperature.flat)
    rate = np.zeros(l)
    for x in range(l):
        cov = np.array(([sH**2,rho*sT*sH],[rho*sT*sH,sT**2]))
        rate[x] = st.multivariate_normal.pdf((habitat.flat[x],temperature.flat[x]),
                                             [muH,muT],
                                             cov)
        if rate[x] < 1e-50:
            rate[x] = 0
    return rate


def create_environment(gridpoints,alpha=2,temperature=.5):
    """
    Args: 
        gridpoints (int): number of gridpoints.
        temperature (float): uniform temperature value.
        alpha (float): autocorrelation parameter for habitat.
    Returns:
       (tuple) two gridpoints^2 matrix: environment and temperature.
    """
    temperature = np.zeros((gridpoints,gridpoints)) + temperature
    habitat = spopdyn.landscape.seq_gaussian(size=gridpoints,alpha=alpha,rmax=1)[0]
    return habitat,temperature



def sticky_experiment(habitat,temperature,species,initial_pop,param):
    """Chain multiple experiments and use the final condition of one as
    the initial conditions of the following.  See the documentation of
    experiment to have a description of the parameters. The only
    difference is that habitat & temperature are list of np.arrays
    (giving the conditions for the consecutive experiments)
    """

    if not os.path.exists(param["name"]):
        os.mkdir(param["name"])
        logger.info("Directory {}/ created.".format(param["name"]))    
    else:
        logger.warning("{}/ exists already.".format(param["name"]))


    

    env = {}
    if len(habitat) == 1:
        env["Habitat"] = habitat[0] 
        habitat = habitat * len(temperature)
    elif len(temperature) == 1:
        env["Temperature"] = temperature[0]
        temperature = temperature * len(habitat)


        
    name = param["name"]
    paths = {"data":name+"/data.pkle"}
    if not os.path.exists(paths["data"]):
        data_file = []
        env_path = []
        for n,(h,t) in enumerate(zip(habitat,temperature)):
            env_path.append((np.mean(h),np.mean(t)))
            param["name"] = name+"/step{}".format(n)
            experiment(h,t,species,initial_pop,param)
            data_file.append(param["name"]+"/data.pkle") 
            initial_pop = spopdyn.extract.get_final_pop(param["name"]+"/PDM/PSSA_trajectory_0_0.txt",param["frames"])
        data = spopdyn.extract.fusion_consecutive(data_file)
        with open(paths["data"],"w") as f:
                logger.info("Saving data as {}.".format(paths["data"]))
                pickle.dump((data,env_path),f)
    else:
        with open(paths["data"],"r") as f:
            logger.info("Loading data from {}.".format(paths["data"]))
            data,env_path = pickle.load(f)

    ts_toplot = ("Regional Species Richness",
                 "Regional Diversity",
                 "Local Species Richness",
                 "Local Diversity",
                 "Local CSI",
                 "Local CTI",
                 "Local Biomass")
    ts = [(n,(data[0][n],data[1][n])) for n in ts_toplot]

    dt = len(ts[0][1][0])/len(env_path)
    events = [dt*(i+1) for i in range(len(env_path[:-1]))]
    spopdyn.display.experimental_report([(k,v) for k,v in env.items()],
                                        species,
                                        ts,
                                        env_path,events)
    plt.savefig("{}/experimental_report.pdf".format(name),pad_inches=0)
    plt.savefig("{}/experimental_report.png".format(name),pad_inches=0)
    plt.clf() 
    
    

def experiment(habitat,temperature,species,initial_pop,param):
    """
    Args:
        habitat (np.array): habitat value in 0,1
        temperature (np.array): temperature value in 0,1
        species (np.array): species properties (muH,muT,sH,sT,rho)
        initial_pop (list of np.array): number of individual by species 
        param (dict): parameters (required are alpha,name)
    """

    #------- CHECK PARAMETERS ----------------#
    if (habitat.shape[1] != temperature.shape[0]
        or habitat.shape[1] != temperature.shape[1]
        or habitat.shape[1] != habitat.shape[0]):
        raise ValueError("Habitat and Temperature should be square matrix of the same size")
    else:
        param["gridpts"] = habitat.shape[0]

    #------- CREATE OUTPUT DIRECTORY ----------------#
    if not os.path.exists(param["name"]):
        os.mkdir(param["name"])
        logger.info("Directory {}/ created.".format(param["name"]))    
    else:
        logger.warning("{}/ exists already.".format(param["name"]))

    paths = {
        "sbml": param["name"]+"/model.sbml", 
        "config": param["name"]+"/config.cfg",
        "data": param["name"]+"/data.pkle",
        "ic": param["name"]+"/initial_conditions.pkle",
    }
    param["paths"] = paths 
    
    #------- WRITE MODEL ----------------#
    if not os.path.exists(paths["sbml"]):
        gr = [growth_rate(temperature, habitat,
                          s[0], s[2], s[1], s[3], s[4])
              for s in species]
        #        null_gr = [n for n,x in enumerate(gr) if x.sum()==0]
        #if len(null_gr):
        #    logger.warning("Species {} dropped (null growth rate in every cell)")
        #    gr.pop(null_gr)
        sbml = spopdyn.sbml_writer.CompetitiveLV(gr,
                                                 param["d"],
                                                 param["m"],
                                                 param["alpha"],
                                                 initial_pop,
                                                 param["gridpts"],
                                             )
        sbml.save(paths["sbml"])
        logger.info("SBML file for model '{}' saved as {}.".format(sbml,
                                                                   paths["sbml"]))
        with open(paths["ic"],'w') as f:
                pickle.dump((species, param, habitat, temperature,initial_pop), f)
    else:
        with open(paths["ic"],'r') as f:
            species, param, habitat, temperature,initial_pop = pickle.load(f)

    #------- WRITE SBML CONFIG FILE ----------------#
    if not os.path.exists(paths["config"]):
        cfg = spopdyn.libpssa_config.libpSSA_config(paths["sbml"],
                                                    param["name"]+"/",
                                                    sbml.species,
                                                    dt=param["dt"], tend=param["tend"],n=param["replicas"],
                                                    spatial=True,
                                                    gridpoints=param["gridpts"],
                                                    K=param["K"]*param["gridpts"]**2)
        with open(paths["config"],"w") as f:
            f.write(cfg)
            logger.info("config file saved as {}.".format(paths["config"]))

    #------- RUN LIBRARY ----------------#
    if not os.path.exists(param["name"]+"/libpssa.log"):
        logger.info("Running pSSAlib...")
        try:
            with open(param["name"]+"/libpssa.err",'w') as f:
                re = subprocess.check_output("{} -c {}".format(PSSALIB, paths["config"]),
                                             shell=True,stderr=f)
        except Exception as e:
            logger.error(e)
            raise e
        else:
            with open(param["name"]+"/libpssa.out",'w') as f:
                f.write(re)
                
    #------- DATA EXTRACTION ----------------#                
    logger.info("Data extraction...")
    if not os.path.exists(paths["data"]):
        files = glob.glob("{}/PDM/PSSA_trajectory_*.txt".format(param["name"]))
        data = spopdyn.extract.extract_multiple_trajectories_to_timeseries(files,param["frames"],species)
        with open(paths["data"],"w") as f:
            logger.info("Saving data as {}.".format(paths["data"]))
            pickle.dump(data,f)

    else:
        with open(paths["data"],"r") as f:
            logger.info("Loading data from {}.".format(paths["data"]))
            data = pickle.load(f)

    env = (("Habitat",habitat),
           ("Temperature",temperature),)

    ts_toplot = ("Regional Species Richness",
                 "Regional Diversity",
                 "Local Species Richness",
                 "Local Diversity",
                 "Local CSI",
                 "Local CTI",
                 "Local Biomass")
    ts = [(n,(data[0][n],data[1][n])) for n in ts_toplot] 

    spopdyn.display.experimental_report(env,
                                        species,
                                        ts)
    plt.savefig("{}/experimental_report.pdf".format(param["name"]),pad_inches=0)
    plt.savefig("{}/experimental_report.png".format(param["name"]),pad_inches=0)
    plt.clf() 


