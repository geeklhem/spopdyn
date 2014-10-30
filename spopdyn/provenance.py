import sys 
import platform
import time
import subprocess

def provenance():
    try:
        version = subprocess.check_output("git describe --always --long",shell=True).strip()
    except Exception:
        version = "Unknown version"
    return (time.asctime(),
            version,
            "|".join(platform.uname()),
            " ".join(sys.argv))

def log(msg):
    with open("log.txt","a+") as f:
        f.write("\t".join(provenance()) +"\t" + msg + "\n" )

def metadata(path):
    ext = path.split(".")[-1]
    if ext == "svg":
        open(path, 'a').write("\n<!-- Metadata: \n{} -->\n".format("\n".join(provenance())))
    elif ext == "eps":
        open(path, 'a').write("\n%Metadata:\n%{}\n".format("\n%".join(provenance())))

    log("Created:{}".format(path))
    
if __name__ == "__main__":
    print provenance()
    subprocess.check_call("touch test.eps",shell=True)
    subprocess.check_call("touch test.svg",shell=True)
    metadata("test.eps")
    metadata("test.svg")

    
