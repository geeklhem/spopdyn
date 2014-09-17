import spopdyn.sbml_writer
import os 
F = {} # Fixtures 
def setup():
    F["compLVargs"] = ([1]*3,2,1,0.5)
    F["file"] = "test.sbml"
def test_lv():
    spopdyn.sbml_writer.CompetitiveLV(*F["compLVargs"]) 

def test_write():
    a= spopdyn.sbml_writer.CompetitiveLV(*F["compLVargs"])
    a.save(F["file"])


def teardown():
    os.remove(F["file"])
