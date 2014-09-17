import logging
import itertools

import libsbml
import numpy as np

logger = logging.getLogger("spopdyn")

class SBMLwriter(object):
    """Create a generic SBML object with an empty model and units.  Ready
    to be extended.

    """
    def __init__(self):
        name = "Generic model"
        filename = None
        try:
            self.document = libsbml.SBMLDocument(2,4)
        except ValueError:
            logger.error("Could not create SBMLDocument object")

        self.model = self.document.createModel()

        # Units 
        self.model.setTimeUnits("second")
        self.model.setSubstanceUnits("item")

        per_second = self.model.createUnitDefinition()
        per_second.setId('per_second')
        unit = per_second.createUnit()
        unit.setKind(libsbml.UNIT_KIND_SECOND)
        unit.setExponent(-1)

        # Compartment
        compartment = self.model.createCompartment()
        compartment.setId('c')
        compartment.setConstant(True)
        compartment.setSize(0)
        compartment.setUnits('volume')

        self.model.setId(self.__repr__())

    def save(self,filename=None,set_filename=True):
        if filename == None:
            if self.filename == None:
                raise ValueError
            else:
                filename = self.filename
        libsbml.writeSBMLToFile(self.document,filename)

        if set_filename == True:
            self.filename = filename

    def __repr__(self):
        return self.name


class CompetitiveLV(SBMLwriter): 
    def __init__(self,r,d,alpha=1,init_prop=0.5):
        """ Cosntructor
        Args:
            r (array): specific growth rate.
            d (array or floar): diffusion constant.
            alpha (matrix or float): interspecific competition strength.
            init_prop (array or float): initial proportion (with respect
                to the carrying capacity) of each species.
        """
        
        self.n = len(r)
        self.species = ["SP{{:0{}}}".format(int(np.log10(self.n)+1)).format(i) for i in range(self.n)]
        self.name = "Competitive Lotka Volterra with {} species".format(self.n)
        SBMLwriter.__init__(self)

        # Turn alpha into a matrix:
        if type(alpha) == int or type(alpha)== float :
            a = np.zeros((self.n,self.n)) + alpha
        # Turn d into an array: 
        if type(d) == float or type(d) == int:
            d = np.zeros(self.n) + d
        # Turn init_prop into an array:
        if type(init_prop) == float or type(init_prop)== int:
            init_prop = np.zeros(self.n) + init_prop

        # Add diffusion constants.
        self.model = self.diffusion(self.model, self.species, d)

        # Reproduction and intraspecific competition 
        for sp,ri,ip in zip(self.species,r,init_prop):
            self.model = self.add_sp(self.model,sp,ri,ip)

        # Interspecific competition. 
        for sp1,sp2 in itertools.combinations(self.species,2):
            n1 = int(sp1[2:])
            n2 = int(sp2[2:])
            self.model = self.inter_spe_comp(self.model, sp1, sp2, a[n1,n2])
            self.model = self.inter_spe_comp(self.model, sp2, sp1, a[n2,n1])

        if self.document.getNumErrors() != 0:
            logger.error(self.document.printErrors())
            raise IOError
        
    def diffusion(self,model,species,d):
        diffusion = ['<libpSSA:diffusion xmlns:libpSSA="uri">']
        for sp,di in zip(species,d):
            diffusion.append('<libpSSA:d{} libpSSA:diffusion="{}"/>'.format(sp,di))
        diffusion.append('</libpSSA:diffusion>')
        model.getListOfSpecies().appendAnnotation("\n".join(diffusion))
        return model 
        
        
    def add_sp(self, model, name, r, initial_ammount):
        #Species
        s = model.createSpecies()
        s.setId(name)
        s.setName(name) 
        s.setCompartment('c') 
        s.setConstant(False)
        s.setInitialAmount(initial_ammount)
        s.setSubstanceUnits('item')

        # Growth rate
        gr = model.createParameter()
        gr.setId('r_{}'.format(name))
        gr.setName('r_{}'.format(name))
        gr.setValue(r)
        gr.setUnits('per_second')

        # Reproduction
        reproduction = model.createReaction()
        reproduction.setId('reproduction_'+name) 
        reproduction.setReversible(False)

        kinetic_law = reproduction.createKineticLaw()
        kinetic_law.setMath(libsbml.parseL3Formula("1*"+gr.getId()))

        reactant = reproduction.createReactant()
        reactant.setSpecies(s.getId())
        reactant.setStoichiometry(1)

        product = reproduction.createProduct()
        product.setSpecies(s.getId())
        product.setStoichiometry(2)

        # Intraspecific competition 
        comp = model.createReaction()
        comp.setId('intra_spe_comp_'+name) 
        comp.setReversible(False)

        kinetic_law = comp.createKineticLaw()
        kinetic_law.setMath(libsbml.parseL3Formula("1*"+gr.getId()))

        reactant = comp.createReactant()
        reactant.setSpecies(s.getId())
        reactant.setStoichiometry(2)

        product = comp.createProduct()
        product.setSpecies(s.getId())
        product.setStoichiometry(1)

        return model

    def inter_spe_comp(self,model,s1,s2,a):
        comp = model.createReaction()
        comp.setId('intra_spe_comp_'+s1+"_"+s2) 
        comp.setReversible(False)

        kinetic_law = comp.createKineticLaw()
        formula = "r_{}*{}".format(s1, a)
        kinetic_law.setMath(libsbml.parseL3Formula(formula))

        reactant = comp.createReactant()
        reactant.setSpecies(s1)
        reactant.setStoichiometry(1)    

        reactant = comp.createReactant()
        reactant.setSpecies(s2)
        reactant.setStoichiometry(1)    

        product = comp.createProduct()
        product.setSpecies(s1)
        product.setStoichiometry(1)    

        return model 
