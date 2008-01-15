package etomica.species;
import java.lang.reflect.Constructor;

import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.simulation.ISimulation;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */

public class SpeciesSpheres extends Species {

    public SpeciesSpheres(ISimulation sim) {
        this(sim, 1);
    }
    public SpeciesSpheres(ISimulation sim, int nA) {
        this(sim, nA, new ElementSimple(sim));
    }
    
    public SpeciesSpheres(ISimulation sim, int nA, Element leafElement) {
        this(sim, nA, leafElement, new ConformationLinear(sim));
    }
    
    public SpeciesSpheres(ISimulation sim, int nA, Element leafElement, Conformation conformation) {
        super();
        setMoleculeFactory(new AtomFactoryHomo(this, sim.getSpace(), nA, conformation));
        AtomTypeSphere atomType = new AtomTypeSphere(leafElement);
        ((AtomFactoryHomo)factory).setChildFactory(sim.isDynamic() ?
                            new AtomFactoryMonoDynamic(sim.getSpace(), atomType) :
                            new AtomFactoryMono(sim.getSpace(), atomType));
    }
    
    public AtomTypeLeaf getLeafType() {
        return (AtomTypeLeaf)getMoleculeType().getChildTypes()[0];
    }
    
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(constructor,new Object[]{new Integer(factory.getNumChildAtoms())});
    }
    
    private static final long serialVersionUID = 1L;
}
