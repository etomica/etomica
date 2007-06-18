package etomica.species;
import java.lang.reflect.Constructor;

import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;

/**
 * Species in which molecules are each made of a single spherical atom.
 * Does not permit multiatomic molecules.  The advantage of this species
 * over the multiatomic version (used with 1 atom), is that one layer of
 * the atom hierarchy is eliminated in SpeciesSpheresMono.  Each atom is
 * the direct child of the species agent (i.e., each atom is at the "molecule"
 * level in the hierarchy, without an intervening AtomGroup).
 * 
 * @author David Kofke
 */
public class SpeciesSpheresMono extends Species {

    /**
     * Constructs instance with a default element
     */
    public SpeciesSpheresMono(ISimulation sim) {
        this(sim, new ElementSimple(sim));
    }
    
    public SpeciesSpheresMono(ISimulation sim, Element element) {
        this(sim, new AtomTypeSphere(element));
    }
    
    private SpeciesSpheresMono(ISimulation sim, AtomTypeSphere atomType) {
        super(sim.isDynamic() ?
                    new AtomFactoryMonoDynamic(sim.getSpace(), atomType) :
                    new AtomFactoryMono(sim.getSpace(), atomType));
    }
    
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class,Element.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{((AtomTypeLeaf)factory.getType()).getElement()});
    }
    
    private static final long serialVersionUID = 1L;
}
