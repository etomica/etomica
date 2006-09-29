package etomica.species;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactoryHetero;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactorySphere;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */

/* History
 * 08/12/03 (DAK) use sim instead of space in AtomFactoryHomo constructor
 */
public class SpeciesSpheresHetero extends Species implements EtomicaElement {

    /**
     * Constructs instance with 0 components and total number of children 
     * equal to 1.  The actual number of components must be set in the factory
     * (AtomFactoryHetero) before use.  The actual number of desired children 
     * can also be set in the factory.
     */
    public SpeciesSpheresHetero(Simulation sim) {
        this(sim,0);
    }
    
    /**
     * Constructs instance with the given number of components and 
     * total number of children equal to 1.  The actual number of desired 
     * desired children can be set in the factory (AtomFactoryHetero) after
     * construction.
     */
    public SpeciesSpheresHetero(Simulation sim, int nComponents) {
        this(sim,nComponents,1);
    }
    
    /**
     * Constructs instance with the given number of child types (components), 
     * total number of children and AtomSequencer.Factory.  The number of 
     * molecules taken from the simulation.
     */
    public SpeciesSpheresHetero(Simulation sim, int nComponents, int nA) {
        this(sim, nComponents, nA, Species.makeAgentType(sim));
    }
    
    private SpeciesSpheresHetero(Simulation sim, int nComponents, int nA, 
            AtomTypeGroup agentType) {
        super(sim, new AtomFactoryHetero(sim, agentType), agentType);
        if (nA < 1) {
            throw new IllegalArgumentException("You must have at least one child atom");
        }
        if (nComponents > 0) {
            AtomFactoryMono[] childFactories = new AtomFactoryMono[nComponents];
            for (int i=0; i<nComponents; i++) {
                AtomTypeSphere atomType = new AtomTypeSphere(sim.getDefaults().atomMass, sim.getDefaults().atomSize);
                atomType.setParentType((AtomTypeGroup)factory.getType());
                childFactories[i] = new AtomFactoryMono(new CoordinateFactorySphere(sim), atomType);
            }
            ((AtomFactoryHetero)factory).setChildFactory(childFactories);
        }
        ((AtomFactoryHetero)factory).setTotalChildren(nA);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
        return info;
    }

    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{Simulation.class,int.class,int.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{new Integer(((AtomFactoryHetero)factory).getNumChildAtoms()),
            new Integer(((AtomFactoryHetero)factory).getChildFactory().length)});
    }
    
}
