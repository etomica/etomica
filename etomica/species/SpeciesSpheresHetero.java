package etomica.species;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryHetero;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
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
        this(sim,makeElements(sim,nComponents));
    }
    
    private static Element[] makeElements(Simulation sim, int nComponents) {
        ElementSimple[] elements = new ElementSimple[nComponents];
        for (int i=0; i<elements.length; i++) {
            elements[i] = new ElementSimple(sim);
        }
        return elements;
    }
    
    /**
     * Constructs instance with the given leaf elements and 
     * total number of children equal to 1.  The actual number of desired 
     * desired children can be set in the factory (AtomFactoryHetero) after
     * construction.
     */
    public SpeciesSpheresHetero(Simulation sim, Element[] leafElements) {
        this(sim, leafElements, 1);
    }
    
    /**
     * Constructs instance with the given number of child types (components), 
     * total number of children and AtomSequencer.Factory.  The number of 
     * molecules taken from the simulation.
     */
    public SpeciesSpheresHetero(Simulation sim, Element[] leafElements, int nA) {
        super(new AtomFactoryHetero(sim));
        if (leafElements.length > 0) {
            AtomFactoryMono[] childFactories = new AtomFactoryMono[leafElements.length];
            for (int i=0; i<leafElements.length; i++) {
                AtomTypeSphere atomType = new AtomTypeSphere(sim, leafElements[i]);
                childFactories[i] = sim.isDynamic() ?
                            new AtomFactoryMonoDynamic(sim.getSpace(), atomType) :
                            new AtomFactoryMono(sim.getSpace(), atomType);
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
            constructor = this.getClass().getConstructor(new Class[]{Simulation.class,Element[].class,int.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        AtomFactory[] childFactories = ((AtomFactoryHetero)factory).getChildFactory();
        Element[] leafElements = new Element[childFactories.length];
        for (int i=0; i<childFactories.length; i++) {
            leafElements[i] = ((AtomTypeLeaf)childFactories[i].getType()).getElement();
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{leafElements, 
            new Integer(((AtomFactoryHetero)factory).getNumChildAtoms())});
    }
    
    private static final long serialVersionUID = 1L;
}
