package etomica.species;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactorySphere;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */

public class SpeciesSpheres extends Species implements EtomicaElement {

    public SpeciesSpheres(Simulation sim) {
        this(sim, 1);
    }
    public SpeciesSpheres(Simulation sim, int nA) {
        this(sim, nA, new ElementSimple(sim));
    }
    
    public SpeciesSpheres(Simulation sim, int nA, Element leafElement) {
        this(sim, nA, leafElement, new ConformationLinear(sim));
    }
    
    public SpeciesSpheres(Simulation sim, int nA, Element leafElement, Conformation conformation) {
        this(sim, nA, leafElement, conformation, Species.makeAgentType(sim));
    }
    
    private SpeciesSpheres(Simulation sim, int nA, Element leafElement, Conformation conformation, 
            AtomTypeGroup agentType) {
        super(sim, new AtomFactoryHomo(sim.space, agentType,
                                nA, conformation), agentType);
        AtomTypeSphere atomType = new AtomTypeSphere(sim, leafElement);
        atomType.setParentType((AtomTypeGroup)factory.getType());
        ((AtomFactoryHomo)factory).setChildFactory(
                new AtomFactoryMono(new CoordinateFactorySphere(sim), atomType));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
        return info;
    }

    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{Simulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{new Integer(((AtomFactoryHomo)factory).getNumChildAtoms())});
    }
    
}
