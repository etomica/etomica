package etomica.species;

import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactoryAngular;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 * 
 */

public class SpeciesSpheresRotating extends Species implements EtomicaElement {
    
    /**
     * Constructs instance with default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresRotating(Simulation sim) {
        this(sim, Species.makeAgentType(sim));
    }
    private SpeciesSpheresRotating(Simulation sim, AtomTypeGroup agentType) {
        super(sim, new AtomFactoryMono(new CoordinateFactoryAngular(sim), 
                new AtomTypeOrientedSphere(agentType,sim.getDefaults().atomMass,sim.getDefaults().atomSize)), agentType);
    }
            
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecules formed from spheres with an attached rotatable direction");
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
        return new SpeciesSignature(getName(),constructor,new Object[]{});
    }
}
