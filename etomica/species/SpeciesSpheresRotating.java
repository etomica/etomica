package etomica.species;

import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.chem.elements.ElementSimple;
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
    
    public SpeciesSpheresRotating(Simulation sim) {
        this(sim, Species.makeAgentType(sim));
    }

    private SpeciesSpheresRotating(Simulation sim, AtomTypeGroup agentType) {
        super(new AtomFactoryMono(new CoordinateFactoryAngular(sim), 
                new AtomTypeOrientedSphere(new ElementSimple(sim),sim.getDefaults().atomSize)), agentType);
        // factory.getType is the AtomTypeOrientedSphere instance we just passed to the AtomFactoryMono
        // we need to finish setting it up by setting its parent
        factory.getType().setParentType(agentType);
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

    private static final long serialVersionUID = 1L;
}
