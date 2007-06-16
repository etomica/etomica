package etomica.species;

import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoAngular;
import etomica.atom.AtomFactoryMonoAngularDynamic;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 * 
 */
public class SpeciesSpheresRotating extends Species implements EtomicaElement {
    
    public SpeciesSpheresRotating(Simulation sim) {
        super(makeAtomFactory(sim, new AtomTypeOrientedSphere(new ElementSimple(sim),1.0)));
    }
    
    private static AtomFactoryMono makeAtomFactory(Simulation sim, AtomTypeOrientedSphere atomType) {
        if (sim.isDynamic()) {
            new AtomFactoryMonoAngularDynamic(sim.getSpace(), atomType);
        }
        return new AtomFactoryMonoAngular(sim.getSpace(), atomType);
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
