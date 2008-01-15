package etomica.species;

import java.lang.reflect.Constructor;

import etomica.EtomicaInfo;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoAngular;
import etomica.atom.AtomFactoryMonoAngularDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.simulation.ISimulation;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 * 
 */
public class SpeciesSpheresRotating extends Species {
    
    public SpeciesSpheresRotating(ISimulation sim) {
        super();
        setMoleculeFactory(new AtomFactoryHomo(this, sim.getSpace(), 1, new ConformationLinear(sim.getSpace(), 1)));
        AtomTypeOrientedSphere atomType = new AtomTypeOrientedSphere(new ElementSimple(sim), 1.0);
        AtomFactory childFactory = null;
        if (sim.isDynamic()) {
            childFactory = new AtomFactoryMonoAngularDynamic(sim.getSpace(), atomType);
        }
        else {
            childFactory = new AtomFactoryMonoAngular(sim.getSpace(), atomType);
        }
        ((AtomFactoryHomo)factory).setChildFactory(childFactory);
    }

    public AtomTypeLeaf getLeafType() {
        return (AtomTypeLeaf)getMoleculeType().getChildTypes()[0];
    }
    
    private static AtomFactoryMono makeAtomFactory(ISimulation sim, AtomTypeOrientedSphere atomType) {
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
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(constructor,new Object[]{});
    }

    private static final long serialVersionUID = 1L;
}
