package etomica.species;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.space.CoordinateFactorySphere;
import etomica.units.Dimension;
import etomica.util.Default;

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

public class SpeciesSpheresMono extends Species implements EtomicaElement {

    private final AtomTypeSphere atomType;
    
    /**
     * Constructs instance with space and AtomSequencer.Factory taken from
     * given simulation, and using default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresMono(Simulation sim) {
        this(sim, sim.potentialMaster.sequencerFactory());
    }

    /**
     * Constructs instance with default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresMono(Simulation sim, AtomSequencerFactory seqFactory) {
        this(sim, seqFactory, Species.makeAgentType(sim));
    }
    
    private SpeciesSpheresMono(Simulation sim, AtomSequencerFactory seqFactory,
                                AtomTypeGroup agentType) {
        super(sim, new AtomFactoryMono(new CoordinateFactorySphere(sim), new AtomTypeSphere(agentType), seqFactory),
                agentType);
        factory.setSpecies(this);
        atomType = (AtomTypeSphere)((AtomFactoryMono)factory).getType();
        nMolecules = Default.MOLECULE_COUNT;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
        return info;
    }
              
    /**
     * The mass of each of the spheres.
     */
    public double getMass() {return atomType.getMass();}
    /**
     * Sets the mass of all spheres to the given value.
     */
    public void setMass(double m) {
        atomType.setMass(m);
    }
    /**
     * @return Dimension.MASS
     */
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    /**
     * The diameter of each of the spheres.
     */
    public double getDiameter() {return atomType.diameter(null);}
    /**
     * Sets the diameter of all spheres to the given value.
     */
    public void setDiameter(double d) {atomType.setDiameter(d);}
    /**
     * @return Dimension.LENGTH
     */
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
    
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
