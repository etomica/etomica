package etomica;

import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.units.Dimension;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 * 
 */

public class SpeciesSpheresRotating extends Species implements EtomicaElement {
    
    public double mass;
    
    public AtomTypeOrientedSphere protoType;
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryMono makeFactory(Space space, AtomSequencerFactory seqFactory,
                                AtomIndexManager indexManager) {
        AtomType type = new AtomTypeOrientedSphere(indexManager, Default.ATOM_MASS, Default.ATOM_SIZE);
        return new AtomFactoryMono(space, type, seqFactory);
    }
    /**
     * Constructs instance with space and AtomSequencer.Factory taken from
     * given simulation, and using default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresRotating(Simulation sim) {
        this(sim, sim.potentialMaster.sequencerFactory());
    }
        
    /**
     * Constructs instance with default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresRotating(Simulation sim, AtomSequencerFactory seqFactory) {
        this(sim, seqFactory, Species.makeAgentType(sim));
    }
    private SpeciesSpheresRotating(Simulation sim, AtomSequencerFactory seqFactory,
                                   AtomType agentType) {
        super(sim, makeFactory(sim.space, seqFactory, 
                agentType.getIndexManager().makeChildManager()), agentType);
        factory.setSpecies(this);
        protoType = (AtomTypeOrientedSphere)((AtomFactoryMono)factory).getType();
        mass = protoType.getMass();
        nMolecules = Default.MOLECULE_COUNT;
    }
            
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecules formed from spheres with an attached rotatable direction");
        return info;
    }
              
    /**
     * The mass of each of the spheres.
     */
    public final double getMass() {return mass;}
    /**
     * Sets the mass of all spheres to the given value.
     */
    public final void setMass(double m) {
        mass = m;
        protoType.setMass(m);
    }
    /**
     * @return Dimension.MASS
     */
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    /**
     * The diameter of each of the spheres.
     */
    public final double getDiameter() {return protoType.diameter(null);}
    /**
     * Sets the diameter of all spheres to the given value.
     */
    public void setDiameter(double d) {protoType.setDiameter(d);}
    /**
     * @return Dimension.LENGTH
     */
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
 
}


