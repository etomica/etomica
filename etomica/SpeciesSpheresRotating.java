package etomica;

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
    private static AtomFactoryMono makeFactory(Space space, AtomSequencerFactory seqFactory) {
        AtomFactoryMono f = new AtomFactoryMono(space, seqFactory);
        AtomType type = new AtomTypeOrientedSphere(f, Default.ATOM_MASS, Default.ATOM_SIZE);
        f.setType(type);
        return f;
    }
    /**
     * Constructs instance with space and AtomSequencer.Factory taken from
     * given simulation, and using default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresRotating(Simulation sim) {
        this(sim.space, sim.potentialMaster.sequencerFactory(), Default.MOLECULE_COUNT);
    }
        
    /**
     * Constructs instance with default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresRotating(Space space, AtomSequencerFactory seqFactory) {
        this(space, seqFactory, Default.MOLECULE_COUNT);
    }
    
    public SpeciesSpheresRotating(Space space, AtomSequencerFactory seqFactory, int nM) {
        super(makeFactory(space, seqFactory));
        factory.setSpecies(this);
        protoType = (AtomTypeOrientedSphere)((AtomFactoryMono)factory).type();
        mass = protoType.getMass();
        nMolecules = nM;
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


