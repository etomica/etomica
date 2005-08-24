package etomica;
import java.lang.reflect.Constructor;

import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.space.CoordinateFactorySphere;
import etomica.units.Dimension;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */

/* History
 * 08/12/03 (DAK) use sim instead of space in AtomFactoryHomo constructor
 */
public class SpeciesSpheres extends Species implements EtomicaElement {

    private double mass;
    private final AtomTypeSphere atomType;
    
    public SpeciesSpheres(Simulation sim) {
        this(sim, 1);
    }
    public SpeciesSpheres(Simulation sim, int nA) {
        this(sim, sim.potentialMaster.sequencerFactory(), nA, new ConformationLinear(sim.space));
    }
    public SpeciesSpheres(Simulation sim, AtomSequencerFactory seqFactory, 
            int nA, Conformation conformation) {
        this(sim, seqFactory, nA, conformation, Species.makeAgentType(sim));
    }
    private SpeciesSpheres(Simulation sim, AtomSequencerFactory seqFactory, 
            int nA, Conformation conformation, AtomTypeGroup agentType) {
        super(sim, new AtomFactoryHomo(sim.space, seqFactory, agentType,
                                nA, conformation), agentType);
        atomType = new AtomTypeSphere((AtomTypeGroup)factory.getType(), Default.ATOM_MASS, Default.ATOM_SIZE);
        ((AtomFactoryHomo)factory).setChildFactory(
                new AtomFactoryMono(new CoordinateFactorySphere(sim), atomType, seqFactory));
        factory.setSpecies(this);
//        ((AtomFactoryHomo)factory).getType().setChildTypes(new AtomType[]{atomType});
        
        nMolecules = Default.MOLECULE_COUNT;
        mass = atomType.getMass();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
        return info;
    }
              
    /**
     * The mass of each of the spheres that form a molecule.
     */
    public final double getMass() {return mass;}
    /**
     * Sets the mass of all spheres in each molecule to the given value.
     */
    public final void setMass(double m) {
        mass = m;
        atomType.setMass(m);
    }
    /**
     * @return Dimension.MASS
     */
    public Dimension getMassDimension() {return Dimension.MASS;}
      
    /**
     * The diameter of each of the spheres that form a molecule.
     */
    public final double getDiameter() {return atomType.diameter(null);}
    /**
     * Sets the diameter of all spheres in each molecule to the given value.
     */
    public void setDiameter(double d) {atomType.setDiameter(d);}
    /**
     * @return Dimension.LENGTH
     */
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
    
    /**
     * Sets the number of spheres in each molecule.  Causes reconstruction
     * of all molecules of this species in all phases.  No action is performed
     * if the given value equals the current value.
     * @param n new number of atoms in each molecule.
     */
    public void setAtomsPerMolecule(final int n) {
        if(n == getAtomsPerMolecule()) return;
        ((AtomFactoryHomo)factory).setAtomsPerGroup(n);
    }
    /**
     * @return the number of spheres in each molecule made by this species.
     */
    public int getAtomsPerMolecule() {
        return ((AtomFactoryHomo)factory).getAtomsPerGroup();
    }

    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{Simulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{new Integer(getAtomsPerMolecule())});
    }
    
}
