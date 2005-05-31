package etomica;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryTree;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeSphere;
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
public class SpeciesTree extends Species implements EtomicaElement {

    private double mass;
    private final AtomTypeSphere atomType;
    
    public SpeciesTree(Simulation sim) {
        this(sim, new int[] {1});
    }
    public SpeciesTree(Simulation sim, int[] nA) {
        this(sim, sim.potentialMaster.sequencerFactory(), nA);
    }
    public SpeciesTree(Simulation sim, AtomSequencerFactory seqFactory, int[] nA) {
        this(sim, seqFactory, nA, Species.makeAgentType(sim));
    }
    
    //TODO extend to permit specification of Conformation[], perhaps AtomSequencerFactory[]
    private SpeciesTree(Simulation sim, AtomSequencerFactory seqFactory, 
            int[] nA, AtomTypeGroup agentType) {
        super(sim, new AtomFactoryTree(sim.space, seqFactory, agentType, nA), agentType);
        atomType = new AtomTypeSphere(((AtomFactoryTree)factory).getLeafType(), Default.ATOM_MASS, Default.ATOM_SIZE);
        ((AtomFactoryTree)factory).setLeafFactory(new AtomFactoryMono(sim.space, atomType, seqFactory));
        factory.setSpecies(this);
        
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
      
}


