package etomica;
import etomica.action.AtomActionAdapter;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
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
public class SpeciesSpheres extends Species implements EtomicaElement {

    private double mass;
    private final AtomTypeSphere atomType;
    
    public SpeciesSpheres(Simulation sim) {
        this(sim, sim.potentialMaster.sequencerFactory(), 1);
    }
    public SpeciesSpheres(Simulation sim, AtomSequencerFactory seqFactory, int nA) {
        this(sim, seqFactory, nA, new ConformationLinear(sim.space));
    }
    public SpeciesSpheres(Simulation sim, AtomSequencerFactory seqFactory, 
            int nA, Conformation conformation) {
        this(sim, seqFactory, nA, conformation, Species.makeAgentType(sim));
    }
    private SpeciesSpheres(Simulation sim, AtomSequencerFactory seqFactory, 
            int nA, Conformation conformation, AtomTypeGroup agentType) {
        super(sim, new AtomFactoryHomo(sim.space, seqFactory, agentType,
                                nA, conformation), agentType);
        atomType = new AtomTypeSphere(factory.getType(), Default.ATOM_MASS, Default.ATOM_SIZE);
        ((AtomFactoryHomo)factory).setChildFactory(new AtomFactoryMono(sim.space, atomType, seqFactory));
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
        allAgents(new AtomActionAdapter() {public void actionPerformed(Atom a) {
            SpeciesAgent agent = (SpeciesAgent)a;
            agent.setNMolecules(agent.getNMolecules());}});
    }
    /**
     * @return the number of spheres in each molecule made by this species.
     */
    public int getAtomsPerMolecule() {
        return ((AtomFactoryHomo)factory).getAtomsPerGroup();
    }

    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
	    IntegratorHard integratorHard1 = new IntegratorHard();
//	    integratorHard1.setTimeStep(0.02);
	    SpeciesSpheres speciesSpheres1 = new SpeciesSpheres(10,6);//10 molecules, 3 atoms per molecule
	    SpeciesSpheres speciesSpheres2 = new SpeciesSpheres(3);
	//    speciesSpheres2.setColor(java.awt.Color.red);
	    final Phase phase = new Phase();
	    P2HardSphere potential = new P2HardSphere();
	    P2HardSphere potential2 = new P2HardSphere();
	    
	    //intermolecular potential
	    Potential2Group potential0 = new Potential2Group();
	    P2HardSphere p0Inter = new P2HardSphere(potential0);
	    
	    //intramolecular potential
	    Potential1Group potential3 = new Potential1Group();
	    P2Tether p2Tether = new P2Tether(potential3);
	    p2Tether.setIterator(new AtomPairIterator(Simulation.instance.space,
	            new AtomIteratorSequential(false),
	            new AtomIteratorBonds()));
	    
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
	    MeterEnergy meterEnergy = new MeterEnergy();
	    meterEnergy.setPhase(phase);
	    DisplayBox displayEnergy = new DisplayBox();
	    displayEnergy.setMeter(meterEnergy);
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
		potential.setSpecies(speciesSpheres1, speciesSpheres2);
        potential2.setSpecies(speciesSpheres2, speciesSpheres2);
        potential0.setSpecies(speciesSpheres1, speciesSpheres1); 
        potential3.setSpecies(speciesSpheres1);
	//    displayPhase1.setColorScheme(integratorHard1.new HighlightColliders());
	    Simulation.makeAndDisplayFrame(Simulation.instance);
	}//end of main
*/
}


