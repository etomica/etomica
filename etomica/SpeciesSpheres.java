package etomica;
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
    public AtomType.Sphere protoType;
    
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryHomo makeFactory(Space space, AtomSequencer.Factory seqFactory, int na, BondInitializer bondInit, Configuration config) {
        AtomFactoryMono f = new AtomFactoryMono(space, seqFactory);//would like to pass this species
        AtomType type = new AtomType.Sphere(f, Default.ATOM_MASS, Default.ATOM_SIZE);
        f.setType(type);
        AtomFactoryHomo fm = new AtomFactoryHomo(space, seqFactory, 
                                f, na, bondInit, config);
        return fm;
 //       return f;
    }
    
    public SpeciesSpheres(Simulation sim) {
        this(sim.space, sim.potentialMaster.sequencerFactory(), Default.MOLECULE_COUNT, 1);
    }
    public SpeciesSpheres(Space space, AtomSequencer.Factory seqFactory, int nM, int nA) {
        this(space, seqFactory, nM, nA, new BondInitializerChain(), new ConfigurationLinear(space));
    }
    public SpeciesSpheres(Space space, AtomSequencer.Factory seqFactory, 
                int nM, int nA, BondInitializer bondInitializer, Configuration config) {
        super(makeFactory(space, seqFactory, nA, bondInitializer, config));
        factory.setSpecies(this);
        protoType = (AtomType.Sphere)((AtomFactoryMono)((AtomFactoryHomo)factory).childFactory()).type();
        nMolecules = nM;
        mass = protoType.getMass();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
        return info;
    }
              
    // Exposed Properties
    public final double getMass() {return mass;}
    public final void setMass(double m) {
        mass = m;
        allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.coord.setMass(mass);}});
    }
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    public final double getDiameter() {return protoType.diameter(null);}
    public void setDiameter(double d) {protoType.setDiameter(d);}
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
                    
    public void setAtomsPerMolecule(final int n) {
        ((AtomFactoryHomo)factory).setAtomsPerGroup(n);
        allAgents(new AtomAction() {public void actionPerformed(Atom a) {
            SpeciesAgent agent = (SpeciesAgent)a;
            agent.setNMolecules(agent.getNMolecules(),true);}});
    }
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


