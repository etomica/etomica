package etomica;
import etomica.units.Dimension;

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

    private double mass;
    public AtomType.Sphere protoType;
    
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryMono makeFactory(Space space, 
                                                AtomSequencer.Factory seqFactory) {
        AtomFactoryMono f = new AtomFactoryMono(space, seqFactory);
        AtomType type = new AtomType.Sphere(f, Default.ATOM_MASS, Default.ATOM_SIZE);
        f.setType(type);
        return f;
    }

    /**
     * Constructs instance with space and AtomSequencer.Factory taken from
     * given simulation, and using default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresMono(Simulation sim) {
        this(sim.space, sim.potentialMaster.sequencerFactory(), Default.MOLECULE_COUNT);
    }
        
    /**
     * Constructs instance with default number of molecules given by
     * Default.MOLECULE_COUNT.
     */
    public SpeciesSpheresMono(Space space, AtomSequencer.Factory seqFactory) {
        this(space, seqFactory, Default.MOLECULE_COUNT);
    }
    
    public SpeciesSpheresMono(Space space, AtomSequencer.Factory seqFactory, int nM) {
        super(makeFactory(space, seqFactory));
        factory.setSpecies(this);
        protoType = (AtomType.Sphere)((AtomFactoryMono)factory).type();
        mass = protoType.getMass();
        nMolecules = nM;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
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

    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
	    IntegratorHard integratorHard1 = new IntegratorHard();
//	    integratorHard1.setTimeStep(0.02);
	    SpeciesSpheres speciesSpheres1 = new SpeciesSpheres(10,3);
	    SpeciesSpheresMono speciesSpheres2 = new SpeciesSpheresMono(3);
	    speciesSpheres2.setColor(java.awt.Color.red);
	    final Phase phase = new Phase();
	    P2HardSphere potential = new P2HardSphere();
	    P2HardSphere potential2 = new P2HardSphere();
	    P2HardSphere potential0 = new P2HardSphere();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
	    MeterEnergy meterEnergy = new MeterEnergy();
	    meterEnergy.setPhase(phase);
	    DisplayBox displayEnergy = new DisplayBox();
	    displayEnergy.setMeter(meterEnergy);
		Simulation.instance.setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(phase);
    //    potentialAgent.setIterator(new AtomPairIterator(phase));
        potentialAgent.setIterator(new AtomPairIterator(phase,
                speciesSpheres1.getAgent(phase).makeLeafAtomIterator(),
                speciesSpheres2.getAgent(phase).makeLeafAtomIterator()));
                
        potentialAgent = (Potential2.Agent)potential2.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase,
                speciesSpheres2.getAgent(phase).makeLeafAtomIterator(),
                speciesSpheres2.getAgent(phase).makeLeafAtomIterator()));
                
        potentialAgent = (Potential2.Agent)potential0.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase,
                speciesSpheres1.getAgent(phase).makeLeafAtomIterator(),
                speciesSpheres1.getAgent(phase).makeLeafAtomIterator()));
	//    displayPhase1.setColorScheme(integratorHard1.new HighlightColliders());
	    Simulation.makeAndDisplayFrame(Simulation.instance);
	}//end of main
*/
}


