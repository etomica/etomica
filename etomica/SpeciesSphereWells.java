package etomica;
import etomica.action.AtomActionAdapter;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomTypeWell;
import etomica.units.Dimension;

 // Dumb, verbatim copy of SpeciesSpheres, with AtomType set to Well instead of Sphere

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */
 
 /* History
  * 10/18/02 Modified makeFactory method to construct AtomFactoryHomo using sim.
  * arguments, instead of sim. 
  * 08/12/03 (DAK) use sim instead of sim.space AtomFactoryMono constructor
  */
public class SpeciesSphereWells extends Species implements EtomicaElement {

    private double mass;
    public AtomTypeSphere protoType;
    
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryHomo makeFactory(Space space, AtomSequencerFactory seqFactory,
                                int na, Configuration config) {
        AtomFactoryMono f = new AtomFactoryMono(space, seqFactory);
        AtomType type = new AtomTypeWell(f, Default.ATOM_MASS, Default.ATOM_SIZE, 1.5);
        f.setType(type);
        AtomFactoryHomo fm = new AtomFactoryHomo(space, seqFactory, f, na, config);
        return fm;
 //       return f;
    }
   
    public SpeciesSphereWells(Simulation sim) {
        this(sim.space, sim.potentialMaster.sequencerFactory(), Default.MOLECULE_COUNT, 1);
    }
    public SpeciesSphereWells(Space space, AtomSequencerFactory seqFactory, int nM, int nA) {
        this(space, seqFactory, nM, nA, new ConfigurationLinear(space));
    }
    public SpeciesSphereWells(Space space, AtomSequencerFactory seqFactory, 
                int nM, int nA, Configuration config) {
        super(makeFactory(space, seqFactory, nA, config));
        factory.setSpecies(this);
        protoType = (AtomTypeSphere)((AtomFactoryMono)((AtomFactoryHomo)factory).childFactory()).type();
        nMolecules = nM;
        mass = protoType.getMass();
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
            agent.setNMolecules(agent.getNMolecules(),true);}});
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
	    speciesSpheres2.setColor(java.awt.Color.red);
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


