package etomica;
import java.awt.Color;
import etomica.units.Dimension;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each disk having the same mass and size (same type).
 * 
 * @author David Kofke
 */
public class SpeciesSpheres extends Species implements EtomicaElement {

    private double mass;
    public AtomType.Disk protoType;
    
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryHomo makeFactory(Simulation sim, int na, BondInitializer bondInit, Configuration config) {
        AtomFactoryMono f = new AtomFactoryMono(sim);
        AtomType type = new AtomType.Disk(f, Default.ATOM_MASS, Default.ATOM_COLOR, Default.ATOM_SIZE);
        f.setType(type);
        AtomFactoryHomo fm = new AtomFactoryHomo(sim,f, na, bondInit, config);
        return fm;
 //       return f;
    }
        
    public SpeciesSpheres() {
        this(Simulation.instance);
    }
    public SpeciesSpheres(int n) {
        this(Simulation.instance, n);
    }
    public SpeciesSpheres(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesSpheres(Simulation sim, int n) {
        this(sim, n, 1);
    }
    public SpeciesSpheres(int nM, int nA) {
        this(Simulation.instance, nM, nA);
    }
    public SpeciesSpheres(Simulation sim, int nM, int nA) {
        this(sim, nM, nA, new BondInitializerChain(), new ConfigurationLinear(sim.space));
    }
    public SpeciesSpheres(Simulation sim, int nM, int nA, BondInitializer bondInitializer, Configuration config) {
        super(sim, makeFactory(sim, nA, bondInitializer, config));
        protoType = (AtomType.Disk)((AtomFactoryMono)((AtomFactoryHomo)factory).childFactory()).type();
        nMolecules = nM;
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
                
    public final double getDiameter() {return protoType.diameter();}
    public void setDiameter(double d) {protoType.setDiameter(d);}
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
	    IntegratorHard integratorHard1 = new IntegratorHard();
//	    integratorHard1.setTimeStep(0.02);
	    SpeciesSpheres speciesDisks1 = new SpeciesSpheres(10,6);//10 molecules, 3 atoms per molecule
	    SpeciesSpheres speciesDisks2 = new SpeciesSpheres(3);
	    speciesDisks2.setColor(java.awt.Color.red);
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
		Simulation.instance.setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
		potential.setSpecies(speciesDisks1, speciesDisks2);
        potential2.setSpecies(speciesDisks2, speciesDisks2);
        potential0.setSpecies(speciesDisks1, speciesDisks1); 
        potential3.setSpecies(speciesDisks1);
	//    displayPhase1.setColorScheme(integratorHard1.new HighlightColliders());
	    Simulation.makeAndDisplayFrame(Simulation.instance);
	}//end of main

}


