package etomica;
import java.awt.Color;
import etomica.units.Dimension;

/**
 * Species in which molecules are made of arbitrary number of disks,
 * with each disk having the same mass and size (same type).
 * 
 * @author David Kofke
 */
public class SpeciesDisks extends Species implements EtomicaElement {

    private double mass;
    public AtomType.Disk protoType;
    
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryHomo makeFactory(Simulation sim, int na) {
        AtomFactoryMono f = new AtomFactoryMono(sim);
        AtomType type = new AtomType.Disk(f, Default.ATOM_MASS, Default.ATOM_COLOR, Default.ATOM_SIZE);
        f.setType(type);
        AtomFactoryHomo fm = new AtomFactoryHomo(sim,f, na);
        return fm;
 //       return f;
    }
        
    public SpeciesDisks() {
        this(Simulation.instance);
    }
    public SpeciesDisks(int n) {
        this(Simulation.instance, n);
    }
    public SpeciesDisks(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesDisks(Simulation sim, int n) {
        this(sim, n, 1);
    }
    public SpeciesDisks(int nM, int nA) {
        this(Simulation.instance, nM, nA);
    }
    public SpeciesDisks(Simulation sim, int nM, int nA) {
        super(sim, makeFactory(sim, nA));
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
	    SpeciesDisks speciesDisks1 = new SpeciesDisks(10,3);
	    SpeciesDisks speciesDisks2 = new SpeciesDisks(3);
	    speciesDisks2.setColor(java.awt.Color.red);
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
		                                    
		potential.setSpecies(speciesDisks1, speciesDisks2);
        potential2.setSpecies(speciesDisks2, speciesDisks2);
        potential0.setSpecies(speciesDisks1, speciesDisks1);        
	//    displayPhase1.setColorScheme(integratorHard1.new HighlightColliders());
	    Simulation.makeAndDisplayFrame(Simulation.instance);
	}//end of main

}


