package etomica;
import etomica.units.Kelvin;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator.
 * Used to evaluate fluid-fluid phase coexistence.  Written to apply to only two phases.
 *
 * @author David Kofke
 */
//need to update to include setPhaseIteratorFactory
public class IntegratorGEMC extends IntegratorMC implements EtomicaElement {
    
    public String version() {return "IntegratorGEMC:01.04.17"+super.version();}
    public Phase secondPhase;
    private final MCMoveAtom atomDisplace1 = new MCMoveAtom(this);
    private final MCMoveAtom atomDisplace2 = new MCMoveAtom(this);
    private final MCMoveVolumeExchange volumeExchange = new MCMoveVolumeExchange(this);
    private final MCMoveMoleculeExchange moleculeExchange = new MCMoveMoleculeExchange(this);
    
    public IntegratorGEMC() {
        this(Simulation.instance);
    }
    public IntegratorGEMC(Simulation sim) {
        super(sim);
        phaseCountMax = 2;
//        super.phaseCountMax = 2;
        phase = new Phase[phaseCountMax];
        atomDisplace1.setAdjustInterval(100);
        atomDisplace2.setAdjustInterval(100);
    }
  
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Gibbs-ensemble Monte Carlo simulation of phase coexistence");
        return info;
    }

  public boolean addPhase(Phase p) {
    for(int i=0; i<phaseCount; i++) {if(phase[i]==p) return false;}  //check that phase is not already registered
    if(!this.wantsPhase()) {return false;}  //if another phase not wanted, return false
    phase[phaseCount] = p;
    phaseCount++;
    if(phaseCount == 1) {
        firstPhase = phase[0];
        atomDisplace1.setPhase(phase[0]);
    }
    if(phaseCount == 2) {
        secondPhase = phase[1];
        volumeExchange.setPhase(phase);
        moleculeExchange.setPhase(phase);
        makeIterators(p.iteratorFactory());
	    p.iteratorFactoryMonitor.addObserver(iteratorFactoryObserver());
        atomDisplace2.setPhase(phase[1]);
    }
    return true;
  }
  
  public Phase[] getPhases() {return phase;}
  public void setPhases(Phase[] phases) {
    phase = phases;
  }
  
  /**
   * Returns the object that performs the volume-exchange move in the GE simulation.
   * Having handle to this object is needed to adjust trial frequency and view acceptance rate.
   */
  public MCMoveVolumeExchange getMCMoveVolumeExchange() {return volumeExchange;}
  /**
   * Returns the object that performs the molecule-exchange move in the GE simulation.
   * Having handle to this object is needed to adjust trial frequency and view acceptance rate.
   */
  public MCMoveMoleculeExchange getMCMoveMoleculeExchange() {return moleculeExchange;}
  /**
   * Returns the object that performs the atom-displacement move in the GE simulation.
   * Having handle to this object is needed to adjust trial frequency and view acceptance rate.
   * @param i indicates request for move for 0th phase (i = 0) or first phase (i = 1)
   */
  public MCMoveAtom getMCMoveAtom(int i) {return (i==0) ? atomDisplace1 : atomDisplace2;}
        
    public static void main(String[] args) {
        Simulation.setUnitSystem(new etomica.units.UnitSystem.LJ());
	    IntegratorGEMC integratorGEMC1 = new IntegratorGEMC();
//	    SpeciesSpheres speciesSphere1 = new SpeciesSphere();
        Species speciesSpheres1 = new SpeciesSpheresMono();
	    Phase phase1 = new Phase();
	    Phase phase2 = new Phase();
//	    P2SquareWell P2SquareWell1 = new P2SquareWell();
	    P2LennardJones P2LennardJones1 = new P2LennardJones();
	    Controller controller1 = new Controller();
	    //Configuration displays for each phase
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    DisplayPhase displayPhase2 = new DisplayPhase();
	    displayPhase1.setPhase(phase1);
	    displayPhase2.setPhase(phase2);
	    //Meters and displays for density in each phase
	    MeterDensity meter1 = new MeterDensity();
	    MeterDensity meter2 = new MeterDensity();
	    DisplayBox box1 = new DisplayBox();
	    DisplayBox box2 = new DisplayBox();
	    box1.setMeter(meter1);
	    box2.setMeter(meter2);
	    //Slider to adjust temperature
	    Modulator modT = new Modulator(integratorGEMC1, "temperature");
	    DeviceSlider temperatureSlider = new DeviceSlider(modT);
	    temperatureSlider.setUnit(new etomica.units.Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);
	    speciesSpheres1.setNMolecules(60);
	    DisplayBox boxT = new DisplayBox(modT);
	    	    
		Simulation.instance.elementCoordinator.go(); 
		 
        P2LennardJones1.setIterator(
            new AtomPairIterator(Simulation.instance.space,
                                 new AtomIteratorSequential(true),
                                 new AtomIteratorSequential(true)));
 /*       potential.set(species.getAgent(phase));
		
        Potential2.Agent potentialAgent = (Potential2.Agent)P2LennardJones1.getAgent(phase1);
        potentialAgent.setIterator(new AtomPairIterator(phase1));
        potentialAgent = (Potential2.Agent)P2LennardJones1.getAgent(phase2);
        potentialAgent.setIterator(new AtomPairIterator(phase2));
*/
	    meter1.setPhase(phase1);
	    meter2.setPhase(phase2);
	    phase1.setIntegrator(integratorGEMC1);
	    phase2.setIntegrator(integratorGEMC1);
	    
	    ColorScheme color1 = new ColorScheme.Simple(java.awt.Color.blue);
	    ColorScheme color2 = new ColorScheme.Simple(java.awt.Color.red);
	    displayPhase1.setColorScheme(color1);
	    displayPhase2.setColorScheme(color2);
	    
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
		
		Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
    
}
