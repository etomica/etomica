package etomica;
import etomica.units.Kelvin;
import java.util.Random;

/**
 * Simple Gibbs-ensemble Monte Carlo integrator.
 * Used to evaluate fluid-fluid phase coexistence.  Written to apply to only two phases.
 *
 * @author David Kofke
 */
//need to update to include setPhaseIteratorFactory
public class IntegratorGEMC extends IntegratorMC {
    
    private final Random rand = new Random();
    public Phase secondPhase;
    private MCMoveAtom atomDisplace1 = new MCMoveAtom();
    private MCMoveAtom atomDisplace2 = new MCMoveAtom();
    private MCMoveVolumeExchange volumeExchange = new MCMoveVolumeExchange(this);
    private MCMoveMoleculeExchange moleculeExchange = new MCMoveMoleculeExchange(this);
    
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
        this.add(atomDisplace1);
        this.add(atomDisplace2);
        this.add(volumeExchange);
        this.add(moleculeExchange);
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
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Simulation.setUnitSystem(new etomica.units.UnitSystem.LJ());
	    IntegratorGEMC integratorGEMC1 = new IntegratorGEMC();
	    SpeciesDisks speciesDisks1 = new SpeciesDisks();
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
	    DeviceSlider temperatureSlider = new DeviceSlider(integratorGEMC1, "temperature");
	    temperatureSlider.setUnit(new etomica.units.Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);
	    speciesDisks1.setNMolecules(60);
	    
//		Simulation.instance.elementCoordinator = new MediatorBasic();
		Simulation.instance.elementCoordinator.go(); 
		 
	    meter1.setPhase(phase1);
	    meter2.setPhase(phase2);
		integratorGEMC1.addIntervalListener(box1);
		integratorGEMC1.addIntervalListener(box2);
		integratorGEMC1.addIntervalListener(displayPhase1);
		integratorGEMC1.addIntervalListener(displayPhase2);
	    phase1.setIntegrator(integratorGEMC1);
	    phase2.setIntegrator(integratorGEMC1);
	    
	    ColorScheme color1 = new ColorScheme.Simple(java.awt.Color.blue);
	    ColorScheme color2 = new ColorScheme.Simple(java.awt.Color.red);
	    displayPhase1.setColorScheme(color1);
	    displayPhase2.setColorScheme(color2);
	    
	    controller1.add(integratorGEMC1);
	    
		Simulation.instance.setBackground(java.awt.Color.yellow);
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
}
