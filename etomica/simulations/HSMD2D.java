package etomica.simulations;
import etomica.Controller;
import etomica.Default;
import etomica.IntegratorHard;
import etomica.IteratorDirective;
import etomica.Mediator;
import etomica.P2HardSphere;
import etomica.Phase;
import etomica.Potential2;
import etomica.Space2D;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.nbr.PotentialCalculationNbrSetup;
import etomica.nbr.PotentialMasterNbr;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D extends SimulationGraphic {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public Potential2 potential;
//    public Controller controller;
    public DisplayPhase display;

    public HSMD2D() {
    	this(new Space2D());
    }
    
    public HSMD2D(Space2D space) {
        super(space, new PotentialMasterNbr(space));
        Default.makeLJDefaults();
  //can't use cell list until integrator is updated for it      setIteratorFactory(new IteratorFactoryCell(this));
//        Default.BOX_SIZE = 30.0;
        Default.ATOM_SIZE = 0.4;
        integrator = new IntegratorHard(this);
        integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
        integrator.setInterval(1);
        integrator.setTimeStep(0.01);
        species = new SpeciesSpheresMono(this);
	    species2 = new SpeciesSpheresMono(this);
	    species.setNMolecules(500);
	    species2.setNMolecules(5);
	    phase = new Phase(this);
	    potential = new P2HardSphere(space);
//	    potential = new P2HardSphere();
	    this.potentialMaster.setSpecies(potential,new Species[]{species,species});
	    this.potentialMaster.setSpecies(potential,new Species[]{species2,species2});
	    this.potentialMaster.setSpecies(potential,new Species[]{species,species2});
//	    controller = new Controller(this);
	    new DeviceTrioControllerButton(this, getController());
	    display = new DisplayPhase(this);
//	    DisplayTimer timer = new DisplayTimer(integrator);
//	    timer.setUpdateInterval(10);
	    ColorSchemeByType.setColor(species, java.awt.Color.red);
	    ColorSchemeByType.setColor(species2, java.awt.Color.blue);
		panel().setBackground(java.awt.Color.yellow);
		elementCoordinator.go();
		this.potentialMaster.calculate(phase,new IteratorDirective(),new PotentialCalculationNbrSetup());
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
		SimulationGraphic.makeAndDisplayFrame(sim);
	//	sim.controller.start();
    }//end of main
    
}