package etomica.simulations;
import etomica.*;
import etomica.graphics.*;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HsMc2d extends SimulationGraphic {
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2HardSphere potential;
    public Controller controller;
    public DisplayPhase display;
    public MeterCycles meterCycles;
    public DisplayBox displayCycles;

    public HsMc2d() {
        super(new Space2D());
 //       super(new etomica.space.continuum.Space(2));
 //       setIteratorFactory(new IteratorFactoryCell(this));
	    integrator = new IntegratorMC(this);
	    mcMoveAtom = new MCMoveAtom(integrator);
	    species = new SpeciesSpheresMono(this);
	    phase = new Phase(this);
	    potential = new P2HardSphere();
	    potential.setSpecies(new Species[] {species});
	    controller = new Controller(this);
	    display = new DisplayPhase(this);
	    meterCycles = new MeterCycles(this);
	    displayCycles = new DisplayBox(this,meterCycles);
		panel().setBackground(java.awt.Color.yellow);
	    
//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);
	    
		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(phase);
//        colorSchemeCell.setLattice(lattice);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HsMc2d sim = new HsMc2d();
		sim.elementCoordinator.go(); 
		SimulationGraphic.makeAndDisplayFrame(sim);
	//	sim.controller.start();
    }//end of main
    
}