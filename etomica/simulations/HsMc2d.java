package etomica.simulations;
import etomica.Controller;
import etomica.DataSourceCountSteps;
import etomica.IntegratorMC;
import etomica.MCMoveAtom;
import etomica.P2HardSphere;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space2D;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HsMc2d extends Simulation {
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public P2HardSphere potential;
    public Controller controller;
    public DataSourceCountSteps meterCycles;

    public HsMc2d() {
        super(new Space2D());
	    integrator = new IntegratorMC(potentialMaster);
	    mcMoveAtom = new MCMoveAtom(potentialMaster);
        integrator.addMCMove(mcMoveAtom);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species2 = new SpeciesSpheresMono(this);
	    phase = new Phase(space);
	    potential = new P2HardSphere();
        potentialMaster.setSpecies(potential, new Species[] {species, species});
        potentialMaster.setSpecies(potential, new Species[] {species, species2});
        potentialMaster.setSpecies(potential, new Species[] {species2, species2});
	    meterCycles = new DataSourceCountSteps();

        phase.speciesMaster.addSpecies(species);
        phase.speciesMaster.addSpecies(species2);
        integrator.addPhase(phase);
        integrator.addIntervalListener(new PhaseImposePbc(phase));

//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);
	    
//		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(phase);
//        colorSchemeCell.setLattice(lattice);
    }
    
}