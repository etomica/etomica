package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationSequential;
import etomica.data.DataSourceCountSteps;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Debug;

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
        super(Space2D.getInstance());
	    integrator = new IntegratorMC(this);
	    mcMoveAtom = new MCMoveAtom(this);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this, integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species2 = new SpeciesSpheresMono(this);
        species.setNMolecules(20);
        species2.setNMolecules(20);
        phase = new Phase(this);
        phase.makeMolecules();
        Debug.setAtoms(phase);
        new ConfigurationSequential(space).initializeCoordinates(phase);
	    potential = new P2HardSphere(this);
        potentialMaster.setSpecies(potential, new Species[] {species, species});
        potentialMaster.setSpecies(potential, new Species[] {species, species2});
        potentialMaster.setSpecies(potential, new Species[] {species2, species2});
	    meterCycles = new DataSourceCountSteps();

        integrator.addPhase(phase);
        integrator.addMCMove(mcMoveAtom);
        integrator.addListener(new PhaseImposePbc(phase));

//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);
	    
//		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(phase);
//        colorSchemeCell.setLattice(lattice);
    }
    
}