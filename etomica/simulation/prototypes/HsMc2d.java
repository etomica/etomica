package etomica.simulation.prototypes;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationLattice;
import etomica.data.DataSourceCountSteps;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.box.Box;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HsMc2d extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species, species2;
    public Box box;
    public P2HardSphere potential;
    public Controller controller;
    public DataSourceCountSteps meterCycles;

    public HsMc2d() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(space);
	    integrator = new IntegratorMC(this, potentialMaster);
	    mcMoveAtom = new MCMoveAtom(this, potentialMaster);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species2 = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        getSpeciesManager().addSpecies(species2);
        box = new Box(this);
        addBox(box);
        box.getAgent(species).setNMolecules(20);
        box.getAgent(species2).setNMolecules(20);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(box);
	    potential = new P2HardSphere(space);
        potentialMaster.addPotential(potential, new Species[] {species, species});
        potentialMaster.addPotential(potential, new Species[] {species, species2});
        potentialMaster.addPotential(potential, new Species[] {species2, species2});
	    meterCycles = new DataSourceCountSteps(integrator);

        integrator.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.addIntervalAction(new BoxImposePbc(box));

//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);
	    
//		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(box);
//        colorSchemeCell.setLattice(lattice);
    }
    
}