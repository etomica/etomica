package etomica.simulation.prototypes;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.config.ConfigurationSequential;
import etomica.integrator.IntegratorHard;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public Potential2 potential;
    public Potential2 potential2;
    public Potential2 potential22;

    public HSMD2D() {
        this(new Default());
    }
    
    public HSMD2D(Default defaults) {
    	this(Space2D.getInstance(), defaults);
    }
    
    private HSMD2D(Space2D space, Default defaults) {
        super(space, true, new PotentialMasterList(space), Default.BIT_LENGTH, defaults);
//        super(space, new PotentialMaster(space));//,IteratorFactoryCell.instance));
        defaults.makeLJDefaults();
        defaults.atomSize = 0.38;

        double neighborRangeFac = 1.6;
        ((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*defaults.atomSize);

        integrator = new IntegratorHard(this);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager();
        nbrManager.setRange(defaults.atomSize*1.6);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
	    species2 = new SpeciesSpheresMono(this);
        species.setNMolecules(512);
        species2.setNMolecules(5);
        potential = new P2HardSphere(this);
        potential2 = new P2HardSphere(this);
        potential22 = new P2HardSphere(this);
//	    this.potentialMaster.setSpecies(potential,new Species[]{species,species});
	    //this.potentialMaster.setSpecies(potential,new Species[]{species2,species2});
	    //this.potentialMaster.setSpecies(potential,new Species[]{species,species2});
        
        NeighborCriterion criterion = new CriterionSimple(this,potential.getRange(),neighborRangeFac*potential.getRange());
        potential.setCriterion(criterion);
        potentialMaster.addPotential(potential,new Species[]{species,species});

        criterion = new CriterionSimple(this,potential.getRange(),neighborRangeFac*potential.getRange());
        potential2.setCriterion(criterion);
        potentialMaster.addPotential(potential2,new Species[]{species2,species2});

        criterion = new CriterionSimple(this,potential.getRange(),neighborRangeFac*potential.getRange());
        potential22.setCriterion(criterion);
        potentialMaster.addPotential(potential22,new Species[]{species2,species});
//        potentialMaster.setSpecies(potential,new Species[]{species,species});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species2});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species});
//        integrator.addIntervalListener(new PhaseImposePbc(phase));

        
//		elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        phase = new Phase(this);
        new ConfigurationSequential(space).initializeCoordinates(phase);
        integrator.setPhase(phase);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
		sim.getController().actionPerformed();
    }//end of main
    
}