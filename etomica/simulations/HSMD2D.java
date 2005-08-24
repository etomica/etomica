package etomica.simulations;
import etomica.Simulation;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.config.ConfigurationSequential;
import etomica.integrator.IntegratorHard;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
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
    	this(Space2D.getInstance());
    }
    
    public HSMD2D(Space2D space) {
        super(space, true, new PotentialMasterNbr(space));
//        super(space, new PotentialMaster(space));//,IteratorFactoryCell.instance));
        Default.makeLJDefaults();
//        Default.BOX_SIZE = 30.0;
        Default.ATOM_SIZE = 0.38;

        double neighborRangeFac = 1.6;
        int nCells = (int)(Default.BOX_SIZE/neighborRangeFac);
        System.out.println("nCells: "+nCells);
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);

        integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        NeighborListManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
        nbrManager.setRange(Default.ATOM_SIZE*1.6);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
	    species2 = new SpeciesSpheresMono(this);
        species.setNMolecules(512);
        species2.setNMolecules(5);
        potential = new P2HardSphere(space);
        potential2 = new P2HardSphere(space);
        potential22 = new P2HardSphere(space);
//	    this.potentialMaster.setSpecies(potential,new Species[]{species,species});
	    //this.potentialMaster.setSpecies(potential,new Species[]{species2,species2});
	    //this.potentialMaster.setSpecies(potential,new Species[]{species,species2});
        
        NeighborCriterion criterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        potential.setCriterion(criterion);
        potentialMaster.setSpecies(potential,new Species[]{species,species});

        criterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        potential2.setCriterion(criterion);
        potentialMaster.setSpecies(potential2,new Species[]{species2,species2});

        criterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        potential22.setCriterion(criterion);
        potentialMaster.setSpecies(potential22,new Species[]{species2,species});
        nbrManager.addCriterion(criterion,new AtomType[]{species.getFactory().getType(),species2.getFactory().getType()});
//        potentialMaster.setSpecies(potential,new Species[]{species,species});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species2});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species});
//        integrator.addIntervalListener(new PhaseImposePbc(phase));

        
//		elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        phase = new Phase(this);
        new ConfigurationSequential(space).initializeCoordinates(phase);
        integrator.addPhase(phase);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
		sim.getController().start();
    }//end of main
    
}