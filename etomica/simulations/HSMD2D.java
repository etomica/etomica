package etomica.simulations;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorHard;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.NeighborManager;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.space2d.Space2D;

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
    	this(Space2D.INSTANCE);
    }
    
    public HSMD2D(Space2D space) {
        super(space, new PotentialMasterNbr(space));
//        super(space, new PotentialMaster(space));//,IteratorFactoryCell.INSTANCE));
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

        NeighborManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
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
        nbrManager.addCriterion(criterion);
        species.getFactory().getType().getNbrManagerAgent().addCriterion(criterion);

        criterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        potential2.setCriterion(criterion);
        potentialMaster.setSpecies(potential2,new Species[]{species2,species2});
        nbrManager.addCriterion(criterion);
        species2.getFactory().getType().getNbrManagerAgent().addCriterion(criterion);

        criterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        potential22.setCriterion(criterion);
        potentialMaster.setSpecies(potential22,new Species[]{species2,species});
        species.getFactory().getType().getNbrManagerAgent().addCriterion(criterion);
        species2.getFactory().getType().getNbrManagerAgent().addCriterion(criterion);
//        potentialMaster.setSpecies(potential,new Species[]{species,species});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species2});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species});
//        integrator.addIntervalListener(new PhaseImposePbc(phase));

        
//		elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        phase = new Phase(this);
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