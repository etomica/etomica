package etomica.simulations;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorHard;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionSimple;
import etomica.nbr.PotentialMasterNbr;
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

    public HSMD2D() {
    	this(Space2D.INSTANCE);
    }
    
    public HSMD2D(Space2D space) {
        super(space, new PotentialMasterNbr(space));
//        super(space, new PotentialMaster(space,IteratorFactoryCell.INSTANCE));
        Default.makeLJDefaults();
//        Default.BOX_SIZE = 30.0;
        Default.ATOM_SIZE = 0.38;

        double neighborRangeFac = 1.6;
        int nCells = (int)(Default.BOX_SIZE/neighborRangeFac);
        System.out.println("nCells: "+nCells);
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);
        ((PotentialMasterNbr)potentialMaster).setMaxNeighborRange(neighborRangeFac);

        integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(false);
        integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
	    species2 = new SpeciesSpheresMono(this);
	    species.setNMolecules(512);
	    species2.setNMolecules(5);
	    phase = new Phase(space);
 //       ((PotentialMasterNbr)potentialMaster).getNbrCellManager(phase);
	    potential = new P2HardSphere(space);
//	    this.potentialMaster.setSpecies(potential,new Species[]{species,species});
	    //this.potentialMaster.setSpecies(potential,new Species[]{species2,species2});
	    //this.potentialMaster.setSpecies(potential,new Species[]{species,species2});
        
 //       integrator.addIntervalListener(new NeighborCellManager(phase,15));
        NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        ((PotentialMasterNbr)potentialMaster).setSpecies(potential,new Species[]{species,species},criterion);
        ((PotentialMasterNbr)potentialMaster).setSpecies(potential,new Species[]{species2,species2},criterion);
        ((PotentialMasterNbr)potentialMaster).setSpecies(potential,new Species[]{species2,species},criterion);
        
//		elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        phase.speciesMaster.addSpecies(species);
        phase.speciesMaster.addSpecies(species2);
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