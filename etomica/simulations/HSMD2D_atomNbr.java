package etomica.simulations;
import etomica.Default;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionSimple;
import etomica.nbratom.PotentialMasterNbr;
import etomica.nbratom.cell.NeighborCellManager;
import etomica.nbratom.cell.PotentialCalculationCellAssign;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.space2d.Space2D;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D_atomNbr extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public Potential2 potential;

    public HSMD2D_atomNbr() {
    	this(Space2D.INSTANCE);
    }
    
    public HSMD2D_atomNbr(Space2D space) {
        super(space, new PotentialMasterNbr(space));
        Default.makeLJDefaults();
        Default.ATOM_SIZE = 0.38;

        double neighborRangeFac = 1.6;
        int nCells = (int)(Default.BOX_SIZE/neighborRangeFac);
        System.out.println("nCells: "+nCells);
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);

        integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(false);
        integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
//	    species2 = new SpeciesSpheresMono(this);
	    species.setNMolecules(512);
//	    species2.setNMolecules(5);
	    phase = new Phase(this);
	    potential = new P2HardSphere(space);
        
        integrator.addIntervalListener(new NeighborCellManager(phase,15));
        NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());

        criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        ((PotentialMasterNbr)potentialMaster).getNeighborManager().addCriterion(criterion);
        potential.setCriterion(criterion);
        potentialMaster.setSpecies(potential,new Species[]{species,species});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species2});
//        potentialMaster.setSpecies(potential,new Species[]{species2,species});
//        integrator.addIntervalListener(new PhaseImposePbc(phase));
        
//		elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
//        phase.speciesMaster.addSpecies(species2);
        integrator.addPhase(phase);
        
        PotentialCalculationCellAssign pc = new PotentialCalculationCellAssign(((PotentialMasterNbr)potentialMaster).getNbrCellManager(phase));
        ((PotentialMasterNbr)potentialMaster).calculate(phase, new IteratorDirective(),pc); 
                
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
      Default.DO_SLEEP = false;
      HSMD2D_atomNbr sim = new etomica.simulations.HSMD2D_atomNbr();
      SimulationGraphic simGraphic = new SimulationGraphic(sim);
      DeviceNSelector nSelector = new DeviceNSelector(sim,sim.phase.getAgent(sim.species));
      simGraphic.add(nSelector);
      simGraphic.makeAndDisplayFrame();
      ColorSchemeByType.setColor(sim.species.getFactory().getType(), java.awt.Color.red);
//      ColorSchemeByType.setColor(sim.species2, java.awt.Color.blue);
      simGraphic.panel().setBackground(java.awt.Color.yellow);
    }//end of main
    
}