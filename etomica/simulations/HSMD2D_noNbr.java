package etomica.simulations;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.potential.P2HardSphere;
import etomica.space2d.Space2D;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D_noNbr extends Simulation {
    
    public ActivityIntegrate activityIntegrate;

    public HSMD2D_noNbr() {
    	this(new Space2D());
    }
    
    public HSMD2D_noNbr(Space2D space) {
        super(space);

        IntegratorHard integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        species.setNMolecules(64);
	    Phase phase = new Phase(space);
	    P2HardSphere potential = new P2HardSphere(space);
	    potentialMaster.setSpecies(potential,new Species[]{species,species});
        
        integrator.addIntervalListener(new PhaseImposePbc(phase));
        phase.speciesMaster.addSpecies(species);
        integrator.addPhase(phase);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D_noNbr sim = new HSMD2D_noNbr();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.activityIntegrate.setDoSleep(true);
		graphic.makeAndDisplayFrame();
    }//end of main
    
}