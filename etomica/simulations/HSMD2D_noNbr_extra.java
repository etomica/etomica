package etomica.simulations;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageSegment;
import etomica.data.AccumulatorAverage.Type;
import etomica.data.meter.MeterPressureHard;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.potential.P2HardSphere;
import etomica.space2d.Space2D;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D_noNbr_extra extends Simulation {
    
    public ActivityIntegrate activityIntegrate;
    public DisplayBox pressureDisplay;

    public HSMD2D_noNbr_extra() {
    	this(new Space2D());
    }
    
    public HSMD2D_noNbr_extra(Space2D space) {
        super(space);

        IntegratorHard integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        species.setNMolecules(64);
	    Phase phase = new Phase(this);
	    P2HardSphere potential = new P2HardSphere(space);
	    potentialMaster.setSpecies(potential,new Species[]{species,species});
        
        integrator.addListener(new PhaseImposePbc(phase));
        integrator.addPhase(phase);
        
        MeterPressureHard pressureMeter = new MeterPressureHard(integrator);
        pressureMeter.setPhase(phase);
        pressureDisplay = new DisplayBox();
        new AccumulatorAverageSegment(pressureMeter, integrator, new Type[] {AccumulatorAverage.AVERAGE}, pressureDisplay);
        
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D_noNbr_extra sim = new HSMD2D_noNbr_extra();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        graphic.add(sim.pressureDisplay);
        sim.activityIntegrate.setDoSleep(true);
		graphic.makeAndDisplayFrame();
    }//end of main
    
}