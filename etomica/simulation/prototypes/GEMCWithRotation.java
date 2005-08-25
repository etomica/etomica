package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.Configuration;
import etomica.config.ConfigurationSequential;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresRotating;

/**
 * Simple Gibbs-ensemble Monte Carlo simulation of rotating molecules.
 */

//in present form uses just a LJ potential, so orientation is irrelevant

public class GEMCWithRotation extends Simulation {
    
    public GEMCWithRotation() {
        this(Space2D.getInstance());
    }
    
    public GEMCWithRotation(Space space) {
        super(space, false, new PotentialMaster(space));
        defaults.atomSize = 1.2;
        defaults.unitSystem = new etomica.units.systems.LJ();
        defaults.temperature = defaults.unitSystem.temperature().toSim(0.420);
        IntegratorGEMC integrator = new IntegratorGEMC(this);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        getController().addAction(activityIntegrate);
	    activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
	    activityIntegrate.setInterval(400);
	    
	    species = new SpeciesSpheresRotating(this);
        species.setNMolecules(200);

	    phase1 = new Phase(this);
        
	    MCMoveRotate mcRotate1 = new MCMoveRotate(potentialMaster, space);
	    mcRotate1.setPhase(new Phase[] {phase1});
        integrator.addMCMove(mcRotate1);

	    phase2 = new Phase(this);
	    //MeterDensity meter2 = new MeterDensity();	    
	    MCMoveRotate mcRotate2 = new MCMoveRotate(potentialMaster, space);
	    mcRotate2.setPhase(new Phase[] {phase2});
        integrator.addMCMove(mcRotate2);

        Configuration config;
        if (space.D() == 2) {
            config = new ConfigurationSequential(space);
        }
        else {
            config = new ConfigurationSequential(space);
        }
        config.initializeCoordinates(phase1);
        config.initializeCoordinates(phase2);
            
        integrator.setPhase(new Phase[] {phase1, phase2});
        
        //MeterDensity meterDensity = new MeterDensity();

	    potential = new P2LennardJones(this);
	    potential.setSigma(species.getDiameter());

        this.potentialMaster.setSpecies(potential,new Species[] {species, species});

        integrator.addListener(new PhaseImposePbc(phase1));
        integrator.addListener(new PhaseImposePbc(phase2));

	    phase2.setDensity(0.1);
    }//end of constructor        
    
    public Phase phase1, phase2;
    public IntegratorGEMC integrator;
    public SpeciesSpheresRotating species;
    public P2LennardJones potential;
}