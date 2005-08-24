package etomica.simulations;
import etomica.Configuration;
import etomica.ConfigurationSequential;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresRotating;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.meter.MeterDensity;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.space2d.Space2D;

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
        Default.ATOM_SIZE = 1.2;
        Default.UNIT_SYSTEM = new etomica.units.systems.LJ();
        Default.TEMPERATURE = Default.UNIT_SYSTEM.temperature().toSim(0.420);
        IntegratorGEMC integrator = new IntegratorGEMC(potentialMaster, space);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
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
	    MeterDensity meter2 = new MeterDensity();	    
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
        
        MeterDensity meterDensity = new MeterDensity();

	    potential = new P2LennardJones(space);
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