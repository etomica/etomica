package etomica.simulation.prototypes;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeSphere;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationSequential;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
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
        defaults.makeLJDefaults();
        defaults.temperature = defaults.unitSystem.temperature().toSim(0.420);
        integrator = new IntegratorGEMC(potentialMaster);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        getController().addAction(activityIntegrate);
	    activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
	    activityIntegrate.setInterval(400);
	    
	    species = new SpeciesSpheresRotating(this);

	    phase1 = new Phase(this);
        phase1.getAgent(species).setNMolecules(200);
        
	    IntegratorMC integratorMC = new IntegratorMC(this);
        integratorMC.setPhase(phase1);
        MCMoveManager moveManager = integratorMC.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, space));
        moveManager.addMCMove(new MCMoveAtom(this));
        integrator.addIntegrator(integratorMC);
        

	    phase2 = new Phase(this);
        phase2.getAgent(species).setNMolecules(200);
        integratorMC = new IntegratorMC(this);
        integratorMC.setPhase(phase1);
        moveManager = integratorMC.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, space));
        moveManager.addMCMove(new MCMoveAtom(this));
        // GEMC integrator adds volume and molecule exchange moves once
        // it has 2 integrators
        integrator.addIntegrator(integratorMC);
        
        Configuration config;
        if (space.D() == 2) {
            config = new ConfigurationSequential(space);
        }
        else {
            config = new ConfigurationLattice(new LatticeCubicFcc());
        }
        config.initializeCoordinates(phase1);
        config.initializeCoordinates(phase2);
            
	    potential = new P2LennardJones(this);
	    potential.setSigma(((AtomTypeSphere)species.getMoleculeType()).getDiameter());

        this.potentialMaster.addPotential(potential,new Species[] {species, species});

        integrator.addListener(new PhaseImposePbc(phase1));
        integrator.addListener(new PhaseImposePbc(phase2));

	    phase2.setDensity(0.1);
    }//end of constructor        
    
    public Phase phase1, phase2;
    public IntegratorGEMC integrator;
    public SpeciesSpheresRotating species;
    public P2LennardJones potential;
}