package etomica.modules.joulethomson;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationSequential;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGear4NPH;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class JouleThomsonSim extends Simulation {
    
    IntegratorGear4NPH integrator;
    IntegratorMD integratorNVE;
    SpeciesSpheresMono species;
    P2LennardJones potential;
    Phase phase;
    IntegratorJT integratorJT;
    ActivityIntegrate activityIntegrate;
    
    public JouleThomsonSim() {this(Space2D.getInstance());}
    public JouleThomsonSim(Space space) {
        super(space);
        int nAtoms = (space.D() < 3) ? 50 : 64;
        
		defaults.historyPeriod = 1000;
 
        //integrator
        integratorNVE = new IntegratorVelocityVerlet(this);
        integrator = new IntegratorGear4NPH(this);
        integrator.setRelaxationRateP(500.);
        integrator.setRelaxationRateH(300.);
        integrator.setTargetP(100.0);
	    
	    //species and potential
	    species = new SpeciesSpheresMono(this);
	    potential = new P2LennardJones(this);
        potentialMaster.addPotential(potential, new Species[]{species, species});
	    phase = new Phase(this);
        phase.getAgent(species).setNMolecules(nAtoms);
        
        Configuration config;
        if (space.D() == 2) {
            config = new ConfigurationSequential(space);
        }
        else {
            config = new ConfigurationLattice(new LatticeCubicFcc());
        }
        config.initializeCoordinates(phase);
        
        integratorJT = new IntegratorJT(potentialMaster, integrator, integratorNVE);
        integratorJT.addListener(new PhaseImposePbc(phase));
        integrator.setPhase(phase);
        integratorNVE.setPhase(phase);

        integrator.setTimeStep(0.005);
        integratorNVE.setTimeStep(0.005);

        activityIntegrate = new ActivityIntegrate(this, integratorJT);
        getController().addAction(activityIntegrate);
    }
    
    public static void main(String[] args) {
        JouleThomsonSim sim = new JouleThomsonSim(Space3D.getInstance());
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(10);
    }
}