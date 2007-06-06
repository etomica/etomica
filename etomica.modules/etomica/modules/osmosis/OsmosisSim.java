package etomica.modules.osmosis;

import java.awt.Color;

import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicSimple;
import etomica.phase.Phase;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Osmosis simulation.
 * @author Jhumpa Adhikari
 * @author Andrew Schultz
 */

public class OsmosisSim extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorHard integrator;
    public SpeciesSpheresMono speciesA,speciesB;
    public Phase phase;
    public P2HardSphere potentialAA,potentialBB,potentialAB;
    public P1HardBoundary boundaryHardTopBottomA, boundaryHardLeftA, boundaryHardRightA;
    public P1HardBoundary boundaryHardB;
    public P1HardWall boundarySemiB;
    public ActivityIntegrate activityIntegrate;
    
    public OsmosisSim() {

        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(this);

        defaults.ignoreOverlap = true;
        final double sigma = defaults.atomSize;

	    speciesA = new SpeciesSpheresMono(this);
	    speciesA.setName("Solvent");
        getSpeciesManager().addSpecies(speciesA);
	    speciesB = new SpeciesSpheresMono(this);
	    speciesB.setName("Solute");
        getSpeciesManager().addSpecies(speciesB);

	    potentialAA = new P2HardSphere(this);
        potentialMaster.addPotential(potentialAA, new Species[]{speciesA, speciesA});
	    potentialBB = new P2HardSphere(this);
        potentialMaster.addPotential(potentialBB, new Species[]{speciesB, speciesB});
	    potentialAB = new P2HardSphere(this);
        potentialMaster.addPotential(potentialAB, new Species[]{speciesA, speciesB});
        
	    boundaryHardTopBottomA = new P1HardBoundary(this);
        potentialMaster.addPotential(boundaryHardTopBottomA, new Species[]{speciesA});
        boundaryHardTopBottomA.setActive(0, true, false);
        boundaryHardTopBottomA.setActive(0, false, false);
	    boundaryHardTopBottomA.setCollisionRadius(0.5*sigma);
        boundaryHardLeftA = new P1HardBoundary(this);
        boundaryHardLeftA.setActive(0, false, false);
        boundaryHardLeftA.setActive(1, true, false);
        boundaryHardLeftA.setActive(1, false, false);
        potentialMaster.addPotential(boundaryHardLeftA, new Species[]{speciesA});
        boundaryHardLeftA.setCollisionRadius(0.5*sigma);
        boundaryHardRightA = new P1HardBoundary(this);
        boundaryHardLeftA.setActive(0, true, false);
        boundaryHardLeftA.setActive(1, true, false);
        boundaryHardLeftA.setActive(1, false, false);
        potentialMaster.addPotential(boundaryHardRightA, new Species[]{speciesA});
        boundaryHardRightA.setCollisionRadius(0.5*sigma);
        
	    boundaryHardB = new P1HardBoundary(this);
        potentialMaster.addPotential(boundaryHardB, new Species[]{speciesB});
	    boundaryHardB.setCollisionRadius(0.5*sigma);
        
	    boundarySemiB = new P1HardWall(this);
        potentialMaster.addPotential(boundarySemiB, new Species[]{speciesB});
	    boundarySemiB.setCollisionRadius(0.5*sigma);
        
        //construct phase
	    phase = new Phase(this);
        addPhase(phase);
        phase.setBoundary(new BoundaryRectangularNonperiodic(space, getRandom()));
        phase.getBoundary().setDimensions(new Vector2D(30.0, 30.0));
        phase.getAgent(speciesA).setNMolecules(30);
        phase.getAgent(speciesB).setNMolecules(10);
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(2, 1.0));
        config.initializeCoordinates(phase);

        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setPhase(phase);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        
        activityIntegrate = new ActivityIntegrate(integrator, false, false);
        getController().addAction(activityIntegrate);
    }

    public static void main(String[] args) {
        OsmosisSim sim = new OsmosisSim();
        sim.register(sim.integrator);
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.speciesA.getMoleculeType(), Color.blue);
        colorScheme.setColor(sim.speciesB.getMoleculeType(), Color.red);
        simGraphic.getDisplayPhase(sim.phase).setColorScheme(colorScheme);
        sim.integrator.setTimeStep(0.05);
        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(1);
    }

} 
