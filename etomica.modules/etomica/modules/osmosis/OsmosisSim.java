package etomica.modules.osmosis;

import java.awt.Color;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeSphere;
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
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
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
    public P1HardBoundary boundaryHardTopBottomB, boundaryHardLeftB, boundaryHardRightB;
    public P1HardWall boundarySemiB;
    public ActivityIntegrate activityIntegrate;

    public OsmosisSim(Space space) {

        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(space);

        final double sigma = 1.0;

	    speciesA = new SpeciesSpheresMono(this);
	    speciesA.setName("Solvent");
        getSpeciesManager().addSpecies(speciesA);
	    speciesB = new SpeciesSpheresMono(this);
	    speciesB.setName("Solute");
        getSpeciesManager().addSpecies(speciesB);
        ((AtomTypeSphere)speciesA.getMoleculeType()).setDiameter(sigma);
        ((AtomTypeSphere)speciesB.getMoleculeType()).setDiameter(sigma);

	    potentialAA = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialAA, new Species[]{speciesA, speciesA});
	    potentialBB = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialBB, new Species[]{speciesB, speciesB});
	    potentialAB = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialAB, new Species[]{speciesA, speciesB});
        
        //boundary for solvent on the top and bottom
	    boundaryHardTopBottomA = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardTopBottomA, new Species[]{speciesA});
        //disable left and right
        boundaryHardTopBottomA.setActive(0, true, false);
        boundaryHardTopBottomA.setActive(0, false, false);
	    boundaryHardTopBottomA.setCollisionRadius(0.5*sigma);
        
	    //left and right boundaries need to be separate so we can measure the
        //force on each and get an osmotic pressure
        //SOLVENT
        // left boundary 
        boundaryHardLeftA = new P1HardBoundary(space, true);
        //disable right
        boundaryHardLeftA.setActive(0, false, false);
        //disable top and bottom
        boundaryHardLeftA.setActive(1, true, false);
        boundaryHardLeftA.setActive(1, false, false);
        if (space.D() == 3) {
            //disable front and back
            boundaryHardLeftA.setActive(2, true, false);
            boundaryHardLeftA.setActive(2, false, false);
        }
        potentialMaster.addPotential(boundaryHardLeftA, new Species[]{speciesA});
        boundaryHardLeftA.setCollisionRadius(0.5*sigma);
        // right boundary
        boundaryHardRightA = new P1HardBoundary(space, true);
        //disable left
        boundaryHardRightA.setActive(0, true, false);
        //disable top and bottom
        boundaryHardRightA.setActive(1, true, false);
        boundaryHardRightA.setActive(1, false, false);
        if (space.D() == 3) {
            //disable front and back
            boundaryHardRightA.setActive(2, true, false);
            boundaryHardRightA.setActive(2, false, false);
        }
        potentialMaster.addPotential(boundaryHardRightA, new Species[]{speciesA});
        boundaryHardRightA.setCollisionRadius(0.5*sigma);
        
        boundaryHardTopBottomA = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardTopBottomA, new Species[]{speciesA});
        //disable left and right
        boundaryHardTopBottomA.setActive(0, true, false);
        boundaryHardTopBottomA.setActive(0, false, false);
        boundaryHardTopBottomA.setCollisionRadius(0.5*sigma);
        
        //left and right boundaries need to be separate so we can measure the
        //force on each and get an osmotic pressure
        //SOLUTE
        // left boundary
        boundaryHardLeftB = new P1HardBoundary(space, true);
        //disable right
        boundaryHardLeftB.setActive(0, false, false);
        //disable top and bottom
        boundaryHardLeftB.setActive(1, true, false);
        boundaryHardLeftB.setActive(1, false, false);
        if (space.D() == 3) {
            //disable front and back
            boundaryHardLeftB.setActive(2, true, false);
            boundaryHardLeftB.setActive(2, false, false);
        }
        potentialMaster.addPotential(boundaryHardLeftB, new Species[]{speciesB});
        boundaryHardLeftB.setCollisionRadius(0.5*sigma);
        // right boundary
        boundaryHardRightB = new P1HardBoundary(space, true);
        //disable left
        boundaryHardRightB.setActive(0, true, false);
        //disable top and bottom
        boundaryHardRightB.setActive(1, true, false);
        boundaryHardRightB.setActive(1, false, false);
        if (space.D() == 3) {
            //disable front and back
            boundaryHardRightB.setActive(2, true, false);
            boundaryHardRightB.setActive(2, false, false);
        }
        potentialMaster.addPotential(boundaryHardRightB, new Species[]{speciesB});
        boundaryHardRightB.setCollisionRadius(0.5*sigma);

        boundaryHardTopBottomB = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardTopBottomB, new Species[]{speciesB});
        //disable left and right
        boundaryHardTopBottomB.setActive(0, true, false);
        boundaryHardTopBottomB.setActive(0, false, false);
        boundaryHardTopBottomB.setCollisionRadius(0.5*sigma);

        //wall in the middle that only applies to the solute
	    boundarySemiB = new P1HardWall(space, sigma);
        potentialMaster.addPotential(boundarySemiB, new Species[]{speciesB});
	    boundarySemiB.setCollisionRadius(0.5*sigma);
        
        //construct phase
	    phase = new Phase(this);
        addPhase(phase);
        phase.setBoundary(new BoundaryRectangularNonperiodic(space, getRandom()));

        ConfigurationLattice config = null;
        if (space instanceof Space2D){ // 2D
            phase.getBoundary().setDimensions(new Vector2D(10.0, 10.0));
            config = new ConfigurationLattice(new LatticeCubicSimple(2, 1.0));
        }
        else if (space instanceof Space3D) { // 3D
            phase.getBoundary().setDimensions(new Vector3D(10.0, 10.0, 10.0));
        	config = new ConfigurationLattice(new LatticeCubicSimple(3, 1.0));
        }
        phase.getAgent(speciesA).setNMolecules(50);
        phase.getAgent(speciesB).setNMolecules(10);
        config.initializeCoordinates(phase);

        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setPhase(phase);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);

        activityIntegrate = new ActivityIntegrate(integrator, false, false);
        getController().addAction(activityIntegrate);
    }

    public static void main(String[] args) {
        OsmosisSim sim = new OsmosisSim(Space2D.getInstance());
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
