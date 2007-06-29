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
import etomica.math.geometry.Plane;
import etomica.box.Box;
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
    protected final static int initialSolute = 10;
    protected final static int initialSolvent = 50;

    public IntegratorHard integrator;
    public SpeciesSpheresMono speciesSolvent,speciesSolute;
    public Box box;
    public P2HardSphere potentialAA,potentialBB,potentialAB;
    public P1HardBoundary boundaryHardA;
    public P1HardBoundary boundaryHardB;
    public P1HardWall boundarySemiB;
    public ActivityIntegrate activityIntegrate;

    public OsmosisSim(Space space) {

        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(space);

        final double sigma = 1.0;

	    speciesSolvent = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(speciesSolvent);
	    speciesSolute = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(speciesSolute);
        ((AtomTypeSphere)speciesSolvent.getMoleculeType()).setDiameter(sigma);
        ((AtomTypeSphere)speciesSolute.getMoleculeType()).setDiameter(sigma);

	    potentialAA = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialAA, new Species[]{speciesSolvent, speciesSolvent});
	    potentialBB = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialBB, new Species[]{speciesSolute, speciesSolute});
	    potentialAB = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialAB, new Species[]{speciesSolvent, speciesSolute});
        
	    //Boundary potential for the solvent
        boundaryHardA = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardA, new Species[]{speciesSolvent});
        boundaryHardA.setCollisionRadius(0.5*sigma);
        
        //Boundary potential for the solute
        boundaryHardB = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardB, new Species[]{speciesSolute});
        boundaryHardB.setCollisionRadius(0.5*sigma);

        //wall in the middle that only applies to the solute
	    boundarySemiB = new P1HardWall(space, sigma);
        potentialMaster.addPotential(boundarySemiB, new Species[]{speciesSolute});
	    boundarySemiB.setCollisionRadius(0.5*sigma);
        
        //construct box
	    box = new Box(this);
        addBox(box);
        box.setBoundary(new BoundaryRectangularNonperiodic(space, getRandom()));

        if (space instanceof Space2D){ // 2D
            box.getBoundary().setDimensions(new Vector2D(10.0, 10.0));
        }
        else if (space instanceof Space3D) { // 3D
            box.getBoundary().setDimensions(new Vector3D(10.0, 10.0, 10.0));
        }
        box.getAgent(speciesSolvent).setNMolecules(initialSolvent);
        box.getAgent(speciesSolute).setNMolecules(initialSolute);

        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setBox(box);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);

        activityIntegrate = new ActivityIntegrate(integrator, false, false);
        getController().addAction(activityIntegrate);
    }

    public static void main(String[] args) {
        final OsmosisSim sim = new OsmosisSim(Space3D.getInstance());
        final ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(3, 1.0));
        config.initializeCoordinates(sim.box);
    	Plane plane = new Plane(sim.getSpace());

        final SimulationGraphic simGraphic = new SimulationGraphic(sim, "Osmosis Sim");
    	((etomica.graphics.DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).addPlane(plane);
    	simGraphic.getController().getReinitButton().setPostAction(new etomica.action.Action () {
    		public void actionPerformed() {
    	        config.initializeCoordinates(sim.box);
    			simGraphic.getDisplayBox(sim.box).repaint();
    		}
    	});
        simGraphic.makeAndDisplayFrame("Osmosis Sim");
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.speciesSolvent.getMoleculeType(), Color.blue);
        colorScheme.setColor(sim.speciesSolute.getMoleculeType(), Color.red);
        simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
        config.initializeCoordinates(sim.box);
        simGraphic.getDisplayBox(sim.box).repaint();
        sim.integrator.setTimeStep(0.05);
        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(1);
    }

} 
