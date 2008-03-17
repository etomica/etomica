package etomica.modules.osmosis;

import java.awt.Color;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicSimple;
import etomica.math.geometry.Plane;
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
    public IBox box;
    public P2HardSphere potentialAA,potentialBB,potentialAB;
    public P1HardBoundary boundaryHardA;
    public P1HardBoundary boundaryHardB;
    public P1HardWall boundarySemiB;
    public ActivityIntegrate activityIntegrate;

    public OsmosisSim(Space _space) {

        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster(space);

        final double sigma = 1.0;

	    speciesSolvent = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(speciesSolvent);
	    speciesSolute = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(speciesSolute);
        ((AtomTypeSphere)speciesSolvent.getLeafType()).setDiameter(sigma);
        ((AtomTypeSphere)speciesSolute.getLeafType()).setDiameter(sigma);

	    potentialAA = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialAA, new AtomType[]{speciesSolvent.getLeafType(), speciesSolvent.getLeafType()});
	    potentialBB = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialBB, new AtomType[]{speciesSolute.getLeafType(), speciesSolute.getLeafType()});
	    potentialAB = new P2HardSphere(space, sigma, true);
        potentialMaster.addPotential(potentialAB, new AtomType[]{speciesSolvent.getLeafType(), speciesSolute.getLeafType()});
        
	    //Boundary potential for the solvent
        boundaryHardA = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardA, new AtomType[]{speciesSolvent.getLeafType()});
        boundaryHardA.setCollisionRadius(0.5*sigma);
        
        //Boundary potential for the solute
        boundaryHardB = new P1HardBoundary(space, true);
        potentialMaster.addPotential(boundaryHardB, new AtomType[]{speciesSolute.getLeafType()});
        boundaryHardB.setCollisionRadius(0.5*sigma);

        //wall in the middle that only applies to the solute
	    boundarySemiB = new P1HardWall(space, sigma);
        potentialMaster.addPotential(boundarySemiB, new AtomType[]{speciesSolute.getLeafType()});
	    boundarySemiB.setCollisionRadius(0.5*sigma);
        
        //construct box
	    box = new Box(this, space);
        addBox(box);
        box.setBoundary(new BoundaryRectangularNonperiodic(space, getRandom()));

        if (space instanceof Space2D){ // 2D
            box.getBoundary().setDimensions(new Vector2D(10.0, 10.0));
        }
        else if (space instanceof Space3D) { // 3D
            box.getBoundary().setDimensions(new Vector3D(10.0, 10.0, 10.0));
        }
        box.setNMolecules(speciesSolvent, initialSolvent);
        box.setNMolecules(speciesSolute, initialSolute);

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setBox(box);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
    }

    public static void main(String[] args) {
    	Space sp = Space3D.getInstance();
        final OsmosisSim sim = new OsmosisSim(sp);
        final ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(sp, 1.0), sp);
        config.initializeCoordinates(sim.box);
    	Plane plane = new Plane(sim.getSpace());

        final SimulationGraphic simGraphic = new SimulationGraphic(sim, "Osmosis Sim", sp);
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
        sim.activityIntegrate.setSleepPeriod(1);
    }

} 
