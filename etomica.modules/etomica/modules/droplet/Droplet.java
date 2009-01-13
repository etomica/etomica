package etomica.modules.droplet;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Mesoscale simulation for Droplet module.
 *
 * @author Andrew Schultz
 */
public class Droplet extends Simulation {

    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono species;
    public final IBox box;
    public final IntegratorDroplet integrator;
    public final ActivityIntegrate activityIntegrate;
    public final P2Cohesion p2;
    public final P1Smash p1Smash;
    public final ConfigurationDroplet config;

    public Droplet(Space _space) {
        super(_space);
        int numAtoms = 2000;
        PotentialMaster potentialMaster = new PotentialMaster();

        //controller and integrator
	    integrator = new IntegratorDroplet(this, potentialMaster, space);
        activityIntegrate = new ActivityIntegrate(integrator);
//        activityIntegrate.setMaxSteps(10);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(0.2);
        integrator.setTemperature(0);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        IAtomTypeLeaf leafType = species.getLeafType();
        
        p2 = new P2Cohesion(space);
        p2.setEpsilon(1.0);
        double vol = 4.0/3.0*Math.PI;
        p2.setDv(vol/numAtoms);
        potentialMaster.addPotential(p2, new IAtomTypeLeaf[]{leafType,leafType});

        p1Smash = new P1Smash(space);
        p1Smash.setG(1.0);
        potentialMaster.addPotential(p1Smash, new IAtomTypeLeaf[]{leafType});

        //construct box
	    box = new Box(new BoundaryRectangularNonperiodic(space), space);
        addBox(box);
        IVectorMutable dim = space.makeVector();
        dim.E(new double[]{10,10,10});
        box.getBoundary().setDimensions(dim);
        box.setNMolecules(species, numAtoms);
        integrator.setBox(box);

        config = new ConfigurationDroplet(random);
        config.initializeCoordinates(box);
    }
    
    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        Droplet sim = new Droplet(space);
        sim.getController().actionPerformed();
    }//end of main
}
