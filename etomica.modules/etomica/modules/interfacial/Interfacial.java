package etomica.modules.interfacial;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation for interfacial tension module.  Simulation itself is just a
 * simple LJ system.
 *
 * @author Andrew Schultz
 */
public class Interfacial extends Simulation {

    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public IBox box;
    public IntegratorVelocityVerlet integrator;
    public ActivityIntegrate activityIntegrate;

    public Interfacial(Space _space) {
        super(_space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, 5, space);

        int N = 250;  //number of atoms

        //controller and integrator
	    integrator = new IntegratorVelocityVerlet(this, potentialMaster, space);
	    integrator.setIsothermal(true);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(0.01);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);//index 1
        getSpeciesManager().addSpecies(species);
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential,4.5);
	    potentialMaster.addPotential(p2Truncated, new IAtomTypeLeaf[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(this, space);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(20);
        box.setDimensions(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(), space).initializeCoordinates(box);
        integrator.setBox(box);

        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));
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
            
        Interfacial sim = new Interfacial(space);
        sim.getController().actionPerformed();
    }//end of main
}
