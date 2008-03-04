package etomica.modules.ljmd;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomType;
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

public class Ljmd extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public IBox box;
    public IntegratorVelocityVerlet integrator;
    public ActivityIntegrate activityIntegrate;
    
    public Ljmd(Space _space) {
        super(_space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, 2.99, space);
        
        int N = 182;  //number of atoms
        
        //controller and integrator
	    integrator = new IntegratorVelocityVerlet(this, potentialMaster, space);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(0.01);
     //   integrator.setDoSleep(false);

	    //species and potentials
	//    SpeciesSphereWells disks = new SpeciesSphereWells(this);//index 1
	    species = new SpeciesSpheresMono(this);//index 1
        getSpeciesManager().addSpecies(species);
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(potential,2.5);
	    potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(this, space);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(15);
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
            
        Ljmd sim = new Ljmd(space);
        sim.getController().actionPerformed();
    }//end of main
    
}
