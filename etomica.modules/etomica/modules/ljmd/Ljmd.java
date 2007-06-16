package etomica.modules.ljmd;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class Ljmd extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public Phase phase;
    public IntegratorVelocityVerlet integrator;
    public ActivityIntegrate activityIntegrate;
    
    public Ljmd(Space space) {
        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        
        int N = 182;  //number of atoms
        
        //controller and integrator
	    integrator = new IntegratorVelocityVerlet(this, potentialMaster);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(0.01);
     //   integrator.setDoSleep(false);

	    //species and potentials
	//    SpeciesSphereWells disks = new SpeciesSphereWells(this);//index 1
	    species = new SpeciesSpheresMono(this);//index 1
        getSpeciesManager().addSpecies(species);
	    species.setName("");
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(potential,2.5);
	    potentialMaster.addPotential(p2Truncated, new Species[]{species, species});
	    
        //construct phase
	    phase = new Phase(this);
        addPhase(phase);
        IVector dim = space.makeVector();
        dim.E(14);
        phase.setDimensions(dim);
        phase.getAgent(species).setNMolecules(N);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(phase);
        integrator.setPhase(phase);
		
        PhaseImposePbc imposePBC = new PhaseImposePbc(phase);
        integrator.addIntervalAction(imposePBC);
        
    }//end of constructor    
    
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
