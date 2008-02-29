package etomica.modules.swmd;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IVector;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2HardSphericalWrapper;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dalton;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.UnitRatio;

public class Swmd extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorHard integrator;
    public Potential2HardSphericalWrapper potentialWrapper;
    public ActivityIntegrate activityIntegrate;
    
    public Swmd(Space _space) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);
        
        int N = space.D() == 3 ? 256 : 100;  //number of atoms
        
        double sigma = 4.0;
        double lambda = 2.0;
        
        //controller and integrator
	    integrator = new IntegratorHard(potentialMaster, getRandom(), 1.0, Kelvin.UNIT.toSim(300), space);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        P1HardPeriodic nullPotential = new P1HardPeriodic(space, sigma*lambda);
        integrator.setNullPotential(nullPotential);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(1);

	    //species and potentials
	    species = new SpeciesSpheresMono(this);//index 1
	    ((ElementSimple)species.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(40));
        getSpeciesManager().addSpecies(species);
        
        //instantiate several potentials for selection in combo-box
	    P2SquareWell potentialSW = new P2SquareWell(space, sigma, lambda, new UnitRatio(Joule.UNIT, Mole.UNIT).toSim(1500), true);
        potentialWrapper = new Potential2HardSphericalWrapper(space,potentialSW);
        potentialMaster.addPotential(potentialWrapper,new AtomType[]{species.getLeafType(),species.getLeafType()});
	    
        //construct box
	    box = new Box(this, space);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(space.D() == 3 ? 30 : 50);
        box.setDimensions(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc() : new LatticeOrthorhombicHexagonal(), space).initializeCoordinates(box);
        integrator.setBox(box);

        integrator.addIntervalAction(new BoxImposePbc(box, space));
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        Swmd sim = new Swmd(space);
        sim.getController().actionPerformed();
    }
}
