package etomica.modules.swmd;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P1HardMoleculeMonatomic;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardMoleculeMonatomic;
import etomica.potential.P2SquareWell;
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
    public IBox box;
    public IntegratorHard integrator;
    public P2HardMoleculeMonatomic potentialWrapper;
    public ActivityIntegrate activityIntegrate;
    
    public Swmd(Space _space) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);
        
        int N = space.D() == 3 ? 256 : 100;  //number of atoms
        
        double sigma = 4.0;
        double lambda = 2.0;
        
        //controller and integrator
	    integrator = new IntegratorHard(this, potentialMaster, getRandom(), 1.0, Kelvin.UNIT.toSim(300), space, true);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        P1HardPeriodic nullPotential = new P1HardPeriodic(space, sigma*lambda);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(1);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);//index 1
	    ((ElementSimple)species.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(space.D() == 3 ? 131 : 40));
        getSpeciesManager().addSpecies(species);
        integrator.setNullPotential(new P1HardMoleculeMonatomic(space, nullPotential), species.getMoleculeType());
        
        //instantiate several potentials for selection in combo-box
	    P2SquareWell potentialSW = new P2SquareWell(space, sigma, lambda, new UnitRatio(Joule.UNIT, Mole.UNIT).toSim(space.D() == 3 ? 1000 : 1500), true);
        potentialWrapper = new P2HardMoleculeMonatomic(space,potentialSW);
        potentialMaster.addPotential(potentialWrapper,new ISpecies[]{species,species});
	    
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
