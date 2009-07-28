package etomica.modules.catalysis;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dalton;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.UnitRatio;

public class Catalysis extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public IBox box;
    public IntegratorHard integrator;
    public ActivityIntegrate activityIntegrate;
    public P2SquareWell potentialSW;
    
    public Catalysis(Space _space) {
        super(_space);
        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this); //List(this, 2.0);
        
        int N = space.D() == 3 ? 256 : 100;  //number of atoms
        
        double sigma = 4.0;
        double lambda = 1.5;
        
        //controller and integrator
	    integrator = new IntegratorHard(this, potentialMaster, space);
	    integrator.setTimeStep(1.0);
	    integrator.setTemperature(Kelvin.UNIT.toSim(300));
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
        integrator.setNullPotential(nullPotential, species.getLeafType());
        
        //instantiate several potentials for selection in combo-box
	    potentialSW = new P2SquareWell(space, sigma, lambda, new UnitRatio(Joule.UNIT, Mole.UNIT).toSim(space.D() == 3 ? 1000 : 1500), true);
        potentialMaster.addPotential(potentialSW,new IAtomType[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(space);
        addBox(box);
        IVectorMutable dim = space.makeVector();
        dim.E(space.D() == 3 ? 30 : 50);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        integrator.setBox(box);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
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
            
        Catalysis sim = new Catalysis(space);
        sim.getController().actionPerformed();
    }
}
