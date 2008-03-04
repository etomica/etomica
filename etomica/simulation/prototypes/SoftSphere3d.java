package etomica.simulation.prototypes;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.AtomTypeLeaf;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
 * @author David Kofke
 */
 
public class SoftSphere3d extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species, species2;
    public IBox box;
    public P2SoftSphere potential;
    public IPotentialMaster potentialMaster;
    public Controller controller;
    public DataSourceCountSteps meterCycles;
    

    public SoftSphere3d(double density, double softness, double temperature) {
        super(Space3D.getInstance());
        potentialMaster = new PotentialMaster(space);
	    integrator = new IntegratorMC(this, potentialMaster);
	    integrator.setTemperature(temperature);
	    
	    
	    mcMoveAtom = new MCMoveAtom(this, potentialMaster);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(10000000);
        getController().addAction(activityIntegrate);
        
        species = new SpeciesSpheresMono(this);
        //species2 = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        //getSpeciesManager().addSpecies(species2);
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, 108);
        box.setDensity(density);
       // box.setNMolecules(species2, 20);
        new ConfigurationLattice(new LatticeCubicFcc(), space).initializeCoordinates(box);
	    potential = new P2SoftSphere(space,1,1,softness);
	    P2SoftSphericalTruncated truncated = new P2SoftSphericalTruncated(potential,box.getBoundary().getDimensions().x(0)/2);
	   // System.out.println("Truncated radius is: " +truncated.getTruncationRadius());
	    
	    AtomTypeLeaf type1 = species.getLeafType();
        //AtomTypeLeaf type2 = species2.getLeafType();
        potentialMaster.addPotential(truncated, new IAtomType[] {type1, type1});
       // potentialMaster.addPotential(potential, new AtomType[] {type1, type2});
        //potentialMaster.addPotential(potential, new AtomType[] {type2, type2});
        
	    meterCycles = new DataSourceCountSteps(integrator);

        integrator.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.addIntervalAction(new BoxImposePbc(box, space));

//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);
	    
//		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(box);
//        colorSchemeCell.setLattice(lattice);
    }
    
    public static void main(String[] args) {
   
    	double density = 1.338;
    	double softness = 0.1;
    	double temperature = 0.1;
    	
    	
        if (args.length > 0) {
            density = Double.parseDouble(args[0]);
        }
        if (args.length > 1) {
            softness = Double.parseDouble(args[1]);
        }
        if (args.length > 1) {
            temperature = Double.parseDouble(args[2]);
        }
    	
    	final SoftSphere3d sim = new SoftSphere3d(density, softness, temperature);
        int numAtoms = 108;
        
        MeterPotentialEnergyFromIntegrator meterEnergy = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        DataPump pump = new DataPump(meterEnergy, null);
        AccumulatorAverageCollapsing accumulator = new AccumulatorAverageCollapsing();
       
        accumulator.setPushInterval(1);
        pump.setDataSink(accumulator);
        sim.integrator.addIntervalAction(pump);
        
        sim.getController().actionPerformed();
        
        
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble)((DataGroup)accumulator.getData()).getData(StatType.STANDARD_DEVIATION.index)).x;
        double energy = ((DataDouble)((DataGroup)accumulator.getData()).getData(StatType.AVERAGE.index)).x;
        Cv /= temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k: "+Cv);
        System.out.println("System Energy: "+energy/numAtoms);
        
   }

    
}