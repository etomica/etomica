package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterWidomInsertion;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MC simulation of soft-sphere fluid model  
 * 
 * @author Tai Boon Tan
 */
public class SimFluidSoftSphere extends Simulation {

    public SimFluidSoftSphere(Space _space, int numAtoms, double density, double temperature, int exponent) {
        super(_space, true);


        potentialMaster = new PotentialMaster();

        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MCMoveAtom move = new MCMoveAtom(this, potentialMaster, space);
        move.setStepSize(0.2);
        move.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
        
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        // activityIntegrate.setMaxSteps(nSteps);

        double L = Math.pow(4.0/density, 1.0/3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{n,n,n};
        boundary = new BoundaryRectangularPeriodic(space, random, n * L);
        basis = new BasisCubicFcc();
        
        Potential2SoftSpherical potential = new P2SoftSphere(space);
        
        double truncationRadius = boundary.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        //potentialMaster.lrcMaster().setEnabled(false); //turn off the long-range correction ::updated 7/4/2008 
        
        IAtomTypeLeaf sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new IAtomTypeLeaf[] {sphereType, sphereType});
        //move.setPotential(pTruncated);

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
       integrator.setBox(box);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 108;
        double density = 2.2335;
        double temperature = 1.0;
        int exponent = 12 ;
        long simSteps =10000;

        // parse arguments
        if (args.length > 0) {
            density = Double.parseDouble(args[0]);
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            nA = Integer.parseInt(args[2]);
        }
        if (args.length > 3) {
            temperature = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
        	exponent = Integer.parseInt(args[4]);
        }

        System.out.println("Running fluid soft sphere simulation");
        System.out.println(nA + " atoms with exponent " + exponent+" and density "+density);
        System.out.println("isotherm temperature at "+temperature);
        System.out.println(simSteps+ " steps");

        // construct simulation
        SimFluidSoftSphere sim = new SimFluidSoftSphere(Space.getInstance(D), nA, density, temperature, exponent);
  
        MeterWidomInsertion meterInsertion = new MeterWidomInsertion(Space.getInstance(D));
        meterInsertion.setIntegrator(sim.integrator);
        meterInsertion.setSpecies(sim.species);
        //meterInsertion.setNInsert();
        
        AccumulatorAverage insertionAverage = new AccumulatorAverageCollapsing();
        DataPump insertionPump = new DataPump(meterInsertion, insertionAverage);
        sim.integrator.addIntervalAction(insertionPump);
        sim.integrator.setActionInterval(insertionPump, 1);
        
        MeterPressure meterPressure = new MeterPressure(sim.space);
        meterPressure.setIntegrator(sim.integrator);
        System.out.println("\nPressure Lattice: "+ meterPressure.getDataAsScalar());
        AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
	    DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
	    sim.integrator.addIntervalAction(pressurePump);
	    sim.integrator.setActionInterval(pressurePump, 10);
        
        
        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergy.setBox(sim.box);
        System.out.println("Lattice Energy per particle: "+ meterEnergy.getDataAsScalar()/nA);
        System.out.println(" ");
        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);
        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.setActionInterval(energyPump, 100);
        
        
        sim.activityIntegrate.setMaxSteps(simSteps/10);  //simSteps/10
        sim.getController().actionPerformed();
        System.out.println("equilibrated");
        
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().reset();
        
        sim.activityIntegrate.setMaxSteps(simSteps);
        sim.getController().actionPerformed();
        
        double insertionScalar = ((DataGroup)insertionAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
        System.out.println("Average insertion scalar: "+ insertionScalar);
        System.out.println("Error insertion scalar: "+ ((DataGroup)insertionAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        System.out.println("Liquid Gibbs free energy: " + -temperature*Math.log(insertionScalar));
        System.out.println(" ");
        
        System.out.println("Average Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index));
        System.out.println("Error Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        System.out.println(" ");
        
        System.out.println("Average Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index));
        System.out.println("Error Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public Primitive primitive;
    public Basis basis;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public PotentialMaster potentialMaster;
}