package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * The original Bennett's Overlapping Sampling Simulation
 * 	- used to check for the computation time
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapSoftSphereTP extends Simulation {

    public SimOverlapSoftSphereTP(Space _space, int numAtoms, double density, double temperature, double[] otherTemperatures, double[] alpha, int exponent, int numAlpha, double alphaSpan, long numSteps) {
        super(_space);
        
        potentialMasterTarget = new PotentialMasterList(this, space);
        
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMasterTarget, getRandom(), temperature);
        atomMove = new MCMoveAtomCoupled(potentialMasterTarget, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        double L = Math.pow(4.0/density, 1.0/3.0);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        primitive = new PrimitiveCubic(space, n*L);
        primitiveUnitCell = new PrimitiveCubic(space, L);
        
        nCells = new int[]{n,n,n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        box.setBoundary(boundary);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
     	if(potentialMasterTarget instanceof PotentialMasterList){
			potential = new P2SoftSphericalTruncated(space, potential, truncationRadius);
		
		} else {
			potential = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
			
		}
        atomMove.setPotential(potential);
        IAtomType sphereType = species.getLeafType();
        potentialMasterTarget.addPotential(potential, new IAtomType[] {sphereType, sphereType });
        
        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitiveUnitCell, box, coordinateDefinition);
        potentialMasterTarget.addPotential(p1Constraint, new IAtomType[] {sphereType});
        potentialMasterTarget.lrcMaster().setEnabled(false);
    
        integrator.setBox(box);

		if (potentialMasterTarget instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList)potentialMasterTarget).setRange(neighborRange);
            ((PotentialMasterList)potentialMasterTarget).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMasterTarget).getNeighborManager(box).reset();
            int potentialCells = ((PotentialMasterList)potentialMasterTarget).getNbrCellManager(box).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }
            if (potentialCells > cellRange*2+1) {
                System.out.println("could probably use a larger truncation radius ("+potentialCells+" > "+(cellRange*2+1)+")");
            }
            //((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundaryTarget.getBoxSize().getX(0));
		}
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy: " + latticeEnergy);
        

        meter = new MeterTargetTP(potentialMasterTarget, species, space, this);
        meter.setCoordinateDefinition(coordinateDefinition);
        meter.setLatticeEnergy(latticeEnergy);
        meter.setTemperature(temperature);
        meter.setOtherTemperatures(otherTemperatures);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps/(numBlocks*interval);
        System.out.println("block size "+blockSize+" interval "+interval);
        accumulator = new AccumulatorAverageFixed(blockSize);
        accumulatorPump = new DataPumpListener(meter, accumulator, interval);
        integrator.getEventManager().addListener(accumulatorPump);
       
        activityIntegrate = new ActivityIntegrate(integrator);
        
        getController().addAction(activityIntegrate);
    }
    
    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();
        System.out.println("equilibration finished");

        accumulator.reset();

        getController().reset();

    }
    
    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereTP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double[] otherTemperatures = params.otherTemperatures;
        double[] alpha = params.alpha;
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        
        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println((numSteps/1000)+" total steps of 1000");

        //instantiate simulation
        final SimOverlapSoftSphereTP sim = new SimOverlapSoftSphereTP(Space.getInstance(3), numMolecules, density, temperature, otherTemperatures, alpha, exponentN, numAlpha, alphaSpan, numSteps);
        
        //start simulation

        sim.initialize(numSteps/20);
        System.out.flush();
        
        final long startTime = System.currentTimeMillis();
       
        sim.activityIntegrate.setMaxSteps(numSteps);
        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        sim.getController().actionPerformed();
        //MeterTargetTP.closeFW();
        
        System.out.println("\nratio averages:\n");
        
        for (int i=0; i<otherTemperatures.length; i++) {
            System.out.println(otherTemperatures[i]);
            double[] iAlpha = sim.meter.getAlpha(i);
            for (int j=0; j<numAlpha; j++) {
                System.out.println("  "+iAlpha[j]+" "+((DataGroup)sim.accumulator.getData()).getData(AccumulatorAverage.StatType.AVERAGE.index).getValue(i*numAlpha+j)
                        +" "+((DataGroup)sim.accumulator.getData()).getData(AccumulatorAverage.StatType.ERROR.index).getValue(i*numAlpha+j)
                        +" "+((DataGroup)sim.accumulator.getData()).getData(AccumulatorAverage.StatType.BLOCK_CORRELATION.index).getValue(i*numAlpha+j));
            }
        }
        
        long endTime = System.currentTimeMillis();
        System.out.println("Time taken: " + (endTime - startTime)/1000.0);
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive, primitiveUnitCell;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterTargetTP meter;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMaster potentialMasterTarget;
    protected double latticeEnergy;
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 500;
        public double density = 1.256;
        public int exponentN = 12;
        public long numSteps = 100000000;
        public double temperature = 1.2;
        public double[] alpha = new double[]{0.011};
        public int numAlpha = 10;
        public double alphaSpan = 1;
        public double[] otherTemperatures = new double[]{1.1};
    }
}
