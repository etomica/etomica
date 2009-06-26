package etomica.models.oneDHardRods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.box.Box;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceScalar;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MeterBoltzmannTarget;
import etomica.normalmode.MeterHarmonicEnergy;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModesFromFile;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WriteS;
import etomica.normalmode.SimOverlapLJModule.MeterPotentialEnergyDifference;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Energy;
import etomica.units.Null;
import etomica.util.ParameterBase;
import etomica.util.RandomNumberGenerator;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * MC simulation of Lennard-Jones system.
 * @author cribbin
 */
public class TestMCMoveChangeSingleMode3DLJ extends Simulation {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimSingleWaveVector";
    Primitive primitive;
    int[] nCells;
    NormalModes nm;
    public BasisMonatomic basis;
    ActivityIntegrate activityIntegrate;
    
    IntegratorMC integrator;
    
    public IBox box;
    public Boundary boundary;
    MCMoveChangeSingleMode changeMove;

    public TestMCMoveChangeSingleMode3DLJ(Space _space, int numAtoms, double density, double 
            temperature, String filename, double harmonicFudge, int awv){
        super(_space, true);
        
//        long seed = 5;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);
        
        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        
        basis = new BasisMonatomic(space);
        
//TARGET
        //Set up target system   - A - 1
        PotentialMasterMonatomic potentialMaster = new 
                PotentialMasterMonatomic(this);
        integrator = new IntegratorMC(this, potentialMaster);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        primitive = new PrimitiveCubic(space, 1.0);
        double v = primitive.unitCell().getVolume();
        primitive.scaleSize(Math.pow(v*density/4, -1.0/3.0));
        int numberOfCells = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{numberOfCells, numberOfCells, numberOfCells};
        boundary = new BoundaryDeformableLattice(primitive, nCells);
        box.setBoundary(boundary);
        Basis basis = new BasisCubicFcc();
        
        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
//        for(int k = 0; k < 32; k++){
//            System.out.println(k + " " +((IAtomPositioned)coordinateDefinitionTarget.getBox().getLeafList().getAtom(k)).getPosition());
//        }
        
        
        
        Potential2SoftSpherical p2 = new P2LennardJones(space, 1.0, 1.0);
        double truncationRadius = boundary.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new 
                P2SoftSphericalTruncatedShifted(space, p2, truncationRadius);
        potentialMaster.addPotential(pTruncated, new IAtomType[]
                {species.getLeafType(), species.getLeafType()});
        
        nm = new NormalModesFromFile(filename, space.D());
        
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        nm.getOmegaSquared(box);
        
        WaveVectorFactory waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        int wvflength = waveVectorFactory.getWaveVectors().length;
        System.out.println("We have " + wvflength +" wave vectors.");
        System.out.println("Wave Vector Coefficients:");
        for(int i = 0; i < wvflength; i++){
            System.out.println(i + " " + waveVectorFactory.getCoefficients()[i]);
        }
        
        changeMove = new MCMoveChangeSingleMode(potentialMaster, random);
        integrator.getMoveManager().addMCMove(changeMove);
        changeMove.setWaveVectors(waveVectorFactory.getWaveVectors());
        changeMove.setWaveVectorCoefficients(
                waveVectorFactory.getCoefficients());
        changeMove.setEigenVectors(nm.getEigenvectors(box));
        changeMove.setCoordinateDefinition(coordinateDefinition);
        changeMove.setBox((IBox)box);
        changeMove.setStepSizeMin(0.001);
        changeMove.setStepSize(0.01);

        
        
        
        
//JOINT
        //Set up the rest of the joint stuff
        setComparedWV(awv);
       
        
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
        
    }

    
    public static void main(String args[]){
        SimOverlapSingleWaveVector3DParam params = new SimOverlapSingleWaveVector3DParam();
        String inputFilename = null;
        if(args.length > 0){
            inputFilename = args[0];
        }
        if(inputFilename != null){
            ReadParameters readParameters = new 
                ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        
        int numMolecules = params.numAtoms;
        double density = params.density;
        int D = params.D;
        double harmonicFudge = params.harmonicFudge;
        String filename = params.filename;
        if(filename.length() == 0){
            filename = "1DHR";
        }
        double temperature = params.temperature;
        int comparedWV = params.comparedWV;
        
        String refFileName = args.length > 0 ? filename+"_ref" : null;
        
        //instantiate simulations!
        TestMCMoveChangeSingleMode3DLJ sim = new TestMCMoveChangeSingleMode3DLJ  (Space.getInstance(D), numMolecules,
                density, temperature, filename, harmonicFudge, comparedWV);
        int numSteps = params.numSteps;
        
        System.out.println("Running Nancy's single " +D+"D Lennard Jones simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println("temperature: " + temperature);
        System.out.println("compared wave vector: " + comparedWV);
        System.out.println("Total steps: "+params.numSteps);
        System.out.println("instantiated");
        
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();


    }
    
    
    public void setComparedWV(int awv){
        changeMove.setHarmonicWV(awv);
    }
    
    public static class SimOverlapSingleWaveVector3DParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 1.3;
        public int D = 3;
        public double harmonicFudge = 1.0;
        public double temperature = 1.0;
        public int comparedWV = 7;
        
        public int numSteps = 40000000;
        
        public String filename = "testMCMoveChange3DLJ";


    }
 
}
