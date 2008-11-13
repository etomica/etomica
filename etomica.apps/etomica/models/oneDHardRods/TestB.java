package etomica.models.oneDHardRods;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.atom.AtomLeaf;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WriteS;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

public class TestB extends Simulation {
    private static final long serialVersionUID = 1L;
    public Boundary boundary;
    IntegratorMC integrator;
    ActivityIntegrate activityIntegrate;
    IBox box;
    CoordinateDefinition coordinateDefinition;
    Primitive primitive;
    int[] nCells;
    SpeciesSpheresMono species;
    NormalModes1DHR nm;
    double[] locations;
    int affectedWV;
    AccumulatorAverageFixed avgOverlap;
    
    private static final String APP_NAME = "TestB";
    

    public TestB(Space _space, int numAtoms, double density, double temperature,
            String filename, double harmonicFudge, int awv){
        super(_space, true);
        
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        p2.setBox(box);
        potentialMaster.addPotential(p2, new IAtomTypeLeaf[]
                {species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundary = new BoundaryRectangularPeriodic(space, getRandom(),
                numAtoms/density);
        nCells = new int[]{numAtoms};
        box.setBoundary(boundary);
        
        coordinateDefinition = new 
                CoordinateDefinitionLeaf(this, box, primitive, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        double neighborRange = 1.01/density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();
        
        integrator = new IntegratorMC(potentialMaster, random, temperature);
        integrator.setBox(box);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        affectedWV = awv;
        
        MCMoveCompareMode convert = new MCMoveCompareMode(potentialMaster, random);
        integrator.getMoveManager().addMCMove(convert);
        convert.setWaveVectors(waveVectorFactory.getWaveVectors());
        convert.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        convert.setOmegaSquared(nm.getOmegaSquared(box), waveVectorFactory.getCoefficients());
        convert.setEigenVectors(nm.getEigenvectors(box));
        convert.setCoordinateDefinition(coordinateDefinition);
        convert.setTemperature(temperature);
        convert.setBox((IBox)box);
        convert.setStepSizeMin(0.001);
        convert.setStepSize(0.01);
        convert.setConvertedWaveVector(affectedWV);
        
        integrator.setBox(box);
        potentialMaster.getNeighborManager(box).reset();
        
        locations = new double[numAtoms];
        IAtomSet leaflist = box.getLeafList();
        for(int i = 0; i < numAtoms; i++){
            //one d is assumed here.
            locations[i] = ( ((AtomLeaf)leaflist.getAtom(i)).getPosition().x(0) );
        }
        
        MeterPotentialEnergy meterAinB = new MeterPotentialEnergy(potentialMaster);
        meterAinB.setBox(box);
       
        MeterConvertModeBrute meterBinB = new MeterConvertModeBrute(potentialMaster,coordinateDefinition,box);
        meterBinB.setCoordinateDefinition(coordinateDefinition);
        meterBinB.setEigenVectors(nm.getEigenvectors(box));
        meterBinB.setOmegaSquared(nm.getOmegaSquared(box));
        meterBinB.setTemperature(temperature);
        meterBinB.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        meterBinB.setWaveVectors(waveVectorFactory.getWaveVectors());
        meterBinB.setConvertedWV(affectedWV);
        
        MeterOverlapTestB meterOverlapB = new MeterOverlapTestB(meterAinB, meterBinB, temperature);
        avgOverlap = new AccumulatorAverageFixed();
        DataPump pumpOverlap = new DataPump(meterOverlapB, avgOverlap);
        integrator.addIntervalAction(pumpOverlap);
        integrator.setActionInterval(pumpOverlap, 100);
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        /*
         * This whole setup defines a set of default parameters
         * in the inner class Sim1DHRParams.  These parameters can be changed
         * individually in an appropriately named file, without affecting
         * the values of the other parameters.  The order of definition in the
         * file is irrelevant.
         * 
         */
        Sim1DHRParams params = new Sim1DHRParams();
        String inputFilename = null;
        if(args.length > 0){
            inputFilename = args[0];
        }
        if(inputFilename != null){
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density;
        long numSteps = params.numSteps;
        int numAtoms = params.numAtoms;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if(filename.length() == 0){
            filename = "normal_modes_1DHR _" + numAtoms;
        }
        String refFileName = args.length>0 ? filename+"_ref" : null;
        
        System.out.println("Running 1D hard rod simulation - TestB");
        System.out.println(numAtoms+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+ " total steps of 1000");
        System.out.println("output data to "+filename);
        System.out.println("affected wave vector = " + params.aWV);
        
        
        //instantiate simulation
        TestB sim = new TestB(Space.getInstance(D), numAtoms, density, 
                temperature, filename, harmonicFudge, params.aWV);
        sim.activityIntegrate.setMaxSteps(numSteps);
        
        MeterNormalMode mnm = new MeterNormalMode();
        mnm.setCoordinateDefinition(sim.coordinateDefinition);
        mnm.setWaveVectorFactory(sim.nm.getWaveVectorFactory());
        mnm.setBox(sim.box);
        mnm.reset();
        
        sim.integrator.addIntervalAction(mnm);
        sim.integrator.setActionInterval(mnm, 2);
        
        
        
        
        
        
        
        
        ((Controller)sim.getController()).actionPerformed();

        
        
        //print out final positions:
        try {
            FileWriter fileWriterE = new FileWriter(filename+".txt");
            for (int i = 0; i<numAtoms; i++) {
                fileWriterE.write(Double.toString(sim.locations[i]));
                fileWriterE.write("\n");
            }
            fileWriterE.write("\n");
            fileWriterE.close();
        }
        catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data "+e);
        }
        

        WriteS sWriter = new WriteS(sim.space);
        sWriter.setFilename(filename);
        sWriter.setMeter(mnm);
        sWriter.setWaveVectorFactory(sim.nm.getWaveVectorFactory());
        sWriter.setOverwrite(true);
        sWriter.actionPerformed();
        
        
        
        System.out.println("overlap: " + sim.avgOverlap.getData().getValue(StatType.AVERAGE.index));
        System.out.println("error: " + sim.avgOverlap.getData().getValue(StatType.ERROR.index));
        System.out.println("Fini.");
    }

    
    /**
     * Inner class for parameters understood by the class's constructor
     */
    public static class Sim1DHRParams extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.5;
        public int D = 1;
        public long numSteps = 1000000;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
        public int aWV = 16;
    }
}
