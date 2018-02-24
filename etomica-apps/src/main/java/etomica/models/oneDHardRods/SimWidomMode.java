/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.*;
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

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimWidomMode extends Simulation {

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimWidomMode";
    public Primitive primitive;
    public IntegratorMC integrator;
    public BasisMonatomic basis;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary bdry;
    public CoordinateDefinition coordinateDefinition;
    int[] nCells;
    NormalModes nm;
    WaveVectorFactory waveVectorFactory;
    MCMoveAtomCoupled mcMoveAtom;
    MCMoveChangeMultipleWV mcMoveMode;
    int harmonicWV;
    MeterWidomModeReal[] realMeter;
    MeterWidomModeImaginary[] imagMeter;
    AccumulatorAverage[] accumulators;

    public SimWidomMode(Space _space, int numAtoms, double density, int blocksize) {
        super(_space);

        System.out.println("THIS CODE IS NOT FINISHED!");
        System.out.println("need to fix this setHarmonicWV");


//        long seed = 3;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        basis = new BasisMonatomic(space);
        bdry = new BoundaryRectangularPeriodic(space, numAtoms / density);
        box = this.makeBox(bdry);
        box.setNMolecules(species, numAtoms);

        Potential2 potential = new P2HardSphere(space, 1.0, true);
        potential = new P2XOrder(space, (Potential2HardSpherical) potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        primitive = new PrimitiveCubic(space, 1.0 / density);
        nCells = new int[]{numAtoms};

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);

        double neighborRange = 1.01 / density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();

        integrator = new IntegratorMC(this, potentialMaster, box);

        nm = new NormalModes1DHR(box.getBoundary(), numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared();

        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);

        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, new MeterPotentialEnergy(potentialMaster), random, space);
        mcMoveAtom.setPotential(potential);
        mcMoveAtom.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtom.setStepSizeMin(0.001);
        mcMoveAtom.setStepSize(0.01);

        mcMoveMode = new MCMoveChangeMultipleWV(potentialMaster, random);
        mcMoveMode.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        mcMoveMode.setCoordinateDefinition(coordinateDefinition);
        mcMoveMode.setEigenVectors(nm.getEigenvectors());
        mcMoveMode.setOmegaSquared(nm.getOmegaSquared());
        mcMoveMode.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        mcMoveMode.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int coordNum = nm.getWaveVectorFactory().getWaveVectors().length * coordinateDim;

        realMeter = new MeterWidomModeReal[coordNum];
        accumulators = new AccumulatorAverageFixed[coordNum * 2];
        DataPump pump;
        IntegratorListenerAction pumpListener;
        for (int i = 0; i < coordNum; i++) {
            String name = new String("widom Meter for mode " + i);
            realMeter[i] = new MeterWidomModeReal(name, potentialMaster,
                    coordinateDefinition, box, i);
            realMeter[i].setEigenVectors(nm.getEigenvectors());
            realMeter[i].setOmegaSquared(nm.getOmegaSquared());
            realMeter[i].setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
            realMeter[i].setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());

            accumulators[i] = new AccumulatorAverageFixed(blocksize);

            pump = new DataPump(realMeter[i], accumulators[i]);
            pumpListener = new IntegratorListenerAction(pump);
            pumpListener.setInterval(blocksize);
            integrator.getEventManager().addListener(pumpListener);
        }

        imagMeter = new MeterWidomModeImaginary[coordNum];
        for (int i = 0; i < coordNum; i++) {
            String name = new String("widom Meter for mode " + i);
            imagMeter[i] = new MeterWidomModeImaginary(name, potentialMaster,
                    coordinateDefinition, box, i);
            imagMeter[i].setEigenVectors(nm.getEigenvectors());
            imagMeter[i].setOmegaSquared(nm.getOmegaSquared());
            imagMeter[i].setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
            imagMeter[i].setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());

            accumulators[i + coordNum] = new AccumulatorAverageFixed(blocksize);

            pump = new DataPump(imagMeter[i], accumulators[i + coordNum]);
            pumpListener = new IntegratorListenerAction(pump);
            pumpListener.setInterval(blocksize);
            integrator.getEventManager().addListener(pumpListener);
        }

        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        SimParam params = new SimParam();
        String inputFilename = null;
        if(args.length > 0) {
            inputFilename = args[0];
        }
        if(inputFilename != null){
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
            inputFilename = params.inputfilename;
        }

        int nA = params.numAtoms;
        double density = params.density;
        int D = params.D;
        double harmonicFudge = params.harmonicFudge;
        String filename = params.filename;
        if(filename.length() == 0){
            filename = "1DHR";
        }
        double temperature = params.temperature;
        int comparedWV = params.comparedWV;
        int nSteps = params.numSteps;
        int bs = params.blockSize;


        SimWidomMode sim = new SimWidomMode(Space.getInstance(D), nA, density, bs);
        System.out.println("Running " + APP_NAME);
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);

        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.setHarmonicWV(comparedWV);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();

        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();

        //After processing...
        int cd = sim.nm.getWaveVectorFactory().getWaveVectors().length *
                sim.coordinateDefinition.getCoordinateDim() * 2;
        double[] results = new double[cd];
        DataGroup group;
        for(int i = 0; i < cd; i++){
            group = (DataGroup)sim.accumulators[i].getData();
            results[i] = ((DataDouble) group.getData(AccumulatorAverage.AVERAGE.index)).x;
        }
        for(int i = 0; i < cd; i++){
            System.out.println(i + "  " + results[i]);
        }

        if(D==1) {
            double AHR = -(nA-1)*Math.log(nA/density-nA)
                + SpecialFunctions.lnFactorial(nA) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }

        System.out.println("Fini.");
    }

    private void setHarmonicWV(int hwv) {
        harmonicWV = hwv;
        System.out.println("THIS CODE IS NOT FINISHED!");
        System.out.println("need to fix this setHarmonicWV");
//        mcMoveMode.setHarmonicWV(hwv);
    }
    
    public static class SimParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.50;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public String inputfilename = "input";
        public double temperature = 1.0;
        public int comparedWV = numAtoms/2;
        
        public int blockSize = 1000;
        public int numSteps = 10000;
    }

}
