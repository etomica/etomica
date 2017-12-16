/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPump;
import etomica.data.DataSplitter;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramSimple;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
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
 * MC simulation
 * 1D hard rods
 * No graphic display
 * Output: histogram files of probability that a mode is zero
 * Calculate free energy of solid
 * 
 * Testing whether or not using modes as degrees of freedom is extensive.
 */

/* 
 * starts in notes 8/09
 */
public class SimExtensiveTest extends Simulation {

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimExtensiveTest";
    public Primitive primitive;
    public IntegratorMC integrator;
    public BasisMonatomic basis;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary bdry;
    public CoordinateDefinition coordinateDefinition;
    int[] nCells;
    NormalModes nm;
    MeterNormalModeCoordinate meternmc;
    WaveVectorFactory waveVectorFactory;
//    MCMoveAtomCoupled mcMoveAtom;
    MCMoveChangeMultipleWV mcMoveMode;
    AccumulatorHistogram[] hists;
    int[] harmonicWV;
    boolean[] skipThisMode;

    public SimExtensiveTest(Space _space, int numAtoms, double density, int blocksize, int nbs) {
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
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
       
        Potential2 potential = new P2HardSphere(space, 1.0, true);
        potential = new P2XOrder(space, (Potential2HardSpherical)potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        primitive = new PrimitiveCubic(space, 1.0/density);
        bdry = new BoundaryRectangularPeriodic(space, numAtoms/density);
        nCells = new int[]{numAtoms};
        box.setBoundary(bdry);
        
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        int coordinateDim = coordinateDefinition.getCoordinateDim();

        double neighborRange = 1.01/density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setBox(box);
        
        nm = new NormalModes1DHR(box.getBoundary(), numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        nm.getOmegaSquared();
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        //Set up skip-these-modes code
        double[] wvc= nm.getWaveVectorFactory().getCoefficients();
        double[][] omega = nm.getOmegaSquared();
        int jump = coordinateDim * nm.getWaveVectorFactory().getWaveVectors().length;
        skipThisMode = new boolean[2*jump];
        for(int i = 0; i < 2*jump; i++){
            skipThisMode[i] = false;
        }
        for(int wvCount = 0; wvCount < wvc.length; wvCount++){
            //Sets up the imaginary modes that should be skipped.
            if(wvc[wvCount] == 0.5) {
                for(int j = 0; j < coordinateDim; j++){
                    skipThisMode[j + coordinateDim*wvCount + jump] = true;
                }
            }
            //Sets up the modes that are center of mass motion to skip
            for(int j = 0; j < omega[wvCount].length; j++){
                if(Double.isInfinite(omega[wvCount][j])){
                    skipThisMode[j + coordinateDim*wvCount] = true;
                    skipThisMode[j + coordinateDim*wvCount + jump] = true;
                }
            }
        }
        
        
//        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, random, space);
//        mcMoveAtom.setPotential(potential);
//        mcMoveAtom.setBox(box);
//        integrator.getMoveManager().addMCMove(mcMoveAtom);
//        mcMoveAtom.setStepSizeMin(0.001);
//        mcMoveAtom.setStepSize(0.01);
        
        mcMoveMode = new MCMoveChangeMultipleWV(potentialMaster, random);
        mcMoveMode.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        mcMoveMode.setCoordinateDefinition(coordinateDefinition);
        mcMoveMode.setEigenVectors(nm.getEigenvectors());
        mcMoveMode.setOmegaSquared(nm.getOmegaSquared());
        mcMoveMode.setWaveVectorCoefficients(nm.getWaveVectorFactory().getCoefficients());
        mcMoveMode.setWaveVectors(nm.getWaveVectorFactory().getWaveVectors());
        
        meternmc = new MeterNormalModeCoordinate(coordinateDefinition, nm.getWaveVectorFactory().getWaveVectors());
        meternmc.setEigenVectors(nm.getEigenvectors());
        meternmc.setOmegaSquared(nm.getOmegaSquared());
        
        int coordNum = nm.getWaveVectorFactory().getWaveVectors().length*coordinateDim*2;
        hists = new AccumulatorHistogram[coordNum];
        DataSplitter splitter = new DataSplitter();
        DataPump pumpFromMeter = new DataPump(meternmc, splitter);
        
        DoubleRange range = new DoubleRange(-1.0, 1.0);
        Histogram template;
        for(int i = 0; i < coordNum; i++){
            template = new HistogramSimple(nbs, range);
            hists[i] = new AccumulatorHistogram(template, nbs);
            splitter.setDataSink(i, hists[i]);
        }
        
        IntegratorListenerAction pumpFromMeterListener = new IntegratorListenerAction(pumpFromMeter);
        pumpFromMeterListener.setInterval(blocksize);
        integrator.getEventManager().addListener(pumpFromMeterListener);
        
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);
        
        
//        IAtomList leaflist = box.getLeafList();
//        double[] locations = new double[numAtoms];
//        System.out.println("starting positions:");
//        for(int i = 0; i < numAtoms; i++){
//            //one d is assumed here.
//            locations[i] = ( ((Atom)leaflist.getAtom(i)).getPosition().x(0) );
//        }
//        
//        for(int i = 0; i < numAtoms; i++){
//            System.out.println(i + "  " + locations[i]);
//        }
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
        int[] comparedWV = params.comparedWV;
        long nSteps = params.numSteps;
        int bs = params.blockSize;
        int nbins = params.nBins;
        String outputfn = params.outputname;


        // construct simulation
        SimExtensiveTest sim = new SimExtensiveTest(Space.getInstance(D), nA, density, bs, nbins);
        System.out.println("Running " + APP_NAME + " "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println(nbins + " starting number of bins");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);
        System.out.println("instantiated");

        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.setHarmonicWV(comparedWV);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();

        int accumulatorLength = sim.hists.length;
        for(int i = 0; i < accumulatorLength; i++){
            if(sim.skipThisMode[i]) {continue;}
            sim.hists[i].reset();
        }

        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();

        /*
         * This loop creates a new write class for each histogram from each
         * AccumulatorHistogram, changes the filename for the histogram output,
         * connects the write class with this histogram, and
         * writes out the results to the file.
         */
        WriteHistograms wh;
        for(int i = 0; i < accumulatorLength; i++){
            if(sim.skipThisMode[i]) {continue;}
            String outputName = new String(outputfn + "_" + i);
            wh = new WriteHistograms(outputName);
            wh.setHistogram(sim.hists[i].getHistograms());
            wh.actionPerformed();
        }

//        IAtomList leaflist = sim.box.getLeafList();
//        double[] locations = new double[nA];
//        System.out.println("final:");
//        for(int i = 0; i < nA; i++){
//            //one d is assumed here.
//            locations[i] = ( ((Atom)leaflist.getAtom(i)).getPosition().x(0) );
//        }
//
//        for(int i = 0; i < 32; i++){
//            System.out.println(i + "  " + locations[i]);
//        }

        if(D==1) {
            double AHR = -(nA-1)*Math.log(nA/density-nA)
                + SpecialFunctions.lnFactorial(nA) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }

        System.out.println("Fini.");
    }

    private void setHarmonicWV(int[] hwv) {
        harmonicWV = hwv;
//        mcMoveMode.setHarmonicWV(hwv);
        System.out.println("need to fix this setHarmonicWV");
    }
    
    public static class SimParam extends ParameterBase {
        public int numAtoms = 500;
        public double density = 0.70;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public String inputfilename = "input";
        public String outputname = "hists";
        public double temperature = 1.0;
        public int[] comparedWV = {499, 497, 495, 493, 491, 489, 487, 485, 483, 
                481, 479, 477, 475, 473, 471, 469, 467, 465, 463, 461, 459, 457,
                455, 453, 451, 449, 447, 445, 443, 441, 439, 437, 435, 433, 431,
                429, 427, 425, 423, 421, 419, 417, 415, 413, 411, 409, 407, 405, 
                403, 401, 399, 397, 395, 393, 391, 389, 387, 385, 383, 381, 379,
                377, 375, 373, 371, 369, 367, 365, 363, 361, 359, 357, 355, 353,
                351, 349, 347, 345, 343, 341, 339, 337, 335, 333, 331, 329, 327,
                325, 323, 321, 319, 317, 315, 313, 311, 309, 307, 305, 303, 301,
                299, 297, 295, 293, 291, 289, 287, 285, 283, 281, 279, 277, 275,
                273, 271, 269, 267, 265, 263, 261, 259, 257, 255, 253, 251, 249,
                247, 245, 243, 241, 239, 237, 235, 233, 231, 229, 227, 225, 223,
                221, 219, 217, 215, 213, 211, 209, 207, 205, 203, 201, 199, 197,
                195, 193, 191, 189, 187, 185, 183, 181, 179, 177, 175, 173, 171,
                169, 167, 165, 163, 161, 159, 157, 155, 153, 151, 149, 147, 145,
                143, 141, 139, 137, 135, 133, 131, 129, 127, 125, 123, 121, 119,
                117, 115, 113, 111, 109, 107, 105, 103, 101, 99, 97, 95, 93, 91,
                89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 
                57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 
                25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1 };
        public int nBins = 200;
        
        public int blockSize = 1000;
        public long numSteps = 100000;
    }

}
