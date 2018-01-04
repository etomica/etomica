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
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.normalmode.*;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * MC simulation
 * 3D hard spheres
 * No graphic display
 * Output: histogram files of probability that a mode is zero
 * Calculate free energy of solid
 * 
 * Treats modes as degrees of freedom
 * 
 */

/*
 * starts in notes 7/09
 */
public class SimDegreeFreedom3D extends Simulation {

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimDegreeFreedom3D";
    public Primitive primitive;
    public IntegratorMC integrator;
    public BasisCubicFcc basis;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public CoordinateDefinition coordinateDefinition;
    int[] nCells;
    NormalModes nm;
    MeterNormalModeCoordinate meternmc;
    WaveVectorFactory waveVectorFactory;
    MCMoveAtomCoupled mcMoveAtom;
    MCMoveChangeMultipleWV mcMoveMode;
    AccumulatorHistogram[] hists;
    int harmonicWV;
    boolean[] skipThisMode;
    
    public SimDegreeFreedom3D(Space _space, int numAtoms, double density, int blocksize, int nbs, String filename) {
        super(_space);
        
        
        System.out.println("THIS CODE IS NOT FINISHED!");
        
//        long seed = 3;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);

        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this);
        
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
        basis = new BasisCubicFcc();
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
       
        Potential2 potential = new P2HardSphere(space, 1.0, true);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0);
        double v = primitive.unitCell().getVolume();
        primitive.scaleSize(Math.pow(v*density/4, -1.0/3.0));
        int numberOfCells = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{numberOfCells, numberOfCells, numberOfCells};
        boundary = new BoundaryDeformableLattice(primitive, nCells);
        box.setBoundary(boundary);
        
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setBox(box);
        
        nm = new NormalModesFromFile(filename, space.D());
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(1.0);
        
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
        
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        System.out.println("Number of wave vectors " + 
                waveVectorFactory.getWaveVectors().length);
        
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
            if(skipThisMode[i]) {continue;}
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
        inputFilename = params.inputfilename;
        if(filename.length() == 0){
            filename = "3DHS";
        }
        String outputName = params.outputfilename;
        double temperature = params.temperature;
        int comparedWV = params.comparedWV;
        long nSteps = params.numSteps;
        int bs = params.blockSize;
        int nbins = params.nBins;

        // construct simulation
        SimDegreeFreedom3D sim = new SimDegreeFreedom3D(Space.getInstance(D), nA, density, bs, nbins, inputFilename);
        System.out.println("Running " + APP_NAME + " "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println(nbins + " starting number of bins");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);

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
            wh = new WriteHistograms(outputName + "_" + i);
            wh.setHistogram(sim.hists[i].getHistograms());
            wh.actionPerformed();
        }
//        System.out.println("number of hits " + sim.hists[0].getHistograms().getCount());


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


        System.out.println("Fini.");
    }

    private void setHarmonicWV(int hwv) {
        harmonicWV = hwv;
        System.out.println("need to fix this setHarmonicWV");

//        mcMoveMode.setHarmonicWV(hwv);
    }
    
    public static class SimParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 1.3;
        public int D = 3;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public String inputfilename = "normal_modes_LJ_3D_32";
        public String outputfilename = "hists";
        public double temperature = 1.0;
        public int comparedWV = 7;
        public int nBins = 200;
        
        public int blockSize = 1000;
        public long numSteps = 10000;
    }

}
