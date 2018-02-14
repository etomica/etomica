/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
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
 * MC simulation
 * 
 * 1D hard rods
 * No graphic display
 * Output: histogram files of probability that a mode is zero
 * Calculate free energy of solid
 * 
 * Treats modes as degrees of freedom; uses a MeterDifferentImage to calculate 
 * what happens when an extra mode is added.
 * 
 */

/*
 * Starts in notes 7/09
 */
public class TestDifferentImage1DHRAdd extends Simulation {

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimDegreeFreedom1DHR";
    public Primitive primitive;
    public IntegratorMC integrator;
    public BasisMonatomic basis;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary bdry;
    public CoordinateDefinition coordinateDefinition;
    int[] nCells;
    NormalModes nm;
    MeterDifferentImageAdd1D meterdi;
    WaveVectorFactory waveVectorFactory;
    MCMoveAtomCoupled mcMoveAtom;
    MCMoveChangeMultipleWV mcMoveMode;
    AccumulatorHistogram[] hists;
    int harmonicWV;
    boolean[] skipThisMode;
    AccumulatorAverageFixed accumulatorDI;


    public TestDifferentImage1DHRAdd(Space _space, int numAtoms, double density, 
            int blocksize, int[] changeable) {
        super(_space);

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
        potential = new P2XOrder(space, (Potential2HardSpherical) potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        primitive = new PrimitiveCubic(space, 1.0 / density);
        bdry = new BoundaryRectangularPeriodic(space, numAtoms / density);
        nCells = new int[]{numAtoms};
        box.setBoundary(bdry);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        int coordinateDim = coordinateDefinition.getCoordinateDim();

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
        mcMoveMode.addChangeableWV(changeable);

        meterdi = new MeterDifferentImageAdd1D("MeterDI", /*potentialMaster,*/ numAtoms, density, this,
                primitive, basis, coordinateDefinition, nm, 1.0);

        accumulatorDI = new AccumulatorAverageFixed(blocksize);
        DataPump pumpFromMeter = new DataPump(meterdi, accumulatorDI);

        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pumpFromMeter);
        pumpListener.setInterval(blocksize);
        integrator.getEventManager().addListener(pumpListener);

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
        int[] changeableWV = params.changeableWV;
        long nSteps = params.numSteps;
        int bs = params.blockSize;
        String outputfn = params.outputname;

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " steps, " + bs + " blocksize");
        System.out.println("input data from " + inputFilename);
        System.out.println("output data to " + filename);

        // construct simulation
        TestDifferentImage1DHRAdd sim = new TestDifferentImage1DHRAdd(Space.getInstance(D),
                nA, density, bs, changeableWV);

        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();

        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();

        //After processing...
        DataGroup group = (DataGroup)sim.accumulatorDI.getData();
        double results = ((DataDouble) group.getData(AccumulatorAverage.AVERAGE.index)).x;
        System.out.println("results: " + results);

        if(D==1) {
            double AHR = -(nA-1)*Math.log(nA/density-nA)
                + SpecialFunctions.lnFactorial(nA) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }

        System.out.println("Fini.");
    }

    public Primitive getPrimitive() {
        return primitive;
    }

    public BasisMonatomic getBasis() {
        return basis;
    }
    
    public static class SimParam extends ParameterBase {
        public int numAtoms = 10;
        public double density = 0.70;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public String inputfilename = "input";
        public String outputname = "hists";
        public double temperature = 1.0;
        public int[] changeableWV = {0, 1, 2, 3, 4};
        public int nBins = 200;
        
        public int blockSize = 1000;
        public long numSteps = 1000000000000L;
    }

}
