/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorScheme;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.*;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

import java.awt.*;
import java.io.*;

/**
 * MC simulation
 * 3D Lennard Jones
 * FCC crystal
 * No graphic display
 * uses single big cell approach
 * Calculate free energy of solid using normal mode insertion method
 * 
 * Treats coupling of modes as overlap variable; 
 * 
 * Generate input files with HarmonicCrystalSoftSphereFCC
 * 
 * Uses overlap sampling.
 */

/*
 * Starts in notes 5/21/10
 */
public class SimDifferentImageSsFccDoubleSize extends Simulation {

    private static final String APP_NAME = "SimDifferentImageFCC";
    public Primitive primitive;
    public Basis basis;
    public ActivityIntegrate activityIntegrate;
    public CoordinateDefinition cDefTarget, cDefRef;
    public IntegratorOverlap integratorSim; //integrator for the whole simulation
    public DataSourceVirialOverlap dsvo;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public Box boxTarget, boxRef;
    public Boundary bdryTarget, bdryRef;
    NormalModes nmRef, nmTarg;
    WaveVectorFactory waveVectorFactoryRef, waveVectorFactoryTarg;
    double bennettParam;       //adjustable parameter - Bennett's parameter
    IntegratorMC[] integrators;
    MeterPotentialEnergy meterTargInTarg, meterRef, meterRefInRef;
    MeterDifferentImageAddDoubleSize meterTargInRef;
    MeterDifferentImageSubtractDoubleSize meterRefInTarg;
    
    double refSumWVC, targSumWVC;
    int targModesCt, refModeCt;
    double constraint; 
    
    
    public SimDifferentImageSsFccDoubleSize(Space _space, int[] nCellsRef, 
            int[] nCellsTarget, double density, double tems, int exponent, String inputFile, double constraint) {
        super(_space);
        System.out.println("Running " + APP_NAME);
        
//        long seed = 0;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);
        
        int targAtoms = 1;
        int refAtoms = 1;
        for(int i = 0; i < space.D(); i++){
            refAtoms *= nCellsRef[i];
            targAtoms *= nCellsTarget[i];
        }
        refAtoms *= 4;     //definitely fcc
        targAtoms *= 4;    //definitely fcc
        
        double temperature = tems;
        this.constraint = constraint;
        String rIn = inputFile + refAtoms;
        String tIn = inputFile + targAtoms;
        
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
        integrators = new IntegratorMC[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
//REFERENCE
        // Set up reference system - B, 0
        boxRef = new Box(space);
        addBox(boxRef);
        boxRef.setNMolecules(species, refAtoms);
        
        double primitiveLength = Math.pow(4.0 / density, 1.0 / 3.0);
        bdryRef = new BoundaryRectangularPeriodic(space, 1.0);
        Vector edges = new Vector3D();
        double[] lengths = new double[3];
        lengths[0] = nCellsRef[0]*primitiveLength;
        lengths[1] = nCellsRef[1]*primitiveLength;
        lengths[2] = nCellsRef[2]*primitiveLength;
        edges.E(lengths);
        bdryRef.setBoxSize(edges);
        boxRef.setBoundary(bdryRef);
        primitive = new PrimitiveOrthorhombic(space, lengths[0], lengths[1],
                lengths[2]);
        
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCellsRef);
        
        cDefRef = new CoordinateDefinitionLeaf(boxRef, primitive, basis, space);
        cDefRef.initializeCoordinates(new int[] {1, 1, 1});
        
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        //Choose the smallest side to define the neighborRange.
        double neighborRange = 0.0;
        if(nCellsRef[0] <= nCellsRef[1] && nCellsRef[0] <= nCellsRef[2]){
            neighborRange = 0.495 * lengths[0];
        } else if(nCellsRef[1] <= nCellsRef[2]) {
            neighborRange = 0.495 * lengths[1];
        }else {
            neighborRange = 0.495 * lengths[2];
        }
        
        if(refAtoms >= 256){
            neighborRange = 2.2;
        }
        
        System.out.println("truncation " + neighborRange);
        Potential2SoftSpherical potentialBase = new P2SoftSphere(space, 1.0, 
                1.0, exponent);
        P2SoftSphericalTruncated potential = new P2SoftSphericalTruncated(
                space, potentialBase, neighborRange);
        potentialMaster.addPotential(potential, new AtomType[]{
                species.getLeafType(), species.getLeafType()});
        potentialMaster.setRange(neighborRange);
        potentialMaster.lrcMaster().setEnabled(false);
        potentialMaster.getNeighborManager(boxRef).reset();
        
        IntegratorMC integratorRef = new IntegratorMC(potentialMaster, random, temperature);
        integratorRef.setBox(boxRef);
        integrators[0] = integratorRef;
        
        nmRef = new NormalModesFromFile(rIn, space.D());
        nmRef.setHarmonicFudge(1.0);
//        nmRef.setTemperature(temperature);  //not needed - deriv based
//        double[][] omega = nmRef.getOmegaSquared();
        waveVectorFactoryRef = nmRef.getWaveVectorFactory();
        waveVectorFactoryRef.makeWaveVectors(boxRef);
//        double[] wvc= nmRef.getWaveVectorFactory().getCoefficients();
        
        System.out.println("We have " + waveVectorFactoryRef.getWaveVectors().length
                +" reference wave vectors.");
        
        meterRefInRef = new MeterPotentialEnergy(potentialMaster);
        meterRefInRef.setBox(boxRef);
        double latticeEnergyRef = meterRefInRef.getDataAsScalar();
        System.out.println("Reference system lattice energy: " +latticeEnergyRef);
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(boxRef);
        MCMoveAtomCoupled mcMoveAtom = new MCMoveAtomCoupled(potentialMaster,
                meterPE, random, space);
        mcMoveAtom.setPotential(potential);
        mcMoveAtom.setDoExcludeNonNeighbors(true);
        mcMoveAtom.setStepSize(0.01);
        integratorRef.getMoveManager().addMCMove(mcMoveAtom);
        integratorRef.setMeterPotentialEnergy(meterRefInRef);
        
        
//TARGET
        // Set up target system
        boxTarget = new Box(space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, targAtoms);
        
        bdryTarget = new BoundaryRectangularPeriodic(space, 1.0);
        edges = new Vector3D();
        lengths = new double[3];
        lengths[0] = nCellsTarget[0]*primitiveLength;
        lengths[1] = nCellsTarget[1]*primitiveLength;
        lengths[2] = nCellsTarget[2]*primitiveLength;
        edges.E(lengths);
        bdryTarget.setBoxSize(edges);
        boxTarget.setBoundary(bdryTarget);
        
        cDefTarget = new CoordinateDefinitionLeaf(boxTarget, primitive, basis, space);
        int[] size = new int[space.D()];
        for(int i=0; i < space.D(); i++){
            size[i] = nCellsTarget[i] / nCellsRef[i];
        }
        cDefTarget.initializeCoordinates(size);
        
        potentialMaster.getNeighborManager(boxTarget).reset();
        
        IntegratorMC integratorTarget = new IntegratorMC(potentialMaster,
                random, temperature);
        integrators[1] = integratorTarget;
        integratorTarget.setBox(boxTarget);
        
        nmTarg = new NormalModesFromFile(tIn, space.D());
        nmTarg.setHarmonicFudge(1.0);
//        nmTarg.setTemperature(temperature);  // notneeded, deriv based
//        omega = nmTarg.getOmegaSquared();
        waveVectorFactoryTarg = nmTarg.getWaveVectorFactory();
        waveVectorFactoryTarg.makeWaveVectors(boxTarget);
//        wvc = nmTarg.getWaveVectorFactory().getCoefficients();
        
        System.out.println("We have " + waveVectorFactoryTarg.getWaveVectors().length 
                +" target wave vectors.");
        
        meterTargInTarg = new MeterPotentialEnergy(potentialMaster);
        meterTargInTarg.setBox(boxTarget);
        double latticeEnergyTarget = meterTargInTarg.getDataAsScalar();
        System.out.println("Target system lattice energy: " +latticeEnergyTarget);
        
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(boxTarget);
        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, meterPE,random, 
                space);
        mcMoveAtom.setPotential(potential);
        mcMoveAtom.setDoExcludeNonNeighbors(true);
        
        mcMoveAtom.setStepSize(0.01);
        integratorTarget.getMoveManager().addMCMove(mcMoveAtom);
        integratorTarget.setMeterPotentialEnergy(meterTargInTarg);
        
        
//JOINT
        //measuring potential of target in reference system
        meterTargInRef = new MeterDifferentImageAddDoubleSize(this,
                space,  temperature, cDefRef, nmRef, cDefTarget, potentialMaster, 
                size, nmTarg, tIn);
        MeterOverlapSameGaussian meterOverlapInRef = new 
                MeterOverlapSameGaussian("MeterOverlapInB", Null.DIMENSION, 
                meterRefInRef, meterTargInRef, temperature);
        meterOverlapInRef.setDsABase(latticeEnergyRef);
        meterOverlapInRef.setDsBBase(latticeEnergyTarget);
        
        //measuring reference potential in target system
        meterRefInTarg = new MeterDifferentImageSubtractDoubleSize(this, space, 
                cDefTarget, nmTarg, cDefRef, potentialMaster, new int[] {1,1,1},
                nmRef, rIn);
        MeterOverlap meterOverlapInTarget = new MeterOverlap("MeterOverlapInA", 
                Null.DIMENSION, meterTargInTarg, meterRefInTarg, temperature);
        meterOverlapInTarget.setDsABase(latticeEnergyTarget);
        meterOverlapInTarget.setDsBBase(latticeEnergyRef);
        
        P1ConstraintNbr nbrConstraint = new P1ConstraintNbr(space, 
                primitiveLength/Math.sqrt(2.0), this, constraint);
        potentialMaster.addPotential(nbrConstraint, new AtomType[]{
                species.getLeafType()});
        nbrConstraint.initBox(boxRef);
        nbrConstraint.initBox(boxTarget);
        nbrConstraint.initBox(meterTargInRef.getBox());
        nbrConstraint.initBox(meterRefInTarg.getBox());
        potentialMaster.getNeighborManager(boxRef).reset();
        potentialMaster.getNeighborManager(boxTarget).reset();
        potentialMaster.getNeighborManager(meterTargInRef.getBox()).reset();
        potentialMaster.getNeighborManager(meterRefInTarg.getBox()).reset();
        
        //Just to be sure!
        potential.setTruncationRadius(3000.0);
        
        meters[1] = meterOverlapInTarget;
        meters[0] = meterOverlapInRef;
        
        //Set up the rest of the joint stuff
        
        integratorSim = new IntegratorOverlap(new IntegratorMC[]{integratorRef,
                integratorTarget});
        
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        
        setBennettParameter(1.0, 30);
        
        activityIntegrate = new ActivityIntegrate(integratorSim, 0, true);
        getController().addAction(activityIntegrate);
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        System.out.println("OverallStart: " + System.currentTimeMillis());
        SimParam params = new SimParam();
        String inputFilename = null;
        if(args.length > 0) {
            inputFilename = args[0];
        }
        if(inputFilename != null){
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }

        double density = params.density;
        int D = params.D;
        String filename = params.filename;
        if(filename.length() == 0){
            filename = "nmi_3DSS_FCC_";
        }
        String inputFile = params.inputFile;
        double temperature = params.temperature;
        int runNumSteps = params.numSteps;
        int subBlockSize = params.subBlockSize;
        int eqNumSteps = params.eqNumSteps;
        int benNumSteps = params.bennettNumSteps;
        int exp = params.exponent;
        double constr = params.constraint;
        boolean first = params.first;
        int[] refCells = params.refShape;
        int[] targCells = params.targShape;
        int nRefA = 1;
        int nTargA = 1;
        for(int i = 0; i < D; i++){
            nRefA *= refCells[i];
            nTargA *= targCells[i];
        }
        nRefA *= 4;     //definitely fcc
        nTargA *= 4;    //definitely fcc

        filename = filename + "_" + nRefA + "_" + nTargA + "_" + temperature;


        //        int numberOfBlocks = params.numberOfBlocks;
//        int runNumSteps = nTargA * runBlockSize * numberOfBlocks * 2;
        int runBlockSize = runNumSteps / nTargA /100;
        System.out.println("RBS "+ runBlockSize);


        // instantiate simulation
        SimDifferentImageSsFccDoubleSize sim = new SimDifferentImageSsFccDoubleSize(
                Space.getInstance(D), refCells, targCells, density,
                temperature, exp, inputFile, constr);
        System.out.println("Dimension " + sim.space.D());
        System.out.println("Temperature " + temperature);
        System.out.println("Constraint " + constr);
        System.out.println("Ref system is " +nRefA + " atoms at density " + density);
        System.out.println("Targ system is " +nTargA + " atoms at density " + density);
        System.out.println("Add scaling: " + sim.meterTargInRef.getScaling());
        System.out.println("Sub scaling: " + sim.meterRefInTarg.getScaling());
        System.out.println(runNumSteps + " steps, " + runBlockSize + " blocksize");
        System.out.println("Target input data from " + inputFile);
        System.out.println("output data to " + filename);
        System.out.println("instantiated");

        if(false) {
//            SimulationGraphic graphic = new SimulationGraphic(sim, sim.space,
//                    sim.getController());
//            graphic.makeAndDisplayFrame();
//            return;


            SimulationGraphic simGraphic = new SimulationGraphic(sim,
                    SimulationGraphic.TABBED_PANE);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if(a.getLeafIndex() == 56) return Color.WHITE;
                    if(a.getLeafIndex() == 57) return Color.YELLOW;
                    if(a.getLeafIndex() == 58) return Color.CYAN;
                    if(a.getLeafIndex() == 59) return Color.BLUE;
                    if(a.getLeafIndex() == 42) return Color.GREEN;
                    if(a.getLeafIndex() == 43) return Color.GRAY;
                    if(a.getLeafIndex() == 50) return Color.GREEN;
                    if(a.getLeafIndex() == 49) return Color.MAGENTA;
                    if(a.getLeafIndex() == 34) return Color.CYAN;
                    if(a.getLeafIndex() == 31) return Color.GRAY;
                    if(a.getLeafIndex() == 29) return Color.MAGENTA;
                    if(a.getLeafIndex() == 21) return Color.YELLOW;
                    if(a.getLeafIndex() == 15) return Color.BLUE;


                    if(true) return Color.red;
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(16*a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.boxRef).setColorScheme(colorScheme);
            simGraphic.getDisplayBox(sim.boxTarget).setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();
            return;
        }

        //Divide out all the steps, so that the subpieces have the proper # of steps
        runNumSteps /= subBlockSize;
        eqNumSteps /= subBlockSize;
        benNumSteps /= subBlockSize;

        System.out.println("run " + runNumSteps);
        System.out.println("ben " + benNumSteps);
        System.out.println(" ea " + eqNumSteps);
//        sim.integratorSim.setAdjustStepFreq(false);
//        sim.integratorSim.setStepFreq0(0.5);

        //start simulation & equilibrate
        sim.integratorSim.getMoveManager().setEquilibrating(true);
        sim.integratorSim.setNumSubSteps(subBlockSize);


        System.out.println("EquilStart: " + System.currentTimeMillis());
        if(first){
            System.out.println("Init Bennett");
            sim.initBennettParameter(filename, benNumSteps, runBlockSize);
            if (Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 ||
                    Double.isInfinite(sim.bennettParam)){
                throw new RuntimeException("Simulation failed to find a valid " +
                        "Bennett parameter");
            }

            System.out.println("equilibrate");
            sim.equilibrate("bennett" , eqNumSteps, runBlockSize);
            if (Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 ||
                    Double.isInfinite(sim.bennettParam)){
                throw new RuntimeException("Simulation failed to find a valid " +
                        "Bennett parameter");
            }
            System.out.println("equilibration finished.");
        } else {
            System.out.println("Init Bennett");
            sim.initBennettParameter("bennett", benNumSteps, runBlockSize);
            System.out.println("equilibrate");
            sim.equilibrate(null, eqNumSteps, runBlockSize);
            System.out.println("equilibration finished.");
        }

        // start simulation
        long t1 = System.currentTimeMillis();
        sim.setAccumulatorBlockSize(runBlockSize);
        sim.integratorSim.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.setMaxSteps(runNumSteps);
        sim.getController().actionPerformed();
        System.out.println("final reference optimal step frequency " +
                sim.integratorSim.getIdealRefStepFraction() + " (actual: " +
                sim.integratorSim.getRefStepFraction() + ")");


        //CALCULATION OF HARMONIC ENERGY
        long t2 = System.currentTimeMillis();
        System.out.println("Calc: " + (t2-t1)/1000.0);
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+", error: "+error);
        DataGroup allYourBase =
            (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("reference ratio average (unscaled): " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1] + " error: " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1]
                .getNBennetPoints() - sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average (unscaled): " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        System.out.println("calculated diff " + (temperature*
                (-Math.log(ratio)
                        - sim.meterTargInRef.getScaling()
                        - 0.5 * sim.space.D() * (nTargA - nRefA) * Math.log(2 * Math.PI * temperature)
                - 0.5 * sim.space.D() * Math.log(nTargA)
                + 0.5 * sim.space.D() * Math.log(nRefA))));

//        System.out.println("new "+ (-Math.log(ratio)
//                - sim.meterTargInRef.getScaling() ));
//        System.out.println("old " +(-Math.log(ratio*)));


        System.out.println("Fini. ");
        System.out.println("End Time: " + System.currentTimeMillis());
    }

    public void setBennettParameter(double benParamCenter, double span) {
        bennettParam = benParamCenter;
        accumulators[0].setBennetParam(benParamCenter, span);
        accumulators[1].setBennetParam(benParamCenter, span);
    }

    public void setBennettParameter(double newBennettParameter) {
        System.out.println("setting ref pref (explicitly) to " +
                newBennettParameter);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
        setBennettParameter(newBennettParameter, 1);

    }

    public void initBennettParameter(String fileName, int initSteps, int initBlockSize) {
        // benParam = -1 indicates we are searching for an appropriate value
        bennettParam = -1.0;
        integratorSim.getMoveManager().setEquilibrating(true);

        if (fileName != null) {
            try {
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String benParamString = bufReader.readLine();
                bennettParam = Double.parseDouble(benParamString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to " + bennettParam);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
                setBennettParameter(bennettParam, 1);
            } catch (IOException e) {
                System.out.println("Bennett parameter not from file");
                // file not there, which is ok.
            }
        }

        if (bennettParam == -1) {

            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps);

            getController().actionPerformed();
            getController().reset();

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize, 41, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize, 41, false), 1);
            setBennettParameter(1, 10);
            activityIntegrate.setMaxSteps(initSteps);

            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            bennettParam = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);

            if (Double.isNaN(bennettParam) || bennettParam == 0 ||
                    Double.isInfinite(bennettParam)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to " + bennettParam);

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11, false), 1);
            setBennettParameter(bennettParam, 2);

            // set benParam back to -1 so that later on we know that we've been looking for
            // the appropriate value
            bennettParam = -1;
            getController().reset();
        }
        integratorSim.getMoveManager().setEquilibrating(false);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage
                                       newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            pumpListener.setInterval(getBox(iBox).getLeafList().getAtomCount());
            integrators[iBox].getEventManager().addListener(pumpListener);
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorSim != null && accumulators[0] != null &&
                accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0], accumulators[1]);
            integratorSim.setReferenceFracSource(dsvo);
        }

    }

    public void setAccumulatorBlockSize(int newBlockSize) {
        for (int i = 0; i < 2; i++) {
            accumulators[i].setBlockSize(newBlockSize);
        }
        try {
            // reset the integrator so that it will re-adjust step frequency
            // and ensure it will take enough data for both ref and target
            integratorSim.reset();
        } catch (ConfigurationOverlapException e) { /* meaningless */ }
    }

    public void equilibrate(String fileName, int initSteps, int initBlockSize) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        integratorSim.getMoveManager().setEquilibrating(true);

        for (int i = 0; i < 2; i++) {
            integrators[i].getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i = 0; i < 2; i++) {
            integrators[i].getMoveManager().setEquilibrating(false);
        }

        if (bennettParam == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            bennettParam = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to " + bennettParam + " (" + newMinDiffLoc + ")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize, 1, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize, 1, false), 1);

            setBennettParameter(bennettParam, 1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(bennettParam) + "\n");
                    bufWriter.close();
                    fileWriter.close();
                } catch (IOException e) {
                    throw new RuntimeException("couldn't write to Bennet parameter file");
                }
            }
        } else {
            dsvo.reset();
        }
        integratorSim.getMoveManager().setEquilibrating(false);
    }
    
    public static class SimParam extends ParameterBase {
        public boolean first = false;
        public int[] refShape = {2, 2, 2};
        public int[] targShape = {2, 2, 4};
        public double density = 1.1964;
        public int D = 3;
        public double harmonicFudge = 1.0;
        public double temperature = 0.01;
        public int exponent = 12;
        public double constraint = 2.0;
        
        public String inputFile = "inputSSDB_DS";
        public String filename = "output";
        
        public int numSteps = 100000000;     //overall # of steps of subintegrators
        public int subBlockSize = 10000;    //# of steps in subintegrator per integrator step
        public int eqNumSteps = 100000;  
        public int bennettNumSteps = 50000;
    }
}
