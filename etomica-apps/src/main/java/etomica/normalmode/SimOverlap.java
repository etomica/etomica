/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P1HardPeriodic;
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
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

import java.io.*;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * @author Andrew Schultz
 */
public class SimOverlap extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public MCMoveHarmonic move;
    public IntegratorHard integratorTarget;
    public PotentialMasterList potentialMasterTarget;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public SimOverlap(Space _space, int numAtoms, double density, double temperature, String filename, double harmonicFudge) {
        super(_space);

        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);

        species.setIsDynamic(true);
        addSpecies(species);

        // TARGET
        potentialMasterTarget = new PotentialMasterList(this, space);
        boxTarget = new Box(space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);

        integratorTarget = new IntegratorHard(this, potentialMasterTarget, space);
        // needs to be irrational so that configurations don't repeat in 1D with 2 atoms
        integratorTarget.setTimeStep(Math.PI);
        integratorTarget.setTemperature(1.0);
        integratorTarget.setBox(boxTarget);

        integratorTarget.setIsothermal(false);
        integrators[1] = integratorTarget;

        Potential2 p2 = new P2HardSphere(space, 1.0, false);
        if (space.D() == 1) {
            // don't need this for the target system, but we do need it for the reference
            p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        }
        potentialMasterTarget.addPotential(p2, new AtomType[]{species.getLeafType(), species.getLeafType()});

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms/density);
            double L = boundaryTarget.getBoxSize().getX(0);
            if (L-(numAtoms-1) > 0.5*L) {
                // at low density, the empty space between two atoms can exceed half the box length
                P1HardPeriodic p1h = new P1HardPeriodic(space);
                p1h.setSigma(1.0);
                integratorTarget.setNullPotential(p1h, species.getLeafType());
            }
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundaryTarget = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
            primitive = new PrimitiveCubic(space, n*L);
        }
        boxTarget.setBoundary(boundaryTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, primitive, basis, space);
        if (space.D() == 1) {
            coordinateDefinitionTarget.initializeCoordinates(nCells);
        }
        else {
            coordinateDefinitionTarget.initializeCoordinates(new int[]{1,1,1});
        }

        double neighborRange;
        if (space.D() == 1) {
            neighborRange = 1.01 / density;
        }
        else {
            //FCC
            double L = Math.pow(4.01/density, 1.0/3.0);
            neighborRange = L / Math.sqrt(2.0);
        }
        potentialMasterTarget.setRange(neighborRange);
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMasterTarget.getNeighborManager(boxTarget).setDoApplyPBC(false);
        potentialMasterTarget.getNeighborManager(boxTarget).reset();

        integratorTarget.setBox(boxTarget);

        // HARMONIC
        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, numAtoms/density);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, n * L);
        }
        boxHarmonic = new Box(boundaryHarmonic, space);
        addBox(boxHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(null, random, 1.0);
        integratorHarmonic.setBox(boxHarmonic);

        integrators[0] = integratorHarmonic;

        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(boxHarmonic, primitive, basis, space);
        if (space.D() == 1) {
            coordinateDefinitionHarmonic.initializeCoordinates(nCells);
            normalModes = new NormalModes1DHR(boundaryTarget, numAtoms);
        }
        else {
            coordinateDefinitionHarmonic.initializeCoordinates(new int[]{1,1,1});
            normalModes = new NormalModesFromFile(filename, space.D());
        }

        normalModes.setHarmonicFudge(harmonicFudge);
        normalModes.setTemperature(temperature);

        move = new MCMoveHarmonic(getRandom());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setTemperature(temperature);

        move.setBox(boxHarmonic);

        integratorHarmonic.getMoveManager().addMCMove(move);

        integratorHarmonic.setBox(boxHarmonic);

        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMasterTarget.getNeighborManager(boxHarmonic).setDoApplyPBC(false);
        potentialMasterTarget.getNeighborManager(boxHarmonic).reset();

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorHarmonic, integratorTarget});
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(new MeterPotentialEnergyFromIntegrator(integratorTarget), meterHarmonicEnergy);
        meterTarget.setTemperature(temperature);
        // lattice energy is 0
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMasterTarget);
        meterHarmonic.setTemperature(temperature);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);

        setRefPref(1.0, 30);

        activityIntegrate = new ActivityIntegrate(integratorOverlap);

        getController().addAction(activityIntegrate);
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlap.SimOverlapParam
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
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
            filename = "normal_modes_HS_"+D+"D_"+numMolecules;
        }
        String refFileName = args.length > 0 ? filename+"_ref" : null;

        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" hard sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        SimOverlap sim = new SimOverlap(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge);

        int totalCells = 1;
        for (int i=0; i<D; i++) {
            totalCells *= sim.nCells[i];
        }
        int basisSize = sim.basis.getScaledCoordinates().length;

        double AHarmonic = CalcHarmonicA.doit(sim.normalModes, D, temperature, basisSize*totalCells);

        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;

        sim.initRefPref(refFileName, numSteps/20);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }

        sim.equilibrate(refFileName, numSteps/10);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }

        System.out.println("equilibration finished");

        IDataSource[] workMeters = new IDataSource[2];

        //Harmonic
        MeterWorkHarmonicPhaseSpace meterWorkHarmonic = new MeterWorkHarmonicPhaseSpace(sim.move, sim.potentialMasterTarget);
        meterWorkHarmonic.setTemperature(temperature);
        meterWorkHarmonic.setLatticeEnergy(0); //Hard Sphere Lattice Energy = 0 (no overlap)
        workMeters[0] = meterWorkHarmonic;

        DataFork dataForkHarmonic = new DataFork();
        DataPumpListener pumpHarmonic = new DataPumpListener(workMeters[0], dataForkHarmonic);

        AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        /*
        AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramCollapsing());
        dataForkHarmonic.addDataSink(histogramHarmonic);
        */
        sim.integrators[0].getEventManager().addListener(pumpHarmonic);
        //

        //Target
        MeterWorkTargetPhaseSpace meterWorkTarget = new MeterWorkTargetPhaseSpace(sim.integratorTarget,sim.meterHarmonicEnergy);
        meterWorkTarget.setLatticeEnergy(0); //Hard Sphere Lattice Energy = 0 (no overlap)
        workMeters[1] = meterWorkTarget;

        DataFork dataForkTarget = new DataFork();
        DataPumpListener pumpTarget = new DataPumpListener(workMeters[1], dataForkTarget, 100);

        AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        dataForkTarget.addDataSink(dataAverageTarget);
        /*
        AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramCollapsing());
        dataForkTarget.addDataSink(histogramTarget);
        */
        sim.integrators[1].getEventManager().addListener(pumpTarget);

        //


        sim.integratorOverlap.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        System.out.println("ideal reference optimal step frequency "+sim.integratorOverlap.getIdealRefStepFraction()+" (actual: "+sim.integratorOverlap.getRefStepFraction()+")");

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: "+(-Math.log(ratio))+", error: "+(error/ratio));
        System.out.println("target free energy: "+(AHarmonic-Math.log(ratio)));
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);
/*
        double betaUab = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
        double betaUba = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);

        double err_betaUab = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
        double err_betaUba = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);

        System.out.println("betaUab: "+betaUab + " ,err: "+err_betaUab);
        System.out.println("betaUba: "+betaUba + " ,err: "+err_betaUba);
  */
        if(D==1) {
            double AHR = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules) + SpecialFunctions.lnFactorial(numMolecules-1) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }

        numMolecules++;
        double AHRNp1 = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules) + SpecialFunctions.lnFactorial(numMolecules-1) ;
        System.out.println("Hard-rod free energy ("+numMolecules+"): "+AHRNp1);

    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter, span);
        accumulators[1].setBennetParam(refPrefCenter, span);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            integrators[iBox].getEventManager().addListener(new IntegratorListenerAction(accumulatorPumps[iBox]));
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0], accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
        }
    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to " + newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
        setRefPref(newRefPref, 1);
    }

    public void initRefPref(String fileName, long initSteps) {
        // refPref = -1 indicates we are searching for an appropriate value
        refPref = -1.0;
        if (fileName != null) {
            try {
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to " + refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
                setRefPref(refPref, 1);
            } catch (IOException e) {
                // file not there, which is ok.
            }
        }

        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps / 2);
            getController().actionPerformed();
            getController().reset();

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41, false), 1);
            setRefPref(1, 40);
            activityIntegrate.setMaxSteps(initSteps);
            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to " + refPref);

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11, false), 1);
            setRefPref(refPref, 5);

            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
            getController().reset();
        }

    }

    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC)
                ((IntegratorMC) integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC)
                ((IntegratorMC) integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to " + refPref + " (" + newMinDiffLoc + ")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
            setRefPref(refPref, 1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref) + "\n");
                    bufWriter.close();
                    fileWriter.close();
                } catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        } else {
            dsvo.reset();
        }
    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 2;
        public double density = 0.7;
        public int D = 1;
        public long numSteps = 100000;
        public double harmonicFudge = 1;
        public String filename = "normal_modes3D_32_130_cubic";
        public double temperature = 1.0;
    }
}
