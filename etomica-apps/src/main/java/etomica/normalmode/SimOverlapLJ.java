/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.*;
import etomica.nbr.list.PotentialMasterList;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.*;
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
 * Simulation to run sampling with the LJ potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * @author Andrew Schultz
 */
public class SimOverlapLJ extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;

    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive, primitiveUnitCell;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public P1Constraint p1Constraint;
    public SimOverlapLJ(Space _space, int numAtoms, double density, double temperature, double harmonicFudge) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        PotentialMaster potentialMasterTarget = new PotentialMasterMonatomic(this);
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        // TARGET
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0 / density);
            boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms / density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            primitive = new PrimitiveCubic(space, n * L);
            primitiveUnitCell = new PrimitiveCubic(space, L);

            nCells = new int[]{n, n, n};
            boundaryTarget = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
        }
        boxTarget = this.makeBox(boundaryTarget);
        boxTarget.setNMolecules(species, numAtoms);

        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget, getRandom(), temperature, boxTarget);
        MCMoveAtomCoupled atomMove = new MCMoveAtomCoupled(potentialMasterTarget, new MeterPotentialEnergy(potentialMasterTarget), getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        integratorTarget.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker) atomMove.getTracker()).setNoisyAdjustment(true);

        integrators[1] = integratorTarget;


        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potential = new P2LennardJones(space, 1.0, 1.0);
        double truncationRadius = boundaryTarget.getBoxSize().getX(0) * 0.45;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();
        potentialMasterTarget.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});
        atomMove.setPotential(pTruncated);


        /*
         * 1-body Potential to Constraint the atom from moving too far away from
         * its lattice-site
         */
        p1Constraint = new P1Constraint(space, primitiveUnitCell.getSize()[0], boxTarget, coordinateDefinitionTarget);
        potentialMasterTarget.addPotential(p1Constraint, new AtomType[]{sphereType});

        if (potentialMasterTarget instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList) potentialMasterTarget).setRange(neighborRange);
            ((PotentialMasterList) potentialMasterTarget).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMasterTarget).getNeighborManager(boxTarget).reset();
            int potentialCells = ((PotentialMasterList) potentialMasterTarget).getNbrCellManager(boxTarget).getLattice().getSize()[0];
            if (potentialCells < cellRange * 2 + 1) {
                throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
            }
            if (potentialCells > cellRange * 2 + 1) {
                System.out.println("could probably use a larger truncation radius (" + potentialCells + " > " + (cellRange * 2 + 1) + ")");
            }
        }

        potentialMasterTarget.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget, boxTarget);
        double latticeEnergy = meterPE.getDataAsScalar();

        // HARMONIC
        boundaryHarmonic = new BoundaryRectangularPeriodic(space);
        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, numAtoms / density);
        } else {
            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, n * L);
        }
        boxHarmonic = this.makeBox(boundaryHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(potentialMasterTarget, random, 1.0, boxHarmonic);

        MCMoveHarmonic move = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;

        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(boxHarmonic, primitive, basis, space);
        coordinateDefinitionHarmonic.initializeCoordinates(new int[]{1, 1, 1});

        String inFile = "LJDB_nA" + numAtoms + "_d0962";
        normalModes = new NormalModesFromFile(inFile, space.D());
        normalModes.setHarmonicFudge(harmonicFudge);
        //normalModes.setTemperature(temperature);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);
        move.setTemperature(temperature);

        move.setBox(boxHarmonic);

        if (potentialMasterTarget instanceof PotentialMasterList) {
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMasterTarget).getNeighborManager(boxHarmonic).reset();
        }

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorHarmonic, integratorTarget});
        MeterHarmonicEnergy meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(new MeterPotentialEnergyFromIntegrator(integratorTarget), meterHarmonicEnergy);
        meterTarget.setTemperature(temperature);
        meterTarget.setLatticeEnergy(latticeEnergy);
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMasterTarget);
        meterHarmonic.setTemperature(temperature);
        meterHarmonic.setLatticeEnergy(latticeEnergy);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);

        setRefPref(1.0, 30);

        this.getController2().addActivity(new ActivityIntegrate2(integratorOverlap));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapLJ.SimOverlapParam
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
            filename = "normal_modes_LJ_3D_" + numMolecules;
        }
        String refFileName = args.length > 0 ? filename + "_ref" : null;
        /*
        if(density == 0.962){
        	inputFile = "LJDB_nA"+numMolecules+"_d0962";
        } else if (density == 1.002){
        	inputFile = "LJDB_nA"+numMolecules+"_d1002";
        } else if (density == 1.211) {
        	inputFile = "LJDB_nA"+numMolecules+"_d1211";
        } else {
        	System.err.println("NO INPUT FILE FOUND!!!!");
        }
        */
        System.out.println("Running " + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal")) + " LJ overlap simulation");
        System.out.println(numMolecules + " atoms at density " + density + " ;temperature " + temperature);
        System.out.println("harmonic fudge: " + harmonicFudge);
        System.out.println((numSteps / 1000) + " total steps of 1000");
        //System.out.println("input File "+inputFile);
        System.out.println("output data to " + filename);

        //instantiate simulation
        SimOverlapLJ sim = new SimOverlapLJ(Space.getInstance(D), numMolecules, density, temperature, harmonicFudge);

        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;

        sim.initRefPref(refFileName, numSteps / 20);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
        System.out.flush();

        sim.equilibrate(refFileName, numSteps / 10);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }

        System.out.println("equilibration finished");
        System.out.flush();

        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorOverlap), numSteps);

        int totalCells = 1;
        for (int i = 0; i < D; i++) {
            totalCells *= sim.nCells[i];
        }
        int basisSize = sim.basis.getScaledCoordinates().length;

        double AHarmonic = CalcHarmonicA.doit(sim.normalModes, D, temperature, basisSize * totalCells);
        System.out.println("Harmonic-reference free energy, A: " + AHarmonic + " " + AHarmonic / numMolecules);
        System.out.println(" ");


        System.out.println("final reference optimal step frequency " + sim.integratorOverlap.getIdealRefStepFraction() + " (actual: " + sim.integratorOverlap.getRefStepFraction() + ")");

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: " + ratio + ", error: " + error);
        System.out.println("free energy difference: " + (-temperature * Math.log(ratio)) + ", error: " + temperature * (error / ratio));
        System.out.println("target free energy: " + (AHarmonic - temperature * Math.log(ratio)));
        System.out.println("target free energy per particle: " + (AHarmonic - temperature * Math.log(ratio)) / numMolecules
                + " ;error: " + temperature * (error / ratio) / numMolecules);
        DataGroup allYourBase = (DataGroup) sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: " + ((DataDoubleArray) allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(sim.accumulators[0].STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[1]);

        allYourBase = (DataGroup) sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        System.out.println("target ratio average: " + ((DataDoubleArray) allYourBase.getData(sim.accumulators[1].AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(sim.accumulators[1].STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(sim.accumulators[1].ERROR.index)).getData()[1]);

        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken: " + (endTime - startTime));

    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
        // needed for Bennett sampling
//        if (accumulators[0].getNBennetPoints() == 1) {
//            ((MeterBoltzmannHarmonic)meters[0]).refPref = refPrefCenter;
//            ((MeterBoltzmannTarget)meters[1]).refPref = refPrefCenter;
//        }
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        newAccumulator.setBlockSize(100);
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            if (iBox == 1) {
                pumpListener.setInterval(boxTarget.getMoleculeList().size());
            }
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
        }
    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
        setRefPref(newRefPref,1);
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
                System.out.println("setting ref pref (from file) to "+refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }

        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            getController2().runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps/2);

            System.out.println("target equilibration finished");

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1,60);
getController2().runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps);


            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
            setRefPref(refPref,5);

            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;

        }

    }

    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        this.getController2().runActivityBlocking(new ActivityIntegrate2(this.integratorOverlap), initSteps);

        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
            setRefPref(refPref,1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        }
        else {
            dsvo.reset();
        }
    }
    //private static String inputFile;

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 32;
        public double density = 0.962;
        public int D = 3;
        public long numSteps = 10000000;
        public double harmonicFudge = 1;
        public String filename = "LJCB_nA32_d0962_T10";
        public double temperature = 1.0;
    }
}
