/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;


import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicSimple;
import etomica.nbr.site.PotentialMasterSite;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;


/**
 * Compute  the dielectric constant of  2D interacting Heisenberg model in conventional way and mapped averaging.
 * Conventional way: get the -bt*bt*<M^2> with M is the total dipole moment;
 * Mapped averaging: get the A_EE secondDerivative of free energy w.r.t electric field E when E is zero
 * Prototype for simulation of a more general magnetic system.
 *
 * @author Weisong Lin
 */
public class Heisenberg extends Simulation {

    private static final String APP_NAME = "Heisenberg";
    public PotentialMasterSite potentialMaster; // difference betweet Pmaster pmastersite
    public Box box;
    public SpeciesSpheresMono spins;
    public P2Spin potential;
    public MCMoveRotate mcMove;
    private IntegratorMC integrator;

    /**
     * 2D heisenberg model in square lattice.
     *
     * @param space           use to define vector/tensor
     * @param nCells          total number of atoms = nCells*nCells
     * @param interactionS    the J in heisenberg energy function: U = J*Cos(theta1-theta2)
     * @param dipoleMagnitude is the strength of heisenberg dipole.
     */
    public Heisenberg(Space space, int nCells, double temperature, double interactionS, double dipoleMagnitude) {
        super(Space2D.getInstance());
//        setRandom(new RandomNumberGenerator(1)); //debug only
//        System.out.println("============================the RandomSeed is one ===========================");

        spins = new SpeciesSpheresRotating(space, new ElementSimple("A"));

        addSpecies(spins);

        potentialMaster = new PotentialMasterSite(this, nCells, space);
        box = new Box(space);
        addBox(box);
        int numAtoms = space.powerD(nCells);

        box.setNMolecules(spins, numAtoms);
        box.getBoundary().setBoxSize(new Vector2D(nCells, nCells));
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(space, 1), space);
        config.initializeCoordinates(box);

        potential = new P2Spin(space, interactionS);
        integrator = new IntegratorMC(this, potentialMaster, box);
        mcMove = new MCMoveRotate(potentialMaster, random, space);
        integrator.getMoveManager().addMCMove(mcMove);
        MCMoveSpinCluster spinMove = new MCMoveSpinCluster(space, random, potentialMaster, integrator, interactionS);
        integrator.getMoveManager().addMCMove(spinMove);
        integrator.setTemperature(temperature);
        this.getController2().addActivity(new ActivityIntegrate2(integrator));
        AtomType type = spins.getLeafType();
//        potentialMaster.addPotential(field, new IAtomType[] {type});
        potentialMaster.addPotential(potential, new AtomType[]{type, type});

    }

    public static void main(String[] args) {
        Param params = new Param();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doConventionalE = false;
            params.doMappingE = false;
            params.doEnergyMF = false;
            params.doCorrelation = false;
            params.formula = 3;
            params.doCV = true;
            params.aEE = true;
            params.mSquare = true;
            params.doVSumMI = true;
            params.doIdeal = true;
            params.doAEEMF = true;
            params.temperature = 5;
            params.nCells = 5;
            params.steps = 10000000;
            params.doGraphic = true;

        }
        final long startTime = System.currentTimeMillis();
        DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Calendar cal = Calendar.getInstance();
        System.out.println("startTime : " + date.format(cal.getTime()));

        double temperature = params.temperature;
        int nCells = params.nCells;
        int numberMolecules = nCells * nCells;
        boolean aEE = params.aEE;
        boolean mSquare = params.mSquare;
        boolean doIdeal = params.doIdeal;
        boolean doPair = params.doPair;
        boolean doVSum = params.doVSum;
        boolean doVSumMI = params.doVSumMI;
        boolean doConventionalE = params.doConventionalE;
        boolean doMappingE = params.doMappingE;
        boolean doCV = params.doCV;
        boolean doCorrelation = params.doCorrelation;
        boolean doDipole = params.doDipole;
        boolean doEnergyMF = params.doEnergyMF;
        boolean doAEEMF = params.doAEEMF;
        int formula = params.formula;
        boolean doGraphic = params.doGraphic;


        int steps = params.steps;
        int nMax = params.nMax;
        double interactionS = params.interactionS;
        double dipoleMagnitude = params.dipoleMagnitude;

        System.out.println("numberMolecules= " + numberMolecules + " steps= " + steps + " nMax= " + nMax
                + " \ntemperature= " + temperature + " interactionS= " + interactionS + " dipoleStrength= " + dipoleMagnitude);

        Space sp = Space2D.getInstance();
        Heisenberg sim = new Heisenberg(sp, nCells, temperature, interactionS, dipoleMagnitude);

        if (doGraphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "XY", numberMolecules);
            ((DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.getSpecies(0).getAtomType(0), 0.2);
            ((DisplayBoxCanvas2D) simGraphic.getDisplayBox(sim.box).canvas).setOrientationColor(Color.BLUE);
            ((DisplayBoxCanvas2D) simGraphic.getDisplayBox(sim.box).canvas).setOrientationLength(5);

            DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integrator, "temperature");
            temperatureSlider.setLabel("Temperature");
            temperatureSlider.setShowBorder(true);
            temperatureSlider.setMinimum(0);
            temperatureSlider.setMaximum(10);
            temperatureSlider.setPrecision(1);
            temperatureSlider.setNMajor(5);
            temperatureSlider.setNMinor(2);
            temperatureSlider.setShowValues(true);
            temperatureSlider.setEditValues(true);
            simGraphic.add(temperatureSlider);

            simGraphic.makeAndDisplayFrame();
            return;
        }


        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), steps / 5);
sim.getController().reset();
        int blockNumber = 100;

        System.out.println("equilibration finished");
        long equilibrationTime = System.currentTimeMillis();
        System.out.println("equilibrationTime: " + (equilibrationTime - startTime) / (1000.0 * 60.0));

        int sampleAtInterval = numberMolecules * 2;
        int samplePerBlock = steps / sampleAtInterval / blockNumber;


        MeterMappedAveragingCV CVMeter = null;
        AccumulatorAverageCovariance CVAccumulator = null;
        if (doCV) {
            CVMeter = new MeterMappedAveragingCV(sim, temperature, interactionS, sim.potentialMaster);
            CVAccumulator = new AccumulatorAverageCovariance(samplePerBlock);
            DataPumpListener CVPump = new DataPumpListener(CVMeter, CVAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(CVPump);
        }


        MeterPotentialEnergy FECon = null;
        AccumulatorAverageFixed FEConAccumulator = null;
        if (doConventionalE) {
            FECon = new MeterPotentialEnergy(sim.potentialMaster);
            FECon.setBox(sim.box);
            FEConAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPumpListener FEConPump = new DataPumpListener(FECon, FEConAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(FEConPump);
        }

        MeterMappedAveragingFreeEnergy FEMap = null;
        AccumulatorAverageFixed FEMapAccumulator = null;

        if (doMappingE) {
            FEMap = new MeterMappedAveragingFreeEnergy(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster);
            FEMapAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPumpListener FEMapPump = new DataPumpListener(FEMap, FEMapAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(FEMapPump);
        }


        MeterSpinMSquare meterMSquare = null;
        AccumulatorAverage dipoleSumSquaredAccumulator = null;
        //mSquare
        if (mSquare) {
            meterMSquare = new MeterSpinMSquare(sim.space, sim.box, dipoleMagnitude);
            dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPump dipolePump = new DataPump(meterMSquare, dipoleSumSquaredAccumulator);
            IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
            dipoleListener.setInterval(sampleAtInterval);
            sim.integrator.getEventManager().addListener(dipoleListener);


        }


        //AEE
        MeterMappedAveragingVSum AEEMeter = null;
        AccumulatorAverageCovariance AEEAccumulator = null;
        if (aEE) {
            AEEMeter = new MeterMappedAveragingVSum(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster, doIdeal, doPair, doVSum, doVSumMI, doAEEMF, nMax);
            AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPumpListener AEEListener = new DataPumpListener(AEEMeter, AEEAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(AEEListener);
        }


        //Correlation
        MeterMappedAveragingCorrelation CorrelationMeter = null;
        AccumulatorAverageFixed CorrelationAccumulator = null;

        if (doCorrelation) {

            CorrelationMeter = new MeterMappedAveragingCorrelation(sim, temperature, interactionS, sim.potentialMaster, formula);
            CorrelationAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPumpListener CorrelationPump = new DataPumpListener(CorrelationMeter, CorrelationAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(CorrelationPump);
        }

        AccumulatorAverageFixed dipoleAccumulator = null;
        AccumulatorAverageFixed meanFieldAccumulator = null;

        if (doDipole) {
            MeterDipoleMoment dipoleMeter = new MeterDipoleMoment(sim.box);
            MeterMeanField meterMeanField = new MeterMeanField(sim.space, sim.box, interactionS, sim.potentialMaster, temperature);
            dipoleAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            meanFieldAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPumpListener pumpDipole = new DataPumpListener(dipoleMeter, dipoleAccumulator, sampleAtInterval);
            DataPumpListener pumpMeanField = new DataPumpListener(meterMeanField, meanFieldAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(pumpDipole);
            sim.integrator.getEventManager().addListener(pumpMeanField);
        }

        AccumulatorAverageCovariance energyMFAccumulator = null;

        if (doEnergyMF) {
            MeterEnergyMeanField meterEnergyMF = new MeterEnergyMeanField(sim.space, sim.box, interactionS, sim.potentialMaster, temperature);
            energyMFAccumulator = new AccumulatorAverageCovariance(samplePerBlock);
            DataPumpListener pumpEnergyMF = new DataPumpListener(meterEnergyMF, energyMFAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(pumpEnergyMF);
        }
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), steps);
        long endTime = System.currentTimeMillis();
        double totalTime = (endTime - startTime) / (1000.0 * 60.0);

        //******************************** simulation start ******************************** //
        //******************************** simulation start ******************************** //
        //******************************** simulation start ******************************** //
        //******************************** simulation start ******************************** //

        double energyMFAvg = 0;
        double energyMFErr = 0;
        if (doEnergyMF) {
            energyMFAvg = energyMFAccumulator.getData(energyMFAccumulator.AVERAGE).getValue(0);
            energyMFErr = energyMFAccumulator.getData(energyMFAccumulator.ERROR).getValue(0);
            double energyMFCor = energyMFAccumulator.getData(energyMFAccumulator.BLOCK_CORRELATION).getValue(0);
            System.out.println("energyMF:\t" + energyMFAvg / (nCells * nCells)
                    + " Err:\t" + energyMFErr / (nCells * nCells) + " Cor:\t " + energyMFCor
                    + " Difficulty:\t" + (energyMFErr * Math.sqrt(totalTime)));

            double cvMFAvg = energyMFAccumulator.getData(energyMFAccumulator.AVERAGE).getValue(1) - energyMFAvg * energyMFAvg / (temperature * temperature);
            double cv1Err = energyMFAccumulator.getData(energyMFAccumulator.ERROR).getValue(1);
            IData covariance = ((DataGroup) energyMFAccumulator.getData()).getData(energyMFAccumulator.BLOCK_COVARIANCE.index);

            double beta2 = 1.0 / (temperature * temperature);
            double cvMFErr = Math.sqrt(cv1Err * cv1Err
                    + 4 * beta2 * beta2 * energyMFAvg * energyMFAvg * energyMFErr * energyMFErr
                    - 4 * beta2 * cv1Err * energyMFAvg * energyMFErr * covariance.getValue(0 * 2 + 1) / Math.sqrt(covariance.getValue(0 * 2 + 0) * covariance.getValue(1 * 2 + 1))
            );
            System.out.println("CvMF:\t" + cvMFAvg / (nCells * nCells)
                    + " Err:\t" + cvMFErr / (nCells * nCells));

        }

        if (doDipole) {
            IData dipoleAvg = dipoleAccumulator.getData(dipoleAccumulator.AVERAGE);
            IData dipoleErr = dipoleAccumulator.getData(dipoleAccumulator.ERROR);
            IData dipoleCor = dipoleAccumulator.getData(dipoleAccumulator.BLOCK_CORRELATION);
            System.out.println("dipoleX:\t" + dipoleAvg.getValue(0)
                    + " Err:\t" + dipoleErr.getValue(0) + " Cor:\t " + dipoleCor.getValue(0)
                    + " Difficulty:\t" + (dipoleErr.getValue(0) * Math.sqrt(totalTime)));
            System.out.println("dipoleY:\t" + dipoleAvg.getValue(1)
                    + " Err:\t" + dipoleErr.getValue(1) + " Cor:\t " + dipoleCor.getValue(1)
                    + " Difficulty:\t" + (dipoleErr.getValue(1) * Math.sqrt(totalTime)));

            IData mapDipoleAvg = meanFieldAccumulator.getData(meanFieldAccumulator.AVERAGE);
            IData mapDipoleErr = meanFieldAccumulator.getData(meanFieldAccumulator.ERROR);
            IData mapDipoleCor = meanFieldAccumulator.getData(meanFieldAccumulator.BLOCK_CORRELATION);
            System.out.println("mapped dipoleX:\t" + mapDipoleAvg.getValue(0)
                    + " Err:\t" + mapDipoleErr.getValue(0) + " Cor:\t " + mapDipoleCor.getValue(0)
                    + " Difficulty:\t" + (mapDipoleErr.getValue(0) * Math.sqrt(totalTime)));
            System.out.println("mapped dipoleY:\t" + mapDipoleAvg.getValue(1)
                    + " Err:\t" + mapDipoleErr.getValue(1) + " Cor:\t " + mapDipoleCor.getValue(1)
                    + " Difficulty:\t" + (mapDipoleErr.getValue(1) * Math.sqrt(totalTime)));
        }

        //******************************** Correlation ******************************** //
        int distance = nCells / 2 + 1;
        int arraySize = -1 + (distance + 1) * distance / 2;
        double[] xValue = new double[arraySize];
        double[] xRow = new double[arraySize];
        double[] xCol = new double[arraySize];
        double[] xErr = new double[arraySize];
        double[] xCor = new double[arraySize];

        if (doCorrelation) {
            System.out.println("We're using formula********** " + formula + " ******************");
            for (int i = 0; i < arraySize; i++) {
                xValue[i] = ((DataGroup) CorrelationAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(i);
                xRow[i] = ((DataGroup) CorrelationAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(i + arraySize);
                xCol[i] = ((DataGroup) CorrelationAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(i + 2 * arraySize);
                xValue[i] = ((DataGroup) CorrelationAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(i);
                xErr[i] = ((DataGroup) CorrelationAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(i);
                xCor[i] = ((DataGroup) CorrelationAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(i);

                System.out.println("x[" + i + "] " + " row:\t" + ((int) xRow[i]) + "\tcol:\t" + ((int) xCol[i])
                        + "\t" + xValue[i] + " Err\t" + xErr[i] + " Cor\t" + xCor[i]);
            }
        }


        //******************************** Heat capacity ******************************** //
        double CV = 0;
        double CVErr = 0;
        double CVCor = 0;
        if (doCV) {
            double CV0 = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            double CV0Err = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            double CV0Cor = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(0);

            double CV1 = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
            double CV1Err = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
            double CV1Cor = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(1);

            System.out.println("CV0Err:\t" + CV0Err + " CV1Err:\t" + CV1Err);

            IData covariance = ((DataGroup) CVAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            CV = CV0 - CV1 * CV1;
            CVErr = Math.sqrt(4 * CV1 * CV1 * CV1Err * CV1Err
                    + CV0Err * CV0Err
                    - 4 * CV1 * CV0Err * CV1Err * covariance.getValue(0 * 2 + 1) / Math.sqrt(covariance.getValue(0 * 2 + 0) * covariance.getValue(1 * 2 + 1))
            );
            CV *= 1.0 / numberMolecules;
            CVErr *= 1.0 / numberMolecules;
            System.out.println("high-T CV/N:\t" + CV + " CVErr/N:\t" + CVErr);
        }


        //******************************** Conventional var[M]******************************** //
        if (mSquare) {
            double dipoleSumSquared = 0;
            double dipoleSumSquaredERR = 0;
            double dipoleSumCor = 0;
            dipoleSumSquared = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
            dipoleSumSquaredERR = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(AccumulatorAverage.ERROR.index)).x;
            dipoleSumCor = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index)).x;
            System.out.println("Conventional:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
                    + " Err:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules) + " Cor:\t " + dipoleSumCor
                    + " Difficulty:\t" + (dipoleSumSquaredERR * Math.sqrt(totalTime) / nCells / nCells));
        }

        //******************************** AEE vSum  ******************************** //

        if (aEE) {
            double AEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            double AEEER = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            double AEECor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(0);
//            System.out.println("PairMapping:\t" + (AEE / numberMolecules)
//                    + " Err:\t" + (AEEER / numberMolecules) + " Cor:\t " + AEECor
//                    + " Difficulty:\t" + (AEEER * Math.sqrt(totalTime) / nCells / nCells));
        }


        if (doIdeal) {
            double sumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
            double errSumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
            double sumIdealCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(1);
            System.out.println("IdealMapping:\t" + (sumIdeal / numberMolecules)
                    + " Err:\t" + (errSumIdeal / numberMolecules) + " Cor:\t " + sumIdealCor
                    + " Difficulty:\t" + (errSumIdeal * Math.sqrt(totalTime) / nCells / nCells));
        }


        //********************************AEE vSumMinusIdeal******************************** //

        if (doVSumMI) {
            double AEEVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(12);
            double AEEERVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(12);
            double AEECorVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(12);
            System.out.println("mappingvSumMI:\t" + (AEEVSumMinusIdeal / numberMolecules)
                    + " Err:\t" + (AEEERVSumMinusIdeal / numberMolecules) + " Cor:\t " + AEECorVSumMinusIdeal
                    + " Difficulty:\t" + (AEEERVSumMinusIdeal * Math.sqrt(totalTime) / nCells / nCells));
        }


        //******************************** Energy part ******************************** //

        double U_ConREF = 0;
        double U_ConREFERR = 0;
        if (doConventionalE) {
            U_ConREF = ((DataDouble) ((DataGroup) FEConAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
            U_ConREFERR = ((DataDouble) ((DataGroup) FEConAccumulator.getData()).getData(AccumulatorAverage.ERROR.index)).x;
            double U_ConREFCor = ((DataDouble) ((DataGroup) FEConAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index)).x;
            System.out.println("U_ConREF:\t" + (U_ConREF / numberMolecules)
                    + " Err:\t" + (U_ConREFERR / numberMolecules) + " Cor:\t " + U_ConREFCor
                    + " Difficulty:\t" + (U_ConREFERR * Math.sqrt(totalTime) / nCells / nCells));
        }

        double UMapREF = 0;
        double UMapREFERR = 0;
        if (doMappingE) {
            UMapREF = ((DataGroup) FEMapAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            UMapREFERR = ((DataGroup) FEMapAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            double UMapREFCor = ((DataGroup) FEMapAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(0);
            System.out.println("UMapREF:\t" + (UMapREF / numberMolecules)
                    + " Err:\t" + (UMapREFERR / numberMolecules) + " Cor:\t " + UMapREFCor
                    + " Difficulty:\t" + (UMapREFERR * Math.sqrt(totalTime) / nCells / nCells));
            System.out.println("U_ConMC and UMapREF is for checking purposes only");
        }

        if (doEnergyMF && (doConventionalE || doMappingE)) {
            System.out.println("Consistency check of MF energy with other formulas for energy (difference/uncertainty):");
            if (doConventionalE) {
                System.out.println("Conventional: " + (U_ConREF - energyMFAvg) / Math.sqrt(U_ConREFERR * U_ConREFERR + energyMFErr * energyMFErr));
            }
            if (doMappingE) {
                System.out.println("  High-T map: " + (UMapREF - energyMFAvg) / Math.sqrt(UMapREFERR * UMapREFERR + energyMFErr * energyMFErr));
            }
        }

        //******************************** Get FE from mapping meter ******************************** //
        if (aEE && doVSumMI) {
            double UMap = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(21);
            double UMapErr = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(21);
            double UMapCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(21);

            System.out.println("UMap:\t\t" + (UMap / numberMolecules)
                    + " Err:\t" + (UMapErr / numberMolecules) + " Cor:\t " + UMapCor
                    + " Difficulty:\t" + (UMapErr * Math.sqrt(totalTime) / nCells / nCells));

            double UMapSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(22);
            double UMapSquareErr = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(22);
            double UMapSquareCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(22);

            System.out.println("UMapSquare:\t" + (UMapSquare / numberMolecules)
                    + " Err:\t" + (UMapSquareErr / numberMolecules) + " Cor:\t " + UMapSquareCor
                    + " Difficulty:\t" + (UMapSquareErr * Math.sqrt(totalTime) / nCells / nCells));

            IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            double varUMap = UMapSquare - UMap * UMap;
            double varUMapErr = Math.sqrt(UMapSquareErr * UMapSquareErr
                    + 4 * UMap * UMap * UMapErr * UMapErr
                    - 4 * UMapSquareErr * UMap * UMapErr * covariance.getValue(21 * 25 + 22) / Math.sqrt(covariance.getValue(21 * 25 + 21) * covariance.getValue(22 * 25 + 22))
            );

            System.out.println("varUMap:\t" + (varUMap / numberMolecules)
                    + " Err:\t" + (varUMapErr / numberMolecules)
                    + " Difficulty:\t" + (varUMapErr * Math.sqrt(totalTime) / nCells / nCells));

        }

        if (doConventionalE) {
            double UCon = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(23);
            double UConErr = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(23);
            double UConCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(23);
            System.out.println("UCon:\t\t" + (UCon / numberMolecules)
                    + " Err:\t" + (UConErr / numberMolecules) + " Cor:\t " + UConCor
                    + " Difficulty:\t" + (UConErr * Math.sqrt(totalTime) / nCells / nCells));
            double UConSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(24);
            double UConSquareErr = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(24);
            double UConSquareCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(24);
            System.out.println("UConSquare:\t" + (UConSquare / numberMolecules)
                    + " Err:\t" + (UConSquareErr / numberMolecules) + " Cor:\t " + UConSquareCor
                    + " Difficulty:\t" + (UConSquareErr * Math.sqrt(totalTime) / nCells / nCells));

            IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            double varUCon = UConSquare - UCon * UCon;
            double varUConErr = Math.sqrt(UConSquareErr * UConSquareErr
                    + 4 * UCon * UCon * UConErr * UConErr
                    - 4 * UConSquareErr * UCon * UConErr * covariance.getValue(21 * 25 + 22) / Math.sqrt(covariance.getValue(21 * 25 + 21) * covariance.getValue(22 * 25 + 22))
            );

            System.out.println("varUCon:\t" + (varUCon / numberMolecules)
                    + " Err:\t" + (varUConErr / numberMolecules)
                    + " Difficulty:\t" + (varUConErr * Math.sqrt(totalTime) / nCells / nCells));

            System.out.println("CVCon:\t\t" + (varUCon / numberMolecules / temperature / temperature)
                    + " Err:\t" + (varUConErr / numberMolecules / temperature / temperature));
        }


        //******************************** print netDipole ******************************** //
//        boolean printDipole = false;
//        if (printDipole) {
//            if (doVSumMI) {
//                double JEMUEx = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(13);
//                double JEMUExError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(13);
//                double JEMUExCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(13);
//                System.out.println("JEMUEx:\t" + (JEMUEx / numberMolecules)
//                        + " Err:\t" + (JEMUExError / numberMolecules) + " Cor:\t " + JEMUExCor
//                        + " Difficulty:\t" + (JEMUExError * Math.sqrt(totalTime) / nCells / nCells));
//            }
//
//            if (doVSumMI) {
//                double Mx = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(19);
//                double MxError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(19);
//                double MxCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(19);
//                System.out.println("Mx:\t" + (Mx / numberMolecules)
//                        + " Err:\t" + (MxError / numberMolecules) + " Cor:\t " + MxCor
//                        + " Difficulty:\t" + (MxError * Math.sqrt(totalTime) / nCells / nCells));
//            }
//
//            if (doVSumMI) {
//                double JEMUEy = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(14);
//                double JEMUEyError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(14);
//                double JEMUEyCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(14);
//                System.out.println("JEMUEy:\t" + (JEMUEy / numberMolecules)
//                        + " Err:\t" + (JEMUEyError / numberMolecules) + " Cor:\t " + JEMUEyCor
//                        + " Difficulty:\t" + (JEMUEyError * Math.sqrt(totalTime) / nCells / nCells));
//            }
//
//
//            if (doVSumMI) {
//                double My = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(20);
//                double MyError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(20);
//                double MyCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(20);
//                System.out.println("My:\t" + (My / numberMolecules)
//                        + " Err:\t" + (MyError / numberMolecules) + " Cor:\t " + MyCor
//                        + " Difficulty:\t" + (MyError * Math.sqrt(totalTime) / nCells / nCells));
//            }
//
//        }

        sim.integrator.getPotentialEnergy();
        endTime = System.currentTimeMillis();
        System.out.println("Total_Time: " + (endTime - startTime) / (1000.0 * 60.0));
    }

    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
        public boolean mSquare = false;
        public boolean aEE = false;
        public boolean doPair = false;
        public boolean doIdeal = false;
        public boolean doVSum = false;
        public boolean doVSumMI = false;
        public boolean doConventionalE = false;
        public boolean doMappingE = false;
        public boolean doCV = false;
        public boolean doCorrelation = false;
        public boolean doDipole = false;
        public boolean doEnergyMF = false;
        public boolean doAEEMF = false;
        public int formula = 3;
        public boolean doGraphic = false;
        public double temperature = 1;
        public double interactionS = 1;
        public double dipoleMagnitude = 1;
        public int nCells = 5;//number of atoms is nCells*nCells
        public int steps = 1000000;
        public int nMax = 1;
    }
}
