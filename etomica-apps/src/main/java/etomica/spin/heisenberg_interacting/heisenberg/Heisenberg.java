/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicSimple;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.site.PotentialMasterSite;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomNumberGenerator;

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
    public final ActivityIntegrate activityIntegrate;
    public PotentialMasterSite potentialMaster; // difference betweet Pmaster pmastersite
    public Box box;
    public SpeciesSpheresMono spins;
    public P2Spin potential;
    public P1MagneticField field;
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

        potentialMaster = new PotentialMasterSite(this, nCells, space);
        box = new Box(space);
        addBox(box);
        int numAtoms = space.powerD(nCells);

        spins = new SpeciesSpheresRotating(space, new ElementSimple("A"));

        addSpecies(spins);
        box.setNMolecules(spins, numAtoms);
        box.getBoundary().setBoxSize(new Vector2D(nCells, nCells));
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(space, 1), space);
        config.initializeCoordinates(box);

        potential = new P2Spin(space, interactionS);
        // the electric field is zero for now, maybe need to add electric field in the future.
        field = new P1MagneticField(space, dipoleMagnitude);
        integrator = new IntegratorMC(this, potentialMaster);
        mcMove = new MCMoveRotate(potentialMaster, random, space);
        integrator.getMoveManager().addMCMove(mcMove);
        integrator.setTemperature(temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        AtomType type = spins.getLeafType();
//        potentialMaster.addPotential(field, new IAtomType[] {type});
        potentialMaster.addPotential(potential, new AtomType[]{type, type});
        integrator.setBox(box);

    }

    public static void main(String[] args) {
        Param params = new Param();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {

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


        int steps = params.steps;
        int nMax = params.nMax;
        double interactionS = params.interactionS;
        double dipoleMagnitude = params.dipoleMagnitude;

        System.out.println("numberMolecules= " + numberMolecules + " steps= " + steps + " nMax= " + nMax
                + " \ntemperature= " + temperature + " interacitonS= " + interactionS + " dipoleStrength= " + dipoleMagnitude);

        Space sp = Space2D.getInstance();
        Heisenberg sim = new Heisenberg(sp, nCells, temperature, interactionS, dipoleMagnitude);

        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "XY", numberMolecules, sim.space, sim.getController());
            ((DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.getSpecies(0).getAtomType(0), 0.4);
            ((DisplayBoxCanvas2D) simGraphic.getDisplayBox(sim.box).canvas).setOrientationColor(Color.BLUE);
            ((DisplayBoxCanvas2D) simGraphic.getDisplayBox(sim.box).canvas).setOrientationLength(2.5);
            simGraphic.makeAndDisplayFrame();
            return;
        }


        sim.activityIntegrate.setMaxSteps(steps / 5);
        sim.getController().actionPerformed();
        sim.getController().reset();
        int blockNumber = 100;

        //        System.out.println("equilibration finished");
        long equilibrationTime = System.currentTimeMillis();
//        System.out.println("equilibrationTime: " + (equilibrationTime - startTime) / (1000.0 * 60.0));

        int sampleAtInterval = numberMolecules;
        int samplePerBlock = steps / sampleAtInterval / blockNumber;


        MeterSpinMSquare meterMSquare = null;
        AccumulatorAverage dipoleSumSquaredAccumulator = null;

        MeterPotentialEnergy PE = null;
        AccumulatorAverageFixed PEAccumulator = null;

        //mSquare
        if (mSquare) {
            meterMSquare = new MeterSpinMSquare(sim.space, sim.box, dipoleMagnitude);
            dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPump dipolePump = new DataPump(meterMSquare, dipoleSumSquaredAccumulator);
            IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
            dipoleListener.setInterval(sampleAtInterval);
            sim.integrator.getEventManager().addListener(dipoleListener);

            PE = new MeterPotentialEnergy(sim.potentialMaster);
            PE.setBox(sim.box);
            PEAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPumpListener PEPump = new DataPumpListener(PE, PEAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(PEPump);
        }


        //AEE
        MeterMappedAveragingVSum AEEMeter = null;
        AccumulatorAverageCovariance AEEAccumulator = null;
        if (aEE) {
            AEEMeter = new MeterMappedAveragingVSum(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster, doIdeal, doPair, doVSum, doVSumMI, nMax);
            AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPumpListener AEEListener = new DataPumpListener(AEEMeter, AEEAccumulator, sampleAtInterval);
            sim.integrator.getEventManager().addListener(AEEListener);
        }

        sim.activityIntegrate.setMaxSteps(steps);
        sim.getController().actionPerformed();


        //******************************** simulation start ******************************** //
        double dipoleSumSquared = 0;
        double dipoleSumSquaredERR = 0;
        double dipoleSumCor = 0;
        if (mSquare) {
            dipoleSumSquared = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
            dipoleSumSquaredERR = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(AccumulatorAverage.ERROR.index)).x;
            dipoleSumCor = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index)).x;
        }

        double U_ConMC = 0;
        double U_ConMCERR = 0;
        double U_ConMCCor = 0;
        if (mSquare) {
            U_ConMC = ((DataDouble) ((DataGroup) PEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
            U_ConMCERR = ((DataDouble) ((DataGroup) PEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index)).x;
            U_ConMCCor = ((DataDouble) ((DataGroup) PEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index)).x;
        }


        double AEE = 0;
        double AEEER = 0;
        double AEECor = 0;
        if (aEE) {
            AEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            AEEER = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            AEECor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(0);

        }
        long endTime = System.currentTimeMillis();

        double totalTime = (endTime - startTime) / (1000.0 * 60.0);

        double sumIdeal = 0;
        double errSumIdeal = 0;
        double sumIdealCor = 0;
        if (doIdeal) {
            sumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
            errSumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
            sumIdealCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(1);
        }


        double AEEVSum = 0;
        double AEEERVSum = 0;
        double AEECorVSum = 0;
        if (doVSum) {
            AEEVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(5);
            AEEERVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(5);
            AEECorVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(5);
        }

        double AEEVSumMinusIdeal = 0;
        double AEEERVSumMinusIdeal = 0;
        double AEECorVSumMinusIdeal = 0;
        if (doVSumMI) {
            AEEVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(12);
            AEEERVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(12);
            AEECorVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(12);
        }
        if (mSquare) {
            System.out.println("Conventional:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
                    + " Err:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules) + " Cor:\t " + dipoleSumCor
                    + " Difficulty:\t" + (dipoleSumSquaredERR * Math.sqrt(totalTime) / nCells / nCells));
        }

        if (doIdeal && !doPair) {//You need to uncomment this line if you make a jar
            //TODO
//        if (doIdeal) {
            System.out.println("IdealMapping:\t" + (sumIdeal / numberMolecules)
                    + " Err:\t" + (errSumIdeal / numberMolecules) + " Cor:\t " + sumIdealCor
                    + " Difficulty:\t" + (errSumIdeal * Math.sqrt(totalTime) / nCells / nCells));
        }
        if (doPair) {
            System.out.println("PairMapping:\t" + (AEE / numberMolecules)
                    + " Err:\t" + (AEEER / numberMolecules) + " Cor:\t " + AEECor
                    + " Difficulty:\t" + (AEEER * Math.sqrt(totalTime) / nCells / nCells));
        }
        if (doVSum) {
            System.out.println("mappingvSum:\t" + (AEEVSum / numberMolecules)
                    + " Err:\t" + (AEEERVSum / numberMolecules) + " Cor:\t " + AEECorVSum
                    + " Difficulty:\t" + (AEEERVSum * Math.sqrt(totalTime) / nCells / nCells));
        }
        if (doVSumMI) {
            System.out.println("mappingvSumMI:\t" + (AEEVSumMinusIdeal / numberMolecules)
                    + " Err:\t" + (AEEERVSumMinusIdeal / numberMolecules) + " Cor:\t " + AEECorVSumMinusIdeal
                    + " Difficulty:\t" + (AEEERVSumMinusIdeal * Math.sqrt(totalTime) / nCells / nCells));
        }


        double U_Map = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(21);
        double U_MapError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(21);
        double U_MapCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(21);


        if (doVSumMI) {
            System.out.println("U_Map:\t" + (U_Map / numberMolecules)
                    + " Err:\t" + (U_MapError / numberMolecules) + " Cor:\t " + U_MapCor
                    + " Difficulty:\t" + (U_MapError * Math.sqrt(totalTime) / nCells / nCells));
        }

        double U_Con = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(23);
        double U_ConError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(23);
        double U_ConCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(23);


        if (doVSumMI) {
            System.out.println("U_Con:\t" + (U_Con / numberMolecules)
                    + " Err:\t" + (U_ConError / numberMolecules) + " Cor:\t " + U_ConCor
                    + " Difficulty:\t" + (U_ConError * Math.sqrt(totalTime) / nCells / nCells));
        }


        if (doVSumMI) {
            System.out.println("U_ConMC:\t" + (U_ConMC / numberMolecules)
                    + " Err:\t" + (U_ConMCERR / numberMolecules) + " Cor:\t " + U_ConMCCor
                    + " Difficulty:\t" + (U_ConMCERR * Math.sqrt(totalTime) / nCells / nCells));
        }


        double U_MapSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(22);
        double U_MapSquareError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(22);
        double U_MapSquareCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(22);

        if (doVSumMI) {
            System.out.println("USquare_Map:\t" + (U_MapSquare / numberMolecules)
                    + " Err:\t" + (U_MapSquareError / numberMolecules) + " Cor:\t " + U_MapSquareCor
                    + " Difficulty:\t" + (U_MapSquareError * Math.sqrt(totalTime) / nCells / nCells));
        }


        double U_ConSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(24);
        double U_ConSquareError = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(24);
        double U_ConSquareCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(24);

        if (doVSumMI) {
            System.out.println("USquare_Con:\t" + (U_ConSquare / numberMolecules)
                    + " Err:\t" + (U_ConSquareError / numberMolecules) + " Cor:\t " + U_ConSquareCor
                    + " Difficulty:\t" + (U_ConSquareError * Math.sqrt(totalTime) / nCells / nCells));
        }

        IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        double CV_Map = U_MapSquare - U_Map * U_Map;
        double CV_MapErr = Math.sqrt(U_MapSquareError * U_MapSquareError
                + 4 * U_Map * U_Map * U_MapError * U_MapError
                - 4 * U_MapSquareError * U_Map  * U_MapError * covariance.getValue(21 * 25 + 22) / Math.sqrt(covariance.getValue(21 * 25 + 21) * covariance.getValue(22 * 25 + 22))
        );

        if (doVSumMI) {
            System.out.println("CV_Map:\t" + (CV_Map / numberMolecules)
                    + " Err:\t" + (CV_MapErr / numberMolecules)) ;
        }




        double CV_Con = U_ConSquare - U_Con * U_Con;
        double CV_ConErr = Math.sqrt(U_ConSquareError * U_ConSquareError
                + 4 * U_Con * U_Con * U_ConError * U_ConError
                - 4 * U_ConSquareError * U_Con  * U_ConError * covariance.getValue(23 * 25 + 24) / Math.sqrt(covariance.getValue(23 * 25 + 23) * covariance.getValue(24 * 25 + 24))
        );

        if (doVSumMI) {
            System.out.println("CV_Con:\t" + (CV_Con / numberMolecules)
                    + " Err:\t" + (CV_ConErr / numberMolecules)) ;
        }




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
        public boolean mSquare = true;
        public boolean aEE = true;
        public boolean doPair = false;
        public boolean doIdeal = true;
        public boolean doVSum = false;
        public boolean doVSumMI = true;
        public double temperature = 10;//Kelvin
        public double interactionS = 1;
        public double dipoleMagnitude = 1;
        public int nCells = 3;//number of atoms is nCells*nCells
        public int steps = 500000;
        public int nMax = 5;
    }
}
