/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.listener.IntegratorListenerAction;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomNumberGenerator;

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
public class Heisenberg3Pair extends Simulation {

    private static final String APP_NAME = "Heisenberg";
    public final ActivityIntegrate activityIntegrate;
    public Box box;
    public SpeciesSpheresMono spins;
    public P2Spin potential;
    public P1MagneticField field;
    public MCMoveRotate3Pair mcMove;
    private IntegratorMC integrator;

    /**
     * 2D heisenberg model in square lattice.
     *
     * @param space           use to define vector/tensor
     * @param interactionS    the J in heisenberg energy function: U = J*Cos(theta1-theta2)
     * @param dipoleMagnitude is the strength of heisenberg dipole.
     */
    public Heisenberg3Pair(Space space, double temperature, double interactionS, double dipoleMagnitude, int numberMolecules) {
        super(Space2D.getInstance());
//        setRandom(new RandomNumberGenerator(1)); //debug only
//        System.out.println("============================the RandomSeed is one ===========================");

        box = new Box(space);
        addBox(box);
        int numAtoms = numberMolecules;

        spins = new SpeciesSpheresRotating(space, new ElementSimple("A"));

        addSpecies(spins);
        box.setNMolecules(spins, numAtoms);

        potential = new P2Spin(space, interactionS);
        // the electric field is zero for now, maybe need to add electric field in the future.
        field = new P1MagneticField(space, dipoleMagnitude);
        integrator = new IntegratorMC(this, null);
        mcMove = new MCMoveRotate3Pair(potential, random, space);
        integrator.getMoveManager().addMCMove(mcMove);
        integrator.setTemperature(temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        AtomType type = spins.getLeafType();
//        potentialMaster.addPotential(field, new IAtomType[] {type});

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
        int numberMolecules = params.numberMolecules;
        boolean aEE = params.aEE;
        boolean mSquare = params.mSquare;
        int steps = params.steps;
        double interactionS = params.interactionS;
        double dipoleMagnitude = params.dipoleMagnitude;

        System.out.println("steps= " + steps);
        System.out.println("numberMolecules= " + numberMolecules);
        System.out.println("temperature= " + temperature);

        Space sp = Space2D.getInstance();
        Heisenberg3Pair sim = new Heisenberg3Pair(sp, temperature, interactionS, dipoleMagnitude, numberMolecules);

        MeterSpinMSquare meterMSquare = null;
        AccumulatorAverage dipoleSumSquaredAccumulator = null;

        sim.activityIntegrate.setMaxSteps(steps / 5);
        sim.getController().actionPerformed();
        sim.getController().reset();
        int blockNumber = 100;


        int sampleAtInterval = numberMolecules;
        int samplePerBlock = steps / sampleAtInterval / blockNumber;
        System.out.println("number of blocks is : " + blockNumber);
        System.out.println("sample per block is : " + samplePerBlock);
        System.out.println("number of molecules are: " + numberMolecules);
        System.out.println("interacitonS= " + interactionS);
        System.out.println("dipoleStrength= " + dipoleMagnitude);

        System.out.println("equilibration finished");
        long equilibrationTime = System.currentTimeMillis();
        System.out.println("equilibrationTime: " + (equilibrationTime - startTime) / (1000.0 * 60.0));

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
        MeterMappedAveragingVSum3Pair AEEMeter = null;
        AccumulatorAverageCovariance AEEAccumulator = null;
        if (aEE) {
            AEEMeter = new MeterMappedAveragingVSum3Pair(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potential);
            AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPump AEEPump = new DataPump(AEEMeter, AEEAccumulator);
            IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
            AEEListener.setInterval(sampleAtInterval);
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
//        if (mSquare) {
//            System.out.println("-<M^2>*bt*bt:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
//                    + " mSquareErr:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules)
//                    + " mSquareDifficulty:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules) * Math.sqrt(totalTime)
//                    + " dipolesumCor= " + dipoleSumCor);
//            System.out.println("mSquare_Time: " + (endTime - startTime) / (1000.0 * 60.0));
//        }
//        if (aEE) {
//            System.out.println("AEE_new:\t" + (AEE / numberMolecules)
//                    + " AEEErr:\t" + (AEEER / numberMolecules)
//                    + " AEEDifficulty:\t" + (AEEER * Math.sqrt(totalTime) / numberMolecules)
//                    + " AEECor= " + AEECor);
//            System.out.println("AEE_Time: " + (endTime - startTime) / (1000.0 * 60.0));
//        }

        double AEEDirect = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(10);
        double errAEEDirect = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(10);

        double AEEDirectJ0 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(13);
        double errAEEDirectJ0 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(13);

        System.out.println(" AEEDirect= " + AEEDirect / numberMolecules + " Err " + errAEEDirect / numberMolecules);
        System.out.println(" AEEDirectJ0= " + AEEDirectJ0 / numberMolecules + " Err " + errAEEDirectJ0 / numberMolecules);

        double JEMUEx = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
        double errJEMUEx = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
        double JEMUEy = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(2);
        double errJEMUEy = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(2);
        double JEMUExSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(8);
        double errJEMUExSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(8);
        double JEMUEySquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(9);
        double errJEMUEySquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(9);
        System.out.println("JEMUEx= " + JEMUEx / numberMolecules + " Err " + errJEMUEx / numberMolecules);
        System.out.println("JEMUEy= " + JEMUEy / numberMolecules + " Err " + errJEMUEy / numberMolecules);
        System.out.println("JEMUExSquare= " + JEMUExSquare / numberMolecules + " Err " + errJEMUExSquare / numberMolecules);
        System.out.println("JEMUEySquare= " + JEMUEySquare / numberMolecules + " Err " + errJEMUEySquare / numberMolecules);


        double JEMUExIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(3);
        double errJEMUExIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(3);
        double JEMUEyIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(4);
        double errJEMUEyIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(4);
        double JEMUExSquareIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(11);
        double errJEMUExSquareIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(11);
        double JEMUEySquareIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(12);
        double errJEMUEySquareIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(12);
        System.out.println("JEMUExIdeal= " + JEMUExIdeal / numberMolecules + " Err " + errJEMUExIdeal / numberMolecules);
        System.out.println("JEMUEyIdeal= " + JEMUEyIdeal / numberMolecules + " Err " + errJEMUEyIdeal / numberMolecules);
        System.out.println("JEMUExSquareIdeal= " + JEMUExSquareIdeal / numberMolecules + " Err " + errJEMUExSquareIdeal / numberMolecules);
        System.out.println("JEMUEySquareIdeal= " + JEMUEySquareIdeal / numberMolecules + " Err " + errJEMUEySquareIdeal / numberMolecules);


//        System.out.println("Conventional:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
//                + " Err:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules));

        double sumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(7);
        double errSumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(7);
        double sumIdealCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(7);

//        System.out.println("IdealMapping:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumIdeal / numberMolecules)
//                + " Err:\t" + (errSumIdeal / numberMolecules));
//
//        System.out.println("mapping:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEE / numberMolecules)
//                + " Err:\t" + (AEEER / numberMolecules));

        double sumJEEMJEJE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(5);
        double errJEEMJEJE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(5);
        System.out.println("JEEMJEJE:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumJEEMJEJE / numberMolecules)
                + " Err:\t" + (errJEEMJEJE / numberMolecules));

        double sumUEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(6);
        double errUEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(6);
        System.out.println("UEE:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumUEE / numberMolecules)
                + " Err:\t" + (errUEE / numberMolecules));


        double AEEVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(17);
        double AEEERVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(17);
        double AEECorVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(17);


        double AEEVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(24);
        double AEEERVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(24);
        double AEECorVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(24);


        System.out.println("Conventional:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
                + " Err:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules) + " Cor:\t " + dipoleSumCor);

        System.out.println("IdealMapping:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumIdeal / numberMolecules)
                + " Err:\t" + (errSumIdeal / numberMolecules) + " Cor:\t " + sumIdealCor);

        System.out.println("mapping:\t\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEE / numberMolecules)
                + " Err:\t" + (AEEER / numberMolecules) + " Cor:\t " + AEECor);

        System.out.println("mappingvSum:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEEVSum / numberMolecules)
                + " Err:\t" + (AEEERVSum / numberMolecules) + " Cor:\t " + AEECorVSum);

        System.out.println("vSumMinusIdeal:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEEVSumMinusIdeal / numberMolecules)
                + " Err:\t" + (AEEERVSumMinusIdeal / numberMolecules) +" Cor:\t " + AEECorVSumMinusIdeal);
//

        System.out.println("These are the result of VSum");
        double sumJEMUExV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(18);
        double errJEMUExV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(18);
        System.out.println("JEMUExV:\t" + " Value:\t" + (sumJEMUExV / numberMolecules)
                + " Err:\t" + (errJEMUExV / numberMolecules));

        double sumJEMUEyV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(19);
        double errJEMUEyV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(19);
        System.out.println("JEMUEyV:\t" + " Value:\t" + (sumJEMUEyV / numberMolecules)
                + " Err:\t" + (errJEMUEyV / numberMolecules));


        double sumJEMUExVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(22);
        double errJEMUExVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(22);
        System.out.println("JEMUExVSquare:\t" + " Value:\t" + (sumJEMUExVSquare / numberMolecules)
                + " Err:\t" + (errJEMUExVSquare / numberMolecules));

        double sumJEMUEyVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(23);
        double errJEMUEyVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(23);
        System.out.println("JEMUEyVSquare:\t" + " Value:\t" + (sumJEMUEyVSquare / numberMolecules)
                + " Err:\t" + (errJEMUEyVSquare / numberMolecules));



        double sumJEEMJEJEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(20);
        double errJEEMJEJEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(20);
        System.out.println("JEEMJEJEV:\t" + " Value:\t" + (sumJEEMJEJEV / numberMolecules)
                + " Err:\t" + (errJEEMJEJEV / numberMolecules));


        double sumUEEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(21);
        double errUEEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(21);
        System.out.println("UEEV:\t" + " Value:\t" + (sumUEEV / numberMolecules)
                + " Err:\t" + (errUEEV / numberMolecules));


        //VSum Minus Ideal
        double sumJEMUExVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(25);
        double errJEMUExVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(25);
        System.out.println("JEMUExVMI:\t" + " Value:\t" + (sumJEMUExVMI / numberMolecules)
                + " Err:\t" + (errJEMUExVMI / numberMolecules));

        double sumJEMUEyVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(26);
        double errJEMUEyVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(26);
        System.out.println("JEMUEyVMI:\t" + " Value:\t" + (sumJEMUEyVMI / numberMolecules)
                + " Err:\t" + (errJEMUEyVMI / numberMolecules));


        double sumJEMUExVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(29);
        double errJEMUExVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(29);
        System.out.println("JEMUExVSquareMI:\t" + " Value:\t" + (sumJEMUExVSquareMI / numberMolecules)
                + " Err:\t" + (errJEMUExVSquareMI / numberMolecules));

        double sumJEMUEyVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(30);
        double errJEMUEyVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(30);
        System.out.println("JEMUEyVSquareMI:\t" + " Value:\t" + (sumJEMUEyVSquareMI / numberMolecules)
                + " Err:\t" + (errJEMUEyVSquareMI / numberMolecules));


        double sumJEEMJEJEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(27);
        double errJEEMJEJEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(27);
        System.out.println("JEEMJEJEVIMI:\t" + " Value:\t" + (sumJEEMJEJEVMI / numberMolecules)
                + " Err:\t" + (errJEEMJEJEVMI / numberMolecules));


        double sumUEEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(28);
        double errUEEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(28);
        System.out.println("UEEVMI:\t\t" + " Value:\t" + (sumUEEVMI / numberMolecules)
                + " Err:\t" + (errUEEVMI / numberMolecules));


        System.out.println("Time: " + (endTime - startTime) / (1000.0 * 60.0));
    }

    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
        public boolean mSquare = true;
        public boolean aEE = true;
        public double temperature = 1.1;// Kelvin
        public double interactionS = 1.2;
        public double dipoleMagnitude = 1.6;
        public int steps = 100000;
        public int numberMolecules = 3;
    }
}
