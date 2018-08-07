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
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.site.PotentialMasterSite;
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
 *
 *
 * @author Weisong Lin
 */
public class HeisenbergPair extends Simulation {

    private static final String APP_NAME = "Heisenberg";
    public final ActivityIntegrate activityIntegrate;
    public Box box;
    public SpeciesSpheresMono spins;
    public P2Spin potential;
    public P1MagneticField field;
    public MCMoveRotatePair mcMove;
    private IntegratorMC integrator;

    /**
     * 2D heisenberg model in square lattice.
     *
     * @param space           use to define vector/tensor
     * @param interactionS    the J in heisenberg energy function: U = J*Cos(theta1-theta2)
     * @param dipoleMagnitude is the strength of heisenberg dipole.
     */
    public HeisenbergPair(Space space, double temperature, double interactionS, double dipoleMagnitude) {
        super(Space2D.getInstance());
     //   setRandom(new RandomNumberGenerator(1)); //debug only
        System.out.println("============================the RandomSeed is one ===========================");

        box = new Box(space);
        addBox(box);
        int numAtoms = 2;

        spins = new SpeciesSpheresRotating(space, new ElementSimple("A"));

        addSpecies(spins);
        box.setNMolecules(spins, numAtoms);

        potential = new P2Spin(space, interactionS);
        // the electric field is zero for now, maybe need to add electric field in the future.
        field = new P1MagneticField(space, dipoleMagnitude);
        integrator = new IntegratorMC(this, null);
        mcMove = new MCMoveRotatePair(potential, random, space);
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

        boolean isGraphic = params.isGraphic;
        double temperature = params.temperature;
        int numberMolecules = 2;
        boolean aEE = params.aEE;
        boolean mSquare = params.mSquare;
        int steps = params.steps;
        double interactionS = params.interactionS;
        double dipoleMagnitude = params.dipoleMagnitude;

        System.out.println("steps= " + steps);
        System.out.println("numberMolecules= " + numberMolecules);
        System.out.println("temperature= " + temperature);

        Space sp = Space2D.getInstance();
        HeisenbergPair sim = new HeisenbergPair(sp, temperature, interactionS, dipoleMagnitude);

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
        System.out.println("number of molecules are: " + 2);
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
        MeterMappedAveragingPair AEEMeter = null;
        AccumulatorAverageCovariance AEEAccumulator = null;
        if (aEE) {
            AEEMeter = new MeterMappedAveragingPair(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potential);
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
        double avgdipole=0;
        double errdipole=0;
        double avgJEEMJEJE=0;
        double errJEEMJEJE=0;
        double avgUEE=0;
        double errUEE=0;
        double avgJEMUE=0;
        double errJEMUE=0;
        double avgJEMUEsquare=0;
        double errJEMUEsquare=0;
        double avgdipoleconv=0;
        double errdipoleconv=0;

        if (aEE) {

//            AEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
//            AEEER = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
//            AEECor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(0);
//            IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
//            covariance.getValue(1);


            avgdipole = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(2);
            errdipole = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(2);
            avgJEMUE=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(4);
            errJEMUE=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(4);
            avgJEMUEsquare=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(3);
            errJEMUEsquare=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(3);
            avgJEEMJEJE=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            errJEEMJEJE=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            avgUEE=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
            errUEE=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
            avgdipoleconv=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(5);
            errdipoleconv=((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(5);
     //       IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
     //       covariance.getValue(1);
            AEE = -avgJEEMJEJE+avgUEE-avgJEMUEsquare+(avgJEMUE*avgJEMUE);
    //        AEEER = Math.sqrt(ERsum0 * ERsum0 + 4 * sum1 * sum1 * ERsum1 * ERsum1 -
    //                2 * ERsum0 * sum1 * 2 * ERsum1 * covariance.getValue(1) / Math.sqrt(covariance.getValue(0) * covariance.getValue(8)));
    //        AEECor = covariance.getValue(1) / Math.sqrt(covariance.getValue(0) * covariance.getValue(3));

        }

        long endTime = System.currentTimeMillis();

        double totalTime = (endTime - startTime) / (1000.0 * 60.0);
        if (mSquare) {
            System.out.println("-<M^2>*bt*bt:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
                    + " mSquareErr:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules)
                    + " mSquareDifficulty:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules) * Math.sqrt(totalTime)
                    + " dipolesumCor= " + dipoleSumCor);
            System.out.println("mSquare_Time: " + (endTime - startTime) / (1000.0 * 60.0));
        }

        if (aEE) {
            System.out.println("AEE_new:\t" + (AEE / numberMolecules));
   //                 + " AEEErr:\t" + (AEEER / numberMolecules)
     //               + " AEEDifficulty:\t" + (AEEER * Math.sqrt(totalTime) / numberMolecules)
     //               + " AEECor= " + AEECor);
            System.out.println("AEE_Time: " + (endTime - startTime) / (1000.0 * 60.0));
            System.out.println("avgdipole: " + avgdipole);
           System.out.println("JEEMJEJE: " + avgJEEMJEJE +"  errJEEMJEJE: " + errJEEMJEJE +"  UEE: " + avgUEE +"  errUEE: " + errUEE +"  JEMUE: " + avgJEMUE +"  errJEMUE: " + errJEMUE + "  JEMUEsq: " + avgJEMUEsquare +"  errJEMUEsq: " + errJEMUEsquare+"  dipoleconv: " + avgdipoleconv+"  errdipoleconv: " + errdipoleconv+"  dipolemap: " + avgdipole+"  errdipolemap: " + errdipole);


        }


    }

    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
        public boolean isGraphic = false;
        public boolean mSquare = true;
        public boolean aEE = true;
        public double temperature = 0.8;// Kelvin
        public double interactionS = 2;
        public double dipoleMagnitude = 10.0;
        public int steps = 15000;
    }
}
