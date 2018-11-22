/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.chem.elements.ElementSimple;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.site.Api1ASite;
import etomica.nbr.site.NeighborSiteManager;
import etomica.nbr.site.PotentialMasterSite;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.systems.LJ;
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
        int steps = params.steps;
        double interactionS = params.interactionS;
        double dipoleMagnitude = params.dipoleMagnitude;

        System.out.println("steps= " + steps);
        System.out.println("numberMolecules= " + numberMolecules);
        System.out.println("temperature= " + temperature);

        Space sp = Space2D.getInstance();
        Heisenberg sim = new Heisenberg(sp, nCells, temperature, interactionS, dipoleMagnitude);

        MeterSpinMSquare meterMSquare = null;
        AccumulatorAverage dipoleSumSquaredAccumulator = null;

        sim.activityIntegrate.setMaxSteps(steps / 5);
        sim.getController().actionPerformed();
        sim.getController().reset();
        int blockNumber = 100;


        int sampleAtInterval = numberMolecules;
        int samplePerBlock = steps / sampleAtInterval / blockNumber;
//        if (samplePerBlock == 0) {
//            samplePerBlock = 1;
//        }
        System.out.println("number of blocks is : " + blockNumber);
        System.out.println("sample per block is : " + samplePerBlock);
        System.out.println("number of molecules are: " + nCells * nCells);
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
        MeterMappedAveraging AEEMeter = null;
        AccumulatorAverageCovariance AEEAccumulator = null;
        if (aEE) {
            AEEMeter = new MeterMappedAveraging(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster);
            AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPump AEEPump = new DataPump(AEEMeter, AEEAccumulator);
            IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
            AEEListener.setInterval(sampleAtInterval);
            sim.integrator.getEventManager().addListener(AEEListener);
        }

        //Mapped averaging using pair excess
        Box box = sim.getBox(0);
//        BoxAgentManager.BoxAgentSource<BoxCellManager> boxAgentSource = new PotentialMasterSite.BoxAgentSiteSource(nCells, sim.space);
//        Api1ASite api = new Api1ASite(sim.getSpace().D(), new BoxAgentManager<BoxCellManager>(boxAgentSource, BoxCellManager.class));
//        api.setDirection(IteratorDirective.Direction.UP);
//        api.setBox(box);
//        IAtomList atoms = box.getLeafList();
//        for(int i=0; i<atoms.getAtomCount(); i++) {
//            IAtom atom = atoms.getAtom(i);
//            api.setTarget(atom);
//            api.reset();
//            for(IAtomList nbr=api.next(); nbr!=null; nbr=api.next()) {
//                System.out.println(nbr);
//            }
//        }
        NeighborListMaker pairList = new NeighborListMaker(nCells * nCells);
        sim.potentialMaster.calculate(box, new IteratorDirective(), pairList);
//        System.out.println("pair count: " + pairList.count);
        MeterMappedAveragingPairExcess[] m2ExList = new MeterMappedAveragingPairExcess[pairList.count];
        AccumulatorAverageCovariance[] m2ExAccumulators = new AccumulatorAverageCovariance[pairList.count];
        for (int i = 0; i < pairList.count; i++) {
            m2ExList[i] = new MeterMappedAveragingPairExcess(pairList.pairs[i], sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster, sim.potential);
            m2ExAccumulators[i] = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPump AEEPump = new DataPump(m2ExList[i], m2ExAccumulators[i]);
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
        double avgdipolex = 0;
        double avgdipoley = 0;
        double errdipoley = 0;

        double errdipolex = 0;
        if (aEE) {
//combined mapping
            double sum0 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            double ERsum0 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            double sum1 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
            double ERsum1 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
            double sum2 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(2);
            double errSum2 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(2);

            double sum3 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(3);
            double ERsum3 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(3);
            double sum4 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(4);
            double ERsum4 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(4);

            IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            covariance.getValue(1);
            AEE = sum0 + sum1 * sum1 + sum2 * sum2 - sum3 * sum3 - sum4 * sum4;//combined mapping
            AEEER = Math.sqrt(ERsum0 * ERsum0 + 4 * sum1 * sum1 * ERsum1 * ERsum1 -
                    2 * ERsum0 * sum1 * 2 * ERsum1 * covariance.getValue(1) / Math.sqrt(covariance.getValue(0) * covariance.getValue(3)));
            AEECor = covariance.getValue(1) / Math.sqrt(covariance.getValue(0) * covariance.getValue(3));

        }

        double aEE2 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(7); //independent-spin mapping
        double aEE2Test = 0;
        double x7 = 0;
        for (int i = 0; i < m2ExAccumulators.length; i++) {
            double m2 = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            double JEMUE2x = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
            double JEMUE2y = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(2);
            double m2Id = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(7);
            aEE2 += m2 + JEMUE2x * JEMUE2x + JEMUE2y * JEMUE2y - m2Id;

            double m2Test = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(6);
            double JEMUE2xIdeal = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(3);
            double JEMUE2yIdeal = ((DataGroup) m2ExAccumulators[i].getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(4);
            aEE2Test += m2Test + JEMUE2xIdeal * JEMUE2xIdeal + JEMUE2yIdeal * JEMUE2yIdeal;
            x7 += m2Id;
        }
        System.out.println("x7 = " + x7 / numberMolecules + " aEE2Test =" + aEE2Test / numberMolecules);
        long endTime = System.currentTimeMillis();

        double totalTime = (endTime - startTime) / (1000.0 * 60.0);
//        if (mSquare) {
//            System.out.println("-<M^2>*bt*bt:\t" + (-dipoleSumSquared / temperature / temperature / nCells / nCells)
//                    + " mSquareErr:\t" + (dipoleSumSquaredERR / temperature / temperature / nCells / nCells)
//                    + " mSquareDifficulty:\t" + (dipoleSumSquaredERR / temperature / temperature / nCells / nCells) * Math.sqrt(totalTime)
//                    + " dipolesumCor= " + dipoleSumCor);
//            System.out.println("mSquare_Time: " + (endTime - startTime) / (1000.0 * 60.0));
//        }
//
//        if (aEE) {
//            System.out.println("AEE_new:\t" + (AEE / nCells / nCells)
//                    + " AEEErr:\t" + (AEEER / nCells / nCells)
//                    + " AEEDifficulty:\t" + (AEEER * Math.sqrt(totalTime) / nCells / nCells)
//                    + " AEECor= " + AEECor);
//            System.out.println("AEE_Time: " + (endTime - startTime) / (1000.0 * 60.0));
////            System.out.println("avgdipolex: " + avgdipolex + "errdipolex: " + errdipolex);
////            System.out.println("avgdipoley: " + avgdipoley + "errdipoley: " + errdipoley);
//        }


        System.out.println("Conventional:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
                + " Err:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules));

        double sumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(7);
        double errSumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(7);
        System.out.println("IdealMapping:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumIdeal / numberMolecules)
                + " Err:\t" + (errSumIdeal / numberMolecules));

        System.out.println("Mapping:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEE / numberMolecules)
                + " Err:\t" + (AEEER / numberMolecules));
        System.out.println("Mapping2:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (aEE2 / numberMolecules));

        double sumJEEMJEJE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(5);
        double errJEEMJEJE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(5);
        System.out.println("JEEMJEJE:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumJEEMJEJE / numberMolecules)
                + " Err:\t" + (errJEEMJEJE / numberMolecules));

        double sumUEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(6);
        double errUEE = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(6);
        System.out.println("UEE:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumUEE / numberMolecules)
                + " Err:\t" + (errUEE / numberMolecules));


        System.out.println("Time: " + (endTime - startTime) / (1000.0 * 60.0));
        System.out.println(" dipolesumCor " + dipoleSumCor);

    }


    static class NeighborListMaker implements PotentialCalculation {
        int count = 0;
        final AtomPair[] pairs;
        public NeighborListMaker(int nAtoms) {
            pairs = new AtomPair[nAtoms * 2];
        }

        public void doCalculation(IAtomList atoms, IPotentialAtomic dummy) {
            System.out.println(atoms);
            pairs[count] = new AtomPair(atoms.getAtom(0), atoms.getAtom(1));
            count++;
        }

    }

    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
        public boolean mSquare = true;
        public boolean aEE = true;
        public double temperature = 1.0;// Kelvin
        public int nCells = 3;//number of atoms is nCells*nCells
        public double interactionS = 1.0;
        public double dipoleMagnitude = 1.0;
        public int steps = 5000;
    }
}
