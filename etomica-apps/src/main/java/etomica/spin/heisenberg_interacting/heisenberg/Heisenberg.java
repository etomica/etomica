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
        MeterMappedAveragingVSum AEEMeter = null;
        AccumulatorAverageCovariance AEEAccumulator = null;
        if (aEE) {
            AEEMeter = new MeterMappedAveragingVSum(sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster);
            AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPump AEEPump = new DataPump(AEEMeter, AEEAccumulator);
            IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
            AEEListener.setInterval(sampleAtInterval);
            sim.integrator.getEventManager().addListener(AEEListener);
        }

        //Mapped averaging using pair excess
//        Box box = sim.getBox(0);
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
//        NeighborListMaker pairList = new NeighborListMaker(nCells * nCells);
//        sim.potentialMaster.calculate(box, new IteratorDirective(), pairList);
////        System.out.println("pair count: " + pairList.count);
//        MeterMappedAveragingPairExcess[] m2ExList = new MeterMappedAveragingPairExcess[pairList.count];
//        AccumulatorAverageCovariance[] m2ExAccumulators = new AccumulatorAverageCovariance[pairList.count];
//        for (int i = 0; i < pairList.count; i++) {
//            m2ExList[i] = new MeterMappedAveragingPairExcess(pairList.pairs[i], sim.space, sim.box, sim, temperature, interactionS, dipoleMagnitude, sim.potentialMaster, sim.potential);
//            m2ExAccumulators[i] = new AccumulatorAverageCovariance(samplePerBlock, true);
//            DataPump AEEPump = new DataPump(m2ExList[i], m2ExAccumulators[i]);
//            IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
//            AEEListener.setInterval(sampleAtInterval);
//            sim.integrator.getEventManager().addListener(AEEListener);
//        }

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
        if (mSquare) {
            System.out.println("-<M^2>*bt*bt:\t" + (-dipoleSumSquared / temperature / temperature / nCells / nCells)
                    + " mSquareErr:\t" + (dipoleSumSquaredERR / temperature / temperature / nCells / nCells)
                    + " mSquareDifficulty:\t" + (dipoleSumSquaredERR / temperature / temperature / nCells / nCells) * Math.sqrt(totalTime)
                    + " dipolesumCor= " + dipoleSumCor);
//            System.out.println("mSquare_Time: " + (endTime - startTime) / (1000.0 * 60.0));
        }

        if (aEE) {
            System.out.println("AEE_new:\t" + (AEE / nCells / nCells)
                    + " AEEErr:\t" + (AEEER / nCells / nCells)
                    + " AEEDifficulty:\t" + (AEEER * Math.sqrt(totalTime) / nCells / nCells)
                    + " AEECor= " + AEECor);
//            System.out.println("AEE_Time: " + (endTime - startTime) / (1000.0 * 60.0));
        }

        double AEEVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(5);
        double AEEERVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(5);
        double AEECorVSum = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(5);
//        System.out.println("AEE_VSUM:\t" + (AEEVSum / nCells / nCells)
//                + " AEEVSumErr:\t" + (AEEERVSum / nCells / nCells)
//                + " AEEVSumDifficulty:\t" + (AEEERVSum * Math.sqrt(totalTime) / nCells / nCells)
//                + " AEEVSumCor= " + AEECorVSum);

        double AEEVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(12);
        double AEEERVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(12);
        double AEECorVSumMinusIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(12);




        double sumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(1);
        double errSumIdeal = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(1);
        double sumIdealCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(1);
//        System.out.println("IdealMapping:\t" + (sumIdeal / numberMolecules)
//                + " Err:\t" + (errSumIdeal / numberMolecules)
//                + " IdealDifficulty:\t" + (errSumIdeal * Math.sqrt(totalTime) / nCells / nCells)
//                + " Cor:\t" + sumIdealCor);

//        double netDipolex = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(2);
//        double netDipolexErr = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(2);
//        double netDipolexCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(2);
//        System.out.println("netDipolex:\t" + (netDipolex / numberMolecules)
//                + " Err:\t" + (netDipolexErr / numberMolecules)
//                + " Cor:\t" + netDipolexCor);
//
//        double netDipoley = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(3);
//        double netDipoleyErr = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(3);
//        double netDipoleyCor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(3);
//        System.out.println("netDipoley:\t" + (netDipoley / numberMolecules)
//                + " Err:\t" + (netDipoleyErr / numberMolecules)
//                + " Cor:\t" + netDipoleyCor);
//        System.out.println("Simulation time " + (endTime - startTime) / (1000.0 * 60.0) + " mins");


        System.out.println("Conventional:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (-dipoleSumSquared / temperature / temperature / numberMolecules)
                + " Err:\t" + (dipoleSumSquaredERR / temperature / temperature / numberMolecules) + " Cor:\t " + dipoleSumCor);

        System.out.println("IdealMapping:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (sumIdeal / numberMolecules)
                + " Err:\t" + (errSumIdeal / numberMolecules) + " Cor:\t " + sumIdealCor);

        System.out.println("mapping:\t\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEE / numberMolecules)
                + " Err:\t" + (AEEER / numberMolecules) + " Cor:\t " + AEECor);

        System.out.println("mappingvSum:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEEVSum / numberMolecules)
                + " Err:\t" + (AEEERVSum / numberMolecules) + " Cor:\t " + AEECorVSum);

        System.out.println("mappingvSumMI:\t" + "bJ\t" + (interactionS / temperature) + " Value:\t" + (AEEVSumMinusIdeal / numberMolecules)
                + " Err:\t" + (AEEERVSumMinusIdeal / numberMolecules) + " Cor:\t " + AEECorVSumMinusIdeal);

//        System.out.println("These are the result of VSum" );
        double sumJEMUExV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(6);
        double errJEMUExV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(6);
        System.out.println("JEMUExV:\t" + " Value:\t" + (sumJEMUExV / numberMolecules)
                + " Err:\t" + (errJEMUExV / numberMolecules));

        double sumJEMUEyV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(7);
        double errJEMUEyV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(7);
        System.out.println("JEMUEyV:\t" + " Value:\t" + (sumJEMUEyV / numberMolecules)
                + " Err:\t" + (errJEMUEyV / numberMolecules));


        double sumJEMUExVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(10);
        double errJEMUExVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(10);
        System.out.println("JEMUExVSquare:\t" + " Value:\t" + (sumJEMUExVSquare / numberMolecules)
                + " Err:\t" + (errJEMUExVSquare / numberMolecules));

        double sumJEMUEyVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(11);
        double errJEMUEyVSquare = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(11);
        System.out.println("JEMUEyVSquare:\t" + " Value:\t" + (sumJEMUEyVSquare / numberMolecules)
                + " Err:\t" + (errJEMUEyVSquare / numberMolecules));


        double sumJEEMJEJEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(8);
        double errJEEMJEJEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(8);
        System.out.println("JEEMJEJEV:\t" + " Value:\t" + (sumJEEMJEJEV / numberMolecules)
                + " Err:\t" + (errJEEMJEJEV / numberMolecules));


        double sumUEEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(9);
        double errUEEV = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(9);
        System.out.println("UEEV:\t\t" + " Value:\t" + (sumUEEV / numberMolecules)
                + " Err:\t" + (errUEEV / numberMolecules));


        //VSum Minus Ideal

        double sumJEMUExVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(13);
        double errJEMUExVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(13);
        System.out.println("JEMUExVMI:\t" + " Value:\t" + (sumJEMUExVMI / numberMolecules)
                + " Err:\t" + (errJEMUExVMI / numberMolecules));

        double sumJEMUEyVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(14);
        double errJEMUEyVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(14);
        System.out.println("JEMUEyVMI:\t" + " Value:\t" + (sumJEMUEyVMI / numberMolecules)
                + " Err:\t" + (errJEMUEyVMI / numberMolecules));


        double sumJEMUExVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(17);
        double errJEMUExVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(17);
        System.out.println("JEMUExVSquareMI:\t" + " Value:\t" + (sumJEMUExVSquareMI / numberMolecules)
                + " Err:\t" + (errJEMUExVSquareMI / numberMolecules));

        double sumJEMUEyVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(18);
        double errJEMUEyVSquareMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(18);
        System.out.println("JEMUEyVSquareMI:\t" + " Value:\t" + (sumJEMUEyVSquareMI / numberMolecules)
                + " Err:\t" + (errJEMUEyVSquareMI / numberMolecules));


        double sumJEEMJEJEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(15);
        double errJEEMJEJEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(15);
        System.out.println("JEEMJEJEVIMI:\t" + " Value:\t" + (sumJEEMJEJEVMI / numberMolecules)
                + " Err:\t" + (errJEEMJEJEVMI / numberMolecules));


        double sumUEEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(16);
        double errUEEVMI = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(16);
        System.out.println("UEEVMI:\t\t" + " Value:\t" + (sumUEEVMI / numberMolecules)
                + " Err:\t" + (errUEEVMI / numberMolecules));

        endTime = System.currentTimeMillis();
        System.out.println("Total_Time: " + (endTime - startTime) / (1000.0 * 60.0));
    }

    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
        public boolean mSquare = true;
        public boolean aEE = true;
        public double temperature = 1.7;// Kelvin
        public double interactionS = 1.3;
        public double dipoleMagnitude = 1.6;
        public int nCells = 3;//number of atoms is nCells*nCells
        public int steps = 10000;
    }
}
