/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.IListener;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class LJPIMC extends Simulation {

    public final PotentialMasterCellPI potentialMaster;
    public final IntegratorMC integrator;
    public final PotentialMasterBonding pmBonding;
    public final MCMoveHOReal2 ringMove;
    public final Box box;
    public final MCMoveMolecule mcMoveTranslate;
    public final PotentialComputeAggregate pmAgg;

    /**
     * Creates simulation with the given parameters
     */
    public LJPIMC(Space space, double mass, int numAtoms, int nBeads, double temperature, double density, double rc, double omega2) {
        super(Space3D.getInstance());

        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(AtomType.simple("A", mass/nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        potentialMaster = new PotentialMasterCellPI(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = false;

        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);

        double beta = 1/temperature;
        double hbar = 1;
        double omegaN = nBeads/(hbar*beta);

        double k2_kin = nBeads == 1 ? 0 : (mass*omegaN*omegaN);

        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,(i+1)%nBeads};
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        P2LennardJones p2lj = new P2LennardJones(1, 1.0/nBeads);
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);

        pmAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);

        integrator = new IntegratorMC(pmAgg, random, temperature, box);

        mcMoveTranslate = new MCMoveMolecule(random, potentialMaster, box);
        if (omega2 == 0) integrator.getMoveManager().addMCMove(mcMoveTranslate);

        ((MCMoveStepTracker)mcMoveTranslate.getTracker()).setNoisyAdjustment(true);

        ringMove = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box);
        if (false) {
            integrator.getMoveManager().addMCMove(ringMove);
        }
        else {
            integrator.getMoveManager().addMCMove(mcMoveTranslate);
            integrator.getMoveManager().addMCMove(new MCMoveMoleculeRotate(random, potentialMaster, box));
            integrator.getMoveManager().addMCMove(new MCMoveAtom(random, pmAgg, box));
        }
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
//            params.steps = 1000000;
            params.nBeads = 16;
            params.isGraphic = true;
        }

        Space space = Space.getInstance(params.D);
        double mass = params.mass;
        int numAtoms = params.numAtoms;
        int nBeads = params.nBeads;
        double temperature = params.temperature;
        double density = params.density;
        double rc = params.rc;
        double omega2 = params.k2/mass;
        boolean isGraphic = params.isGraphic;

        LJPIMC sim = new LJPIMC(space, mass, numAtoms, nBeads, temperature, density, rc, omega2);
        long steps = params.steps;
        int interval = numAtoms;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("mass: " + mass);
        System.out.println("k2: " + params.k2);
        System.out.println("N: " + numAtoms);
        System.out.println("nBeads: " + nBeads);
        System.out.println("T: " + temperature);
        System.out.println("density: " + density);
        System.out.println("steps: " + steps);
        System.out.println("rc: " + rc);

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;
                public Color getAtomColor(IAtom a) {
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
                    return allColors[(768 * a.getIndex() / (nBeads))];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            MeterMSDHO meterMSD = new MeterMSDHO(sim.box);
            AccumulatorHistory historyMSD = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyMSD.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
            DataPumpListener pumpMSD = new DataPumpListener(meterMSD, historyMSD, interval);
            sim.integrator.getEventManager().addListener(pumpMSD);
            DisplayPlotXChart plotMSD = new DisplayPlotXChart();
            plotMSD.setLabel("MSD");
            plotMSD.setDoLegend(false);
            historyMSD.addDataSink(plotMSD.makeSink("MSD history"));
            simGraphic.add(plotMSD);

            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
            DataPumpListener pumpPE = new DataPumpListener(meterPE, historyPE, interval);
            sim.integrator.getEventManager().addListener(pumpPE);
            MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
            AccumulatorHistory historyPE2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE2.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
            DataPumpListener pumpPE2 = new DataPumpListener(meterPE2, historyPE2, interval);
            sim.integrator.getEventManager().addListener(pumpPE2);
            DisplayPlotXChart plotPE = new DisplayPlotXChart();
            plotPE.setLabel("PE");
//            plotPE.setDoLegend(false);
            historyPE.addDataSink(plotPE.makeSink("PE history"));
            plotPE.setLegend(new DataTag[]{meterPE.getTag()}, "Integrator PE");
            historyPE2.addDataSink(plotPE.makeSink("PE2 history"));
            plotPE.setLegend(new DataTag[]{meterPE2.getTag()}, "recompute PE");
            simGraphic.add(plotPE);

            Vector[] latticePositions = space.makeVectorArray(numAtoms);
            Vector COM0 = space.makeVector();
            for (IAtom atom : sim.box.getLeafList()) {
                if (atom.getIndex() == 0) {
                    latticePositions[atom.getParentGroup().getIndex()].E(atom.getPosition());
                }
                COM0.PE(atom.getPosition());
            }
            COM0.TE(1.0/sim.box.getLeafList().size());

            DataSourceScalar meterCOM = new DataSourceScalar("COM", Length.DIMENSION) {
                @Override
                public double getDataAsScalar() {
                    Vector COM = sim.box.getSpace().makeVector();
                    Vector dr = sim.box.getSpace().makeVector();
                    for (IAtom atom : sim.box.getLeafList()) {
                        Vector R = latticePositions[atom.getParentGroup().getIndex()];
                        dr.Ev1Mv2(atom.getPosition(), R);
                        sim.box.getBoundary().nearestImage(dr);
                        dr.PE(R);
                        COM.PE(dr);
                    }
                    COM.TE(1.0/sim.box.getLeafList().size());
                    return Math.sqrt(COM.squared());
                }
            };
            AccumulatorHistory historyCOM = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyCOM.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
            DataPumpListener pumpCOM = new DataPumpListener(meterCOM, historyCOM, interval);
            sim.integrator.getEventManager().addListener(pumpCOM);
            DisplayPlotXChart plotCOM = new DisplayPlotXChart();
            plotCOM.setLabel("COM");
            historyCOM.addDataSink(plotCOM.makeSink("COM"));
            plotCOM.setLegend(new DataTag[]{meterCOM.getTag()}, "COM");
            simGraphic.add(plotCOM);

            MeterAcceptance meterAcceptance = new MeterAcceptance(sim.ringMove);
            AccumulatorHistory historyAcceptance = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyAcceptance.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
            DataPumpListener pumpAcceptance = new DataPumpListener(meterAcceptance, historyAcceptance, interval);
            sim.integrator.getEventManager().addListener(pumpAcceptance);
            historyAcceptance.addDataSink(plotCOM.makeSink("acceptance"));
            plotCOM.setLegend(new DataTag[]{meterAcceptance.getTag()}, "acceptance");
            sim.integrator.getMoveEventManager().addListener(meterAcceptance);


            simGraphic.makeAndDisplayFrame(" PIMC ");

            return;
        }

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        System.out.println("equilibration finished");

        // <U>
//        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener pump = new DataPumpListener(meterPE, acc, interval);
//        sim.integrator.getEventManager().addListener(pump);

//            public MeterPIPrim(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, int nBeads, double betaN, Box box) {

        MeterPIPrim meterPrim = new MeterPIPrim(sim.pmBonding, sim.potentialMaster, nBeads, sim.ringMove.betaN, sim.box);
        AccumulatorAverageCovariance accumulatorPrim = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpPrim);

        MeterPIVir meterVir = new MeterPIVir(sim.potentialMaster, sim.ringMove.betaN, nBeads, sim.box);
        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpVir);

        MeterPICentVir meterCentVir = new MeterPICentVir(sim.potentialMaster, 1/temperature, nBeads, sim.box);
        AccumulatorAverageCovariance accumulatorCentVir = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpCentVir = new DataPumpListener(meterCentVir, accumulatorCentVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpCentVir);

        MeterPIHMAc meterHMAc = new MeterPIHMAc(sim.potentialMaster, sim.ringMove.betaN, nBeads, sim.box);
        AccumulatorAverageCovariance accumulatorHMAc = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpHMAc = new DataPumpListener(meterHMAc, accumulatorHMAc, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMAc);

        MeterPIHMAReal2 meterHMAReal2 = new MeterPIHMAReal2(sim.pmBonding, sim.potentialMaster, sim.ringMove.beta, sim.ringMove);
        AccumulatorAverageCovariance accumulatorHMAReal2 = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpHMAReal2 = new DataPumpListener(meterHMAReal2, accumulatorHMAReal2, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMAReal2);

        //Run ...
        sim.integrator.resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        long t2 = System.currentTimeMillis();

        System.out.println();
        System.out.println("acceptance probability: " + sim.ringMove.getTracker().acceptanceProbability());
        System.out.println("translate step size: " + sim.mcMoveTranslate.getStepSize());

        //<U>
//        DataGroup dataPE = (DataGroup) acc.getData();
//        double avg = dataPE.getValue(acc.AVERAGE.index) / numAtoms;
//        double err = dataPE.getValue(acc.ERROR.index) / numAtoms;
//        double cor = dataPE.getValue(acc.BLOCK_CORRELATION.index);
//
//        System.out.println();
//        double En_sim = avg + sim.box.getSpace().D()/2.0*temperature;
//        System.out.println("energy avg: " + En_sim + "  err: " + err + "  cor: " + cor);

        //Prim
        DataGroup dataPrim = (DataGroup) accumulatorPrim.getData();
        double avgEnPrim = dataPrim.getValue(accumulatorPrim.AVERAGE.index)/numAtoms;
        double errEnPrim = dataPrim.getValue(accumulatorPrim.ERROR.index)/numAtoms;
        double corEnPrim = dataPrim.getValue(accumulatorPrim.BLOCK_CORRELATION.index);
        System.out.println(" En_prim: " + avgEnPrim + " +/- " + errEnPrim + " cor: " + corEnPrim);

        //Vir
        DataGroup dataVir = (DataGroup) accumulatorVir.getData();
        double avgEnVir = dataVir.getValue(accumulatorVir.AVERAGE.index)/numAtoms;
        double errEnVir = dataVir.getValue(accumulatorVir.ERROR.index)/numAtoms;
        double corEnVir = dataVir.getValue(accumulatorVir.BLOCK_CORRELATION.index);
        System.out.println(" En_vir: " + avgEnVir + " +/- " + errEnVir + " cor: " + corEnVir);

        //CentVir
        DataGroup dataCentVir = (DataGroup) accumulatorCentVir.getData();
        double avgEnCentVir = dataCentVir.getValue(accumulatorCentVir.AVERAGE.index)/numAtoms;
        double errEnCentVir = dataCentVir.getValue(accumulatorCentVir.ERROR.index)/numAtoms;
        double corEnCentVir = dataCentVir.getValue(accumulatorCentVir.BLOCK_CORRELATION.index);
        System.out.println(" En_cvir: " + avgEnCentVir + " +/- " + errEnCentVir + " cor: " + corEnCentVir);

        //HMAc
        DataGroup dataHMAc = (DataGroup) accumulatorHMAc.getData();
        double avgEnHMAc = dataHMAc.getValue(accumulatorHMAc.AVERAGE.index)/numAtoms;
        double errEnHMAc = dataHMAc.getValue(accumulatorHMAc.ERROR.index)/numAtoms;
        double corEnHMAc = dataHMAc.getValue(accumulatorHMAc.BLOCK_CORRELATION.index);
        System.out.println(" En_hmac: " + avgEnHMAc + " +/- " + errEnHMAc + " cor: " + corEnHMAc);

        //HMA-EC-staging
        DataGroup dataHMAReal2 = (DataGroup) accumulatorHMAReal2.getData();
        double avgEnHMAReal2 = dataHMAReal2.getValue(accumulatorHMAReal2.AVERAGE.index)/numAtoms;
        double errEnHMAReal2 = dataHMAReal2.getValue(accumulatorHMAReal2.ERROR.index)/numAtoms;
        double corEnHMAReal2 = dataHMAReal2.getValue(accumulatorHMAReal2.BLOCK_CORRELATION.index);
        System.out.println(" En_hma2: " + avgEnHMAReal2 + " +/- " + errEnHMAReal2 + " cor: " + corEnHMAReal2);

        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public int D = 3;
        public int nBeads = 4;
        public double k2 = 219.231319;
        public long steps = 100000;
        public double density = 1.0;
        public double temperature = 0.5;
        public int numAtoms = 108;
        public double mass = 100.0;
        public double rc = 2.5;
        public boolean isGraphic = false;
    }

    public static class MeterAcceptance extends DataSourceScalar implements IListener<MCMoveEvent> {

        protected final MCMove move;
        protected double chiSum = 0;
        protected int numTrials = 0;

        public MeterAcceptance(MCMove move) {
            super("acceptance", Null.DIMENSION);
            this.move = move;
        }

        @Override
        public double getDataAsScalar() {
            double avg = chiSum / numTrials;
            chiSum = 0;
            numTrials = 0;
            return avg;
        }

        @Override
        public void actionPerformed(MCMoveEvent event) {
            if (!(event instanceof MCMoveTrialCompletedEvent) || (move != null && event.getMCMove() != move)) return;
            chiSum += Math.min(((MCMoveTrialCompletedEvent)event).chi, 1);
            numTrials++;
        }
    }
}