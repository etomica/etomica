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
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
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
    public LJPIMC(Space space, double mass, int numAtoms, int nBeads, double temperature, double density, double rc) {
        super(Space3D.getInstance());

        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(AtomType.simple("A", mass/nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        potentialMaster = new PotentialMasterCellPI(getSpeciesManager(), box, 2, BondingInfo.noBonding());

        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);

        double beta = 1/temperature;
        double hbar = 1;
        double omegaN = nBeads/(hbar*beta);

        double k2_kin = nBeads == 1 ? 0 : (mass*omegaN*omegaN/nBeads);

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
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);

        pmAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);

        integrator = new IntegratorMC(pmAgg, random, temperature, box);

        double omega2 = .5;
        mcMoveTranslate = new MCMoveMolecule(random, potentialMaster, box);
        if (omega2 == 0) integrator.getMoveManager().addMCMove(mcMoveTranslate);

        ((MCMoveStepTracker)mcMoveTranslate.getTracker()).setNoisyAdjustment(true);

        System.out.println("omega2: "+omega2);
        ringMove = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box);
        integrator.getMoveManager().addMCMove(ringMove);
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
            params.steps = 100000;
        }

        Space space = Space.getInstance(params.D);
        double mass = params.mass;
        int numAtoms = params.numAtoms;
        int nBeads = params.nBeads;
        double temperature = params.temperature;
        double density = params.density;
        double rc = params.rc;

        LJPIMC sim = new LJPIMC(space, mass, numAtoms, nBeads, temperature, density, rc);
        long steps = params.steps;
        int interval = numAtoms;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperature);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        if (false) {
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

            simGraphic.makeAndDisplayFrame(" PIMC ");


            return;
        }

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        System.out.println("equilibration finished");

        // data collection
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(meterPE, acc, interval);
        sim.integrator.getEventManager().addListener(pump);

        // data collection
        MeterMSDHO meterMSD = new MeterMSDHO(sim.box);
        AccumulatorAverageFixed accMSD = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpMSD = new DataPumpListener(meterMSD, accMSD, interval);
        sim.integrator.getEventManager().addListener(pumpMSD);

        sim.integrator.resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        long t2 = System.currentTimeMillis();

        System.out.println();
        System.out.println("acceptance probability: " + sim.ringMove.getTracker().acceptanceProbability());
        System.out.println("translate step size: " + sim.mcMoveTranslate.getStepSize());

        DataGroup dataPE = (DataGroup) acc.getData();
        double avg = dataPE.getValue(acc.AVERAGE.index) / numAtoms;
        double err = dataPE.getValue(acc.ERROR.index) / numAtoms;
        double cor = dataPE.getValue(acc.BLOCK_CORRELATION.index);

        DataGroup dataMSD = (DataGroup) accMSD.getData();
        double avgMSD = dataMSD.getValue(accMSD.AVERAGE.index);
        double errMSD = dataMSD.getValue(accMSD.ERROR.index);
        double corMSD = dataMSD.getValue(accMSD.BLOCK_CORRELATION.index);

        System.out.println();
        System.out.println("energy avg: " + avg + "  err: " + err + "  cor: " + cor);
        System.out.println("MSD avg: " + avgMSD + "  err: " + errMSD + "  cor: " + corMSD);
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public int D = 3;
        public int nBeads = 16;
        public long steps = 1000000;
        public double density = 1;
        public double temperature = 0.5;
        public int numAtoms = 108;
        public double mass = 1000;
        public double rc = 2.5;
    }
}
