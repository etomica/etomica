/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.reproducibility;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class MethaneNPTMC extends Simulation {

    public final PotentialMasterCell potentialMaster;
    public final IntegratorMC integrator;
    public MCMoveVolume mcMoveVolume;
    public MCMoveAtom mcMoveAtom;

    /**
     * Creates simulation with the given parameters
     */
    public MethaneNPTMC(int numAtoms, double density, double temperature, double rc, double pressure, boolean shift, boolean lrc) {
        super(Space3D.getInstance());

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        Box box = new Box(space);
        addBox(box);
        potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = lrc;

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        double epsilon = Kelvin.UNIT.toSim(148);
        double sigma = 3.73;
        P2LennardJones p2lj = new P2LennardJones(sigma, epsilon);
        IPotential2 p2 = shift ? new P2SoftSphericalTruncatedShifted(p2lj, rc) : new P2SoftSphericalTruncated(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);

        integrator = new IntegratorMC(potentialMaster, random, temperature, box);

        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
//        integrator.getMoveManager().setFrequency(mcMoveAtom, 1);

        if (pressure >= 0) {
            mcMoveVolume = new MCMoveVolume(integrator, random, pressure);
            integrator.getMoveManager().addMCMove(mcMoveVolume);
        }

        getController().addActivity(new ActivityIntegrate(integrator));
    }

    public static void main(String[] args) {

        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
            params.numAtoms = 900;
            params.steps = 120 * 1000 * params.numAtoms;
        }

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressureKPa = params.pressureKPa;
        double pressure = pUnit.toSim(pressureKPa);

        MethaneNPTMC sim = new MethaneNPTMC(params.numAtoms, params.density, temperature, params.rc, pressure, params.shift, params.lrc);
        // TraPPE says CH4 is 16.04
        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(16.04/Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

        if (true) {
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            ((DiameterHashByType)graphic.getDisplayBox(sim.box()).getDiameterHash()).setDiameter(sim.species().getLeafType(), 3.73);

            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator);
            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPE.setTimeDataSource(timeSource);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 10);
            sim.integrator.getEventManager().addListener(pumpPE);

            DisplayPlot historyPE = new DisplayPlot();
            accPE.setDataSink(historyPE.getDataSet().makeDataSink());
            historyPE.setLabel("PE");
            historyPE.setUnit(new SimpleUnit(Null.DIMENSION, params.numAtoms, "1/N", "1/N", false));
            graphic.add(historyPE);

            MeterDensity meterDensity = new MeterDensity(sim.box());
            AccumulatorHistory accDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
            accDensity.setTimeDataSource(timeSource);
            DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, 10);
            sim.integrator.getEventManager().addListener(pumpDensity);

            DisplayPlot historyDensity = new DisplayPlot();
            accDensity.setDataSink(historyDensity.getDataSet().makeDataSink());
            historyDensity.setLabel("Density");
            historyDensity.setUnit(dUnit);
            graphic.add(historyDensity);

            MeterPressure meterPressure = new MeterPressure(sim.box(), sim.potentialMaster);
            meterPressure.setTemperature(temperature);
            AccumulatorHistory accPressure = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPressure.setTimeDataSource(timeSource);
            DataPumpListener pumpPressure = new DataPumpListener(meterPressure, accPressure, 2*params.numAtoms);
            sim.integrator.getEventManager().addListener(pumpPressure);

            DisplayPlotXChart historyPressure = new DisplayPlotXChart();
            accPressure.setDataSink(historyPressure.getDataSet().makeDataSink());
            historyPressure.setLabel("Pressure");
            historyPressure.setUnit(pUnit);
            graphic.add(historyPressure);

            graphic.makeAndDisplayFrame();
            return;
        }

        long steps = params.steps;
        int interval = 10;
        int blocks = 100;
        long blockSize = Math.max(steps / (interval * blocks), 1);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + temperatureK);
        System.out.println("P: " + pressureKPa);
        System.out.println("initial density: " + dUnit.fromSim(params.density));
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        System.out.println("equilibration finished\n");
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.mcMoveVolume.getTracker().reset();
        sim.mcMoveAtom.getTracker().reset();

        // data collection
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed accPE = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, interval);
        sim.integrator.getEventManager().addListener(pumpPE);

        MeterDensity meterDensity = new MeterDensity(sim.box());
        AccumulatorAverageFixed accDensity = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, interval);
        sim.integrator.getEventManager().addListener(pumpDensity);

        sim.integrator.resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        long t2 = System.currentTimeMillis();

        MCMoveStepTracker volumeTracker = (MCMoveStepTracker) sim.mcMoveVolume.getTracker();
        System.out.println("volume change fraction: "+((double)volumeTracker.nTrials)/params.steps);
        System.out.println("volume change step size: "+sim.mcMoveVolume.getStepSize());
        System.out.println("volume change acceptance: "+volumeTracker.acceptanceProbability());

        MCMoveStepTracker displacementTracker = (MCMoveStepTracker) sim.mcMoveAtom.getTracker();
        System.out.println("displacement fraction: "+((double)displacementTracker.nTrials)/params.steps);
        System.out.println("displacement step size: "+sim.mcMoveAtom.getStepSize());
        System.out.println("displacement acceptance: "+displacementTracker.acceptanceProbability());
        System.out.println();

        DataGroup dataPE = (DataGroup) accPE.getData();
        int numAtoms = sim.getBox(0).getLeafList().size();
        double avg = dataPE.getValue(accPE.AVERAGE.index) / numAtoms;
        double err = dataPE.getValue(accPE.ERROR.index) / numAtoms;
        double cor = dataPE.getValue(accPE.BLOCK_CORRELATION.index);

        DataGroup dataDensity = (DataGroup) accDensity.getData();
        double avgDensity = dataDensity.getValue(accPE.AVERAGE.index);
        double errDensity = dataDensity.getValue(accPE.ERROR.index);
        double corDensity = dataDensity.getValue(accPE.BLOCK_CORRELATION.index);

        System.out.println("energy avg: " + avg + "  err: " + err + "  cor: " + cor);
        System.out.println("density avg: " + avgDensity + "  err: " + errDensity + "  cor: " + corDensity);
        System.out.println("density avg (g/cm^3): " + dUnit.fromSim(avgDensity) + "  err: " + dUnit.fromSim(errDensity) + "  cor: " + corDensity);
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class SimParams extends ParameterBase {
        public long steps = 100000000;
        public double density = 0.01;
        public double temperatureK = 140;
        public int numAtoms = 1000;
        public double rc = 14;
        public double pressureKPa = 1318;
        public boolean shift = false;
        public boolean lrc = true;
    }
}
