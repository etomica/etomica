/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.IAction;
import etomica.action.SimulationDataAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPressureHardFasterer;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.*;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */

public class HSMD2D_noNbrFasterer extends Simulation {

    public AccumulatorAverageCollapsing pressureAverage;
    public AccumulatorHistory pressureHistory;
    public AccumulatorAverageCollapsing temperatureAverage;
    public AccumulatorHistory temperatureHistory;
    public Box box;
    public SpeciesGeneral species;
    public IntegratorHardFasterer integrator;
    public DataPump pressurePump;
    public DataPump temperaturePump;
    public MeterPressureHardFasterer meterPressure;

    public HSMD2D_noNbrFasterer() {
        super(Space2D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));

        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);
        PotentialComputePair potentialMaster = new PotentialComputePair(this, box, neighborManager);
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);

        box.getBoundary().setBoxSize(Vector.of(new double[]{10, 10}));
        this.getController().addActivity(new ActivityIntegrate(integrator));
        box.setNMolecules(species, 64);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        P2HardGeneric potential = P2HardSphere.makePotential(1.0);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), potential);
        P1HardBoundary potentialBoundary = new P1HardBoundary(space, false, box);
        pcField.setFieldPotential(species.getLeafType(), potentialBoundary);
//        potentialBoundary.setActive(0,true,true);
//        potentialBoundary.setActive(1,true,true);
//        potentialBoundary.setActive(0,false,true);
//        potentialBoundary.setActive(1,false,true);
        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), IntegratorHardFasterer.extractFieldPotentials(pcField), neighborManager, random, 0.05, 1, box, getSpeciesManager(), null);
        integrator.setIsothermal(true);

        meterPressure = new MeterPressureHardFasterer(integrator);
        pressureAverage = new AccumulatorAverageCollapsing();
        pressurePump = new DataPump(meterPressure, pressureAverage);
//        IntervalActionAdapter pressureAction = new IntervalActionAdapter(pressurePump, integrator);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pressurePump));

        pressureHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
        pressureAverage.addDataSink(pressureHistory, new AccumulatorAverage.StatType[]{pressureAverage.AVERAGE});

        MeterTemperature meterTemperature = new MeterTemperature(box, space.D());
        temperatureAverage = new AccumulatorAverageCollapsing();
        temperaturePump = new DataPump(meterTemperature, temperatureAverage);
        integrator.getEventManager().addListener(new IntegratorListenerAction(temperaturePump));


        temperatureHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
        temperatureAverage.addDataSink(temperatureHistory, new AccumulatorAverage.StatType[]{temperatureAverage.AVERAGE});
        DataSourceCountTimeFasterer timeCounter = new DataSourceCountTimeFasterer(integrator);

        temperatureHistory.setTimeDataSource(timeCounter);
        pressureHistory.setTimeDataSource(timeCounter);
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final String APP_NAME = "HSMD2D no Nbr";

        final HSMD2D_noNbrFasterer sim = new HSMD2D_noNbrFasterer();
        final SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        sim.getController().setSleepPeriod(10);
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
        pressureDisplay.setAccumulator(sim.pressureAverage);
        DisplayPlot pressurePlot = new DisplayPlot();
        pressurePlot.setLabel("Pressure");
        sim.pressureHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());

        DisplayTextBoxesCAE temperatureDisplay = new DisplayTextBoxesCAE();
        temperatureDisplay.setAccumulator(sim.temperatureAverage);
        DisplayPlot temperaturePlot = new DisplayPlot();
        temperaturePlot.setLabel("Temp");
        sim.temperatureHistory.setDataSink(temperaturePlot.getDataSet().makeDataSink());

        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        SimulationRestart simr = new SimulationRestart(sim);
        SimulationDataAction ra = simr.getDataResetAction();
        ra.getDataStreamPumps().add(sim.temperaturePump);
        ra.getDataStreamPumps().add(sim.pressurePump);

        IAction resetmeter = new IAction() {
            @Override
            public void actionPerformed() {
                sim.meterPressure.reset();
            }
        };

        simr.setPostAction(resetmeter);
        nSelector.setResetAction(simr);
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);
        IAction repaintAction = graphic.getPaintAction(sim.box);

        nSelector.setPostAction(repaintAction);
        graphic.getController().getReinitButton().setPostAction(repaintAction);
        graphic.getController().getReinitButton().setAction(simr);
        graphic.getController().getDataStreamPumps().add(sim.temperaturePump);
        graphic.getController().getDataStreamPumps().add(sim.pressurePump);

        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(repaintAction));

        DeviceThermoSlider thermo = new DeviceThermoSlider(sim.getController(), sim.integrator);
        thermo.setMinimum(0.0);
        thermo.setMaximum(600.0);
        thermo.setSliderMajorValues(3);
        thermo.setAdiabatic();

        graphic.add(nSelector);
        graphic.add(thermo);
        graphic.add(pressureDisplay);
        graphic.add(pressurePlot);
        graphic.add(temperatureDisplay);
        graphic.add(temperaturePlot);
        graphic.makeAndDisplayFrame(APP_NAME);
    }//end of main

}
