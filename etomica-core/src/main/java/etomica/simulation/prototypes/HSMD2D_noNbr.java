/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.*;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */

public class HSMD2D_noNbr extends Simulation {

    private static final long serialVersionUID = 1L;
    public ActivityIntegrate activityIntegrate;
    public AccumulatorAverage pressureAverage;
    public AccumulatorHistory pressureHistory;
    public AccumulatorAverageCollapsing temperatureAverage;
    public AccumulatorHistory temperatureHistory;
    public Box box;
    public SpeciesSpheresMono species;
    public IntegratorHard integrator;

    public HSMD2D_noNbr(Space _space) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box = new Box(new BoundaryRectangularNonperiodic(space), space);
        addBox(box);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{10, 10}));
        box.setNMolecules(species, 64);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        P2HardSphere potential = new P2HardSphere(space);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});
        P1HardBoundary potentialBoundary = new P1HardBoundary(space);
        potentialMaster.addPotential(potentialBoundary, new AtomType[]{species.getLeafType()});
//        potentialBoundary.setActive(0,true,true);
//        potentialBoundary.setActive(1,true,true);
//        potentialBoundary.setActive(0,false,true);
//        potentialBoundary.setActive(1,false,true);

        integrator.setBox(box);
//        integrator.setIsothermal(true);

//        MeterPressureHard meterPressure = new MeterPressureHard(integrator);
//        meterPressure.setBox(box);
//        pressureAverage = new AccumulatorAverage();
//        DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
//        IntervalActionAdapter pressureAction = new IntervalActionAdapter(pressurePump, integrator);
//
//        pressureHistory = new AccumulatorHistory();
//        pressureAverage.makeDataPusher(
//          new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE}).
//                                      addDataSink(pressureHistory);

        MeterTemperature meterTemperature = new MeterTemperature(box, space.D());
        temperatureAverage = new AccumulatorAverageCollapsing();
        DataPump temperaturePump = new DataPump(meterTemperature, temperatureAverage);
        integrator.getEventManager().addListener(new IntegratorListenerAction(temperaturePump));

//        pressureHistory = new AccumulatorHistory();
//        pressureAverage.makeDataPusher(
//          new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE}).
//                                      addDataSink(pressureHistory);
        temperatureHistory = new AccumulatorHistory();
        temperatureAverage.addDataSink(temperatureHistory, new AccumulatorAverage.StatType[]{AccumulatorAverage.AVERAGE});
        DataSourceCountTime timeCounter = new DataSourceCountTime(integrator);
        temperatureHistory.setTimeDataSource(timeCounter);
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final String APP_NAME = "HSMD2D no Nbr";

        Space sp = Space2D.getInstance();
        final HSMD2D_noNbr sim = new HSMD2D_noNbr(sp);
        final SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, sim.space, sim.getController());
        sim.activityIntegrate.setSleepPeriod(10);
//        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
//        pressureDisplay.setAccumulator(sim.pressureAverage);
//        DisplayPlot pressurePlot = new DisplayPlot();
//        sim.pressureHistory.addDataSink(pressurePlot.makeDataSink());
        DisplayTextBoxesCAE temperatureDisplay = new DisplayTextBoxesCAE();
        temperatureDisplay.setAccumulator(sim.temperatureAverage);
        DisplayPlot temperaturePlot = new DisplayPlot();
        temperaturePlot.setLabel("Temp");
        sim.temperatureHistory.setDataSink(temperaturePlot.getDataSet().makeDataSink());
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim, sp, sim.getController()));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);
        IAction repaintAction = graphic.getPaintAction(sim.box);

        nSelector.setPostAction(repaintAction);
        graphic.getController().getReinitButton().setPostAction(repaintAction);

        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(repaintAction));

        DeviceThermoSlider thermo = new DeviceThermoSlider(sim.getController(), sim.integrator);
        thermo.setMinimum(0.0);
        thermo.setMaximum(600.0);
        thermo.setSliderMajorValues(3);
        thermo.setAdiabatic();

        graphic.add(nSelector);
        graphic.add(thermo);
//        graphic.add(pressureDisplay);
//        graphic.add(pressurePlot);
        graphic.add(temperatureDisplay);
        graphic.add(temperaturePlot);
        graphic.makeAndDisplayFrame(APP_NAME);
    }//end of main

}
