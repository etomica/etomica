/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.nucleation;

import etomica.action.IAction;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.modules.vle.DeviceThermoSliderGEMC;
import etomica.space.Space;
import etomica.units.Pixel;

import java.util.List;

public class SWVLE extends SimulationGraphic {

    private final static String APP_NAME = "SW / VLE";
    private final static int REPAINT_INTERVAL = 200;
    public boolean showNumMoleculesPlots = true;

    public SWVLE(final SWVLESim sim, Space _space) {
        super(sim, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        getDisplayBox(sim.boxLiquid).setPixelUnit(new Pixel(8));
        getPanel().tabbedPane.setTitleAt(0, "Liquid");
        getDisplayBox(sim.boxVapor).setPixelUnit(new Pixel(8));
        getPanel().tabbedPane.setTitleAt(1, "Vapor");

        final IAction resetMCMoves = new IAction() {
            public void actionPerformed() {
                List<MCMove> moves = sim.integratorLiquid.getMoveManager().getMCMoves();
                for (int i = 0; i < moves.size(); i++) {
                    if (moves.get(i).getTracker() instanceof MCMoveStepTracker) {
                        ((MCMoveStepTracker) moves.get(i).getTracker()).resetAdjustStep();
                    }
                }

                moves = sim.integratorVapor.getMoveManager().getMCMoves();
                for (int i = 0; i < moves.size(); i++) {
                    if (moves.get(i).getTracker() instanceof MCMoveStepTracker) {
                        ((MCMoveStepTracker) moves.get(i).getTracker()).resetAdjustStep();
                    }
                }

                moves = sim.integratorGEMC.getMoveManager().getMCMoves();
                for (int i = 0; i < moves.size(); i++) {
                    if (moves.get(i).getTracker() instanceof MCMoveStepTracker) {
                        ((MCMoveStepTracker) moves.get(i).getTracker()).resetAdjustStep();
                    }
                }
            }
        };

        DeviceThermoSliderGEMC thermoSlider = new DeviceThermoSliderGEMC(sim.getController());
        thermoSlider.setPrecision(2);
        thermoSlider.setIntegrators(sim.integratorLiquid, sim.integratorVapor, sim.integratorGEMC);
        thermoSlider.setIsothermal();
        thermoSlider.setMinimum(0);
        thermoSlider.setMaximum(2);
        thermoSlider.doUpdate();
        thermoSlider.setIsothermalButtonsVisibility(false);
        thermoSlider.setPostAction(resetMCMoves);
        add(thermoSlider);

        MeterDensity meterDensityLiquid = new MeterDensity(sim.getSpace());
        meterDensityLiquid.setBox(sim.boxLiquid);
        DataFork fork = new DataFork();
        DataPump pumpLiquidDensity = new DataPump(meterDensityLiquid, fork);
        getController().getDataStreamPumps().add(pumpLiquidDensity);
        AccumulatorAverageCollapsing avgLiquidDensity = new AccumulatorAverageCollapsing();
        avgLiquidDensity.setPushInterval(5);
        fork.addDataSink(avgLiquidDensity);
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integratorGEMC);
        AccumulatorHistory historyDensityLiquid = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyDensityLiquid.setTimeDataSource(stepCounter);
        fork.addDataSink(historyDensityLiquid);
        IntegratorListenerAction pumpLiquidDensityListener = new IntegratorListenerAction(pumpLiquidDensity);
        sim.integratorLiquid.getEventManager().addListener(pumpLiquidDensityListener);
        pumpLiquidDensityListener.setInterval(100);

        MeterDensity meterDensityVapor = new MeterDensity(sim.getSpace());
        meterDensityVapor.setBox(sim.boxVapor);
        fork = new DataFork();
        DataPump pumpVaporDensity = new DataPump(meterDensityVapor, fork);
        getController().getDataStreamPumps().add(pumpVaporDensity);
        AccumulatorAverageCollapsing avgVaporDensity = new AccumulatorAverageCollapsing();
        avgVaporDensity.setPushInterval(5);
        fork.addDataSink(avgVaporDensity);
        AccumulatorHistory historyDensityVapor = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyDensityVapor.setTimeDataSource(stepCounter);
        fork.addDataSink(historyDensityVapor);
        IntegratorListenerAction pumpVaporDensityListener = new IntegratorListenerAction(pumpVaporDensity);
        sim.integratorVapor.getEventManager().addListener(pumpVaporDensityListener);
        pumpVaporDensityListener.setInterval(100);

        DisplayTextBoxesCAE liquidDensityDisplay = new DisplayTextBoxesCAE();
        liquidDensityDisplay.setAccumulator(avgLiquidDensity);
        add(liquidDensityDisplay);
        DisplayPlotXChart liquidDensityHistoryPlot = new DisplayPlotXChart();
        liquidDensityHistoryPlot.setLabel("Liquid Density");
        add(liquidDensityHistoryPlot);
        historyDensityLiquid.addDataSink(liquidDensityHistoryPlot.getDataSet().makeDataSink());
        liquidDensityHistoryPlot.setLegend(new DataTag[]{historyDensityLiquid.getTag()}, "Density (mol/L)");
        historyDensityLiquid.setPushInterval(1);

        DisplayTextBoxesCAE vaporDensityDisplay = new DisplayTextBoxesCAE();
        vaporDensityDisplay.setAccumulator(avgVaporDensity);
        add(vaporDensityDisplay);
        DisplayPlotXChart vaporDensityHistoryPlot = new DisplayPlotXChart();
        vaporDensityHistoryPlot.setLegend(new DataTag[]{historyDensityVapor.getTag()}, "Density (mol/L)");
        vaporDensityHistoryPlot.setLabel("Vapor Density");
        add(vaporDensityHistoryPlot);
        historyDensityVapor.addDataSink(vaporDensityHistoryPlot.getDataSet().makeDataSink());
        historyDensityVapor.setPushInterval(1);

        if (showNumMoleculesPlots) {
            MeterNMolecules meterNMoleculesLiquid = new MeterNMolecules();
            meterNMoleculesLiquid.setBox(sim.boxLiquid);
            AccumulatorHistory historyNMoleculesLiquid = new AccumulatorHistory(new HistoryCollapsingAverage(100));
            historyNMoleculesLiquid.setTimeDataSource(stepCounter);
            DataPump pumpNMoleculesLiquid = new DataPump(meterNMoleculesLiquid, historyNMoleculesLiquid);
            IntegratorListenerAction pumpNMoleculesLiquidListener = new IntegratorListenerAction(pumpNMoleculesLiquid);
            sim.integratorLiquid.getEventManager().addListener(pumpNMoleculesLiquidListener);
            pumpNMoleculesLiquidListener.setInterval(100);
            getController().getDataStreamPumps().add(pumpNMoleculesLiquid);
            MeterNMolecules meterNMoleculesVapor = new MeterNMolecules();
            meterNMoleculesVapor.setBox(sim.boxVapor);
            AccumulatorHistory historyNMoleculesVapor = new AccumulatorHistory(new HistoryCollapsingAverage(100));
            historyNMoleculesVapor.setTimeDataSource(stepCounter);
            DataPump pumpNMoleculesVapor = new DataPump(meterNMoleculesVapor, historyNMoleculesVapor);
            IntegratorListenerAction pumpNMoleculesVaporListener = new IntegratorListenerAction(pumpNMoleculesVapor);
            sim.integratorVapor.getEventManager().addListener(pumpNMoleculesVaporListener);
            pumpNMoleculesVaporListener.setInterval(100);
            getController().getDataStreamPumps().add(pumpNMoleculesVapor);

            DisplayPlotXChart nMoleculesLiquidHistoryPlot = new DisplayPlotXChart();
            nMoleculesLiquidHistoryPlot.setLabel("# of Liquid Atoms");
            add(nMoleculesLiquidHistoryPlot);
            historyNMoleculesLiquid.addDataSink(nMoleculesLiquidHistoryPlot.getDataSet().makeDataSink());
            DisplayPlotXChart nMoleculesVaporHistoryPlot = new DisplayPlotXChart();
            nMoleculesVaporHistoryPlot.setLabel("# of Vapor Atoms");
            add(nMoleculesVaporHistoryPlot);
            historyNMoleculesVapor.addDataSink(nMoleculesVaporHistoryPlot.getDataSet().makeDataSink());
        }

    }

    public static void main(String[] args) {
        int D = 2;
        if (args.length > 0) D = Integer.parseInt(args[0]);
        SWVLESim sim = new SWVLESim(D);
        SWVLE SWVLE = new SWVLE(sim, sim.getSpace());
        SWVLE.makeAndDisplayFrame();
    }
}
