/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import java.util.List;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataTag;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressure;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.units.Angstrom;
import etomica.units.CompoundUnit;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pascal;
import etomica.units.Pixel;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.data.history.HistoryCollapsingAverage;

public class VLE extends SimulationGraphic {

    private final static String APP_NAME = "Virial / VLE";
    private final static int REPAINT_INTERVAL = 200;
    public boolean showNumMoleculesPlots = false;

    public VLE(final VLESim sim, Space _space) {
        super(sim, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, sim.getController());

        getDisplayBox(sim.boxLiquid).setPixelUnit(new Pixel(8));
        getPanel().tabbedPane.setTitleAt(0, "Liquid");
        getDisplayBox(sim.boxVapor).setPixelUnit(new Pixel(8));
        getPanel().tabbedPane.setTitleAt(1, "Vapor");

        final IAction resetMCMoves = new IAction() {
            public void actionPerformed() {
                List<MCMove> moves = sim.integratorLiquid.getMoveManager().getMCMoves();
                for (int i=0; i<moves.size(); i++) {
                    if (moves.get(i).getTracker() instanceof MCMoveStepTracker) {
                        ((MCMoveStepTracker)moves.get(i).getTracker()).resetAdjustStep();
                    }
                }

                moves = sim.integratorVapor.getMoveManager().getMCMoves();
                for (int i=0; i<moves.size(); i++) {
                    if (moves.get(i).getTracker() instanceof MCMoveStepTracker) {
                        ((MCMoveStepTracker)moves.get(i).getTracker()).resetAdjustStep();
                    }
                }
                
                moves = sim.integratorGEMC.getMoveManager().getMCMoves();
                for (int i=0; i<moves.size(); i++) {
                    if (moves.get(i).getTracker() instanceof MCMoveStepTracker) {
                        ((MCMoveStepTracker)moves.get(i).getTracker()).resetAdjustStep();
                    }
                }
            }
        };

        DeviceThermoSliderGEMC thermoSlider = new DeviceThermoSliderGEMC(sim.getController());
        thermoSlider.setUnit(Kelvin.UNIT);
        thermoSlider.setPrecision(1);
        thermoSlider.setIntegrators(sim.integratorLiquid, sim.integratorVapor);
        thermoSlider.setIsothermal();
        thermoSlider.setMinimum(200);
        thermoSlider.setMaximum(400);
        thermoSlider.doUpdate();
        thermoSlider.setIsothermalButtonsVisibility(false);
        thermoSlider.setPostAction(resetMCMoves);
        add(thermoSlider);
        
        final DeviceSlider sigmaSlider = new DeviceSlider(sim.getController(), new ModifierGeneral(sim, "sigma"));
        sigmaSlider.setPrecision(2);
        sigmaSlider.setMinimum(2);
        sigmaSlider.setMaximum(5);
        sigmaSlider.setShowValues(true);
        sigmaSlider.setEditValues(true);
        sigmaSlider.setShowBorder(true);
        sigmaSlider.setLabel("sigma (A)");
        add(sigmaSlider);
        sigmaSlider.doUpdate();
        final DiameterHashByType diameterManager = (DiameterHashByType)getDisplayBox(sim.boxLiquid).getDiameterHash();
        getDisplayBox(sim.boxVapor).setDiameterHash(diameterManager);
        diameterManager.setDiameter(sim.species.getLeafType(), sim.getSigma());
        sigmaSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                diameterManager.setDiameter(sim.species.getLeafType(), sigmaSlider.getValue());
                resetMCMoves.actionPerformed();
                getDisplayBox(sim.boxLiquid).repaint();
                getDisplayBox(sim.boxVapor).repaint();
            }
        });
        
        DeviceSlider epsilonSlider = new DeviceSlider(sim.getController(), new ModifierGeneral(sim, "epsilon"));
        epsilonSlider.setUnit(Kelvin.UNIT);
        epsilonSlider.setPrecision(1);
        epsilonSlider.setMinimum(100);
        epsilonSlider.setMaximum(250);
        epsilonSlider.setShowValues(true);
        epsilonSlider.setEditValues(true);
        epsilonSlider.setShowBorder(true);
        epsilonSlider.setLabel("epsilon (K)");
        epsilonSlider.setPostAction(resetMCMoves);
        add(epsilonSlider);
        epsilonSlider.doUpdate();

        DeviceSlider momentSlider = new DeviceSlider(sim.getController());
        momentSlider.setUnit(new CompoundUnit(new Unit[]{Debye.UNIT, Angstrom.UNIT}, new double[]{1,1}));
        momentSlider.setMinimum(3);
        momentSlider.setMaximum(6);
        momentSlider.setPrecision(2);
        momentSlider.setModifier(new ModifierGeneral(sim, "moment"));
        momentSlider.setShowValues(true);
        momentSlider.setEditValues(true);
        momentSlider.setShowBorder(true);
        momentSlider.setLabel("Quadrupole Moment (Debye A)");
        momentSlider.doUpdate();
        momentSlider.setPostAction(resetMCMoves);
        add(momentSlider);

        MeterDensity meterDensityLiquid = new MeterDensity(sim.getSpace());
        meterDensityLiquid.setBox(sim.boxLiquid);
        DataFork fork = new DataFork();
        DataPump pumpLiquidDensity = new DataPump(meterDensityLiquid, fork);
        getController().getDataStreamPumps().add(pumpLiquidDensity);
        AccumulatorAverageCollapsing avgLiquidDensity = new AccumulatorAverageCollapsing();
        avgLiquidDensity.setPushInterval(5);
        fork.addDataSink(avgLiquidDensity);
//        AccumulatorHistogram histogramLiquidDensity = new AccumulatorHistogram(new HistogramExpanding(0.01));
//        fork.addDataSink(histogramLiquidDensity);
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
//        AccumulatorHistogram histogramVaporDensity = new AccumulatorHistogram(new HistogramExpanding(0.001));
//        fork.addDataSink(histogramVaporDensity);
        AccumulatorHistory historyDensityVapor = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyDensityVapor.setTimeDataSource(stepCounter);
        fork.addDataSink(historyDensityVapor);
        IntegratorListenerAction pumpVaporDensityListener = new IntegratorListenerAction(pumpVaporDensity);
        sim.integratorVapor.getEventManager().addListener(pumpVaporDensityListener);
        pumpVaporDensityListener.setInterval(100);
        

//        DisplayPlot liquidDensityHistogramPlot = new DisplayPlot();
//        liquidDensityHistogramPlot.setLabel("Liquid Density");
//        add(liquidDensityHistogramPlot);
//        histogramLiquidDensity.addDataSink(liquidDensityHistogramPlot.getDataSet().makeDataSink());
//        histogramLiquidDensity.setPushInterval(100);
        DisplayTextBoxesCAE liquidDensityDisplay = new DisplayTextBoxesCAE();
        liquidDensityDisplay.setUnit(new UnitRatio(Mole.UNIT, Liter.UNIT));
        liquidDensityDisplay.setAccumulator(avgLiquidDensity);
        liquidDensityDisplay.setLabel("Liquid Density (mol/L)");
        add(liquidDensityDisplay);
        DisplayPlot liquidDensityHistoryPlot = new DisplayPlot();
        liquidDensityHistoryPlot.setUnit(new UnitRatio(Mole.UNIT, Liter.UNIT));
        liquidDensityHistoryPlot.setLabel("Liquid Density");
        add(liquidDensityHistoryPlot);
        historyDensityLiquid.addDataSink(liquidDensityHistoryPlot.getDataSet().makeDataSink());
        liquidDensityHistoryPlot.setLegend(new DataTag[]{historyDensityLiquid.getTag()}, "Density (mol/L)");
        historyDensityLiquid.setPushInterval(1);

//        DisplayPlot vaporDensityHistogramPlot = new DisplayPlot();
//        vaporDensityHistogramPlot.setLabel("Vapor Density");
//        add(vaporDensityHistogramPlot);
//        histogramVaporDensity.addDataSink(vaporDensityHistogramPlot.getDataSet().makeDataSink());
//        histogramVaporDensity.setPushInterval(100);
        DisplayTextBoxesCAE vaporDensityDisplay = new DisplayTextBoxesCAE();
        vaporDensityDisplay.setUnit(new UnitRatio(Mole.UNIT, Liter.UNIT));
        vaporDensityDisplay.setAccumulator(avgVaporDensity);
        vaporDensityDisplay.setLabel("Vapor Density (mol/L)");
        add(vaporDensityDisplay);
        DisplayPlot vaporDensityHistoryPlot = new DisplayPlot();
        vaporDensityHistoryPlot.setUnit(new UnitRatio(Mole.UNIT, Liter.UNIT));
        vaporDensityHistoryPlot.setLegend(new DataTag[]{historyDensityVapor.getTag()}, "Density (mol/L)");
        vaporDensityHistoryPlot.setLabel("Vapor Density");
        add(vaporDensityHistoryPlot);
        historyDensityVapor.addDataSink(vaporDensityHistoryPlot.getDataSet().makeDataSink());
        historyDensityVapor.setPushInterval(1);
//      historyDensityLiquid.addDataSink(densityHistoryPlot.getDataSet().makeDataSink());
//      historyDensityLiquid.setPushInterval(1);
//      historyDensityVapor.addDataSink(densityHistoryPlot.getDataSet().makeDataSink());
//      historyDensityVapor.setPushInterval(1);
//      densityHistoryPlot.setLegend(new DataTag[]{meterDensityLiquid.getTag()}, "Liquid");
//      densityHistoryPlot.setLegend(new DataTag[]{meterDensityVapor.getTag()}, "Vapor");

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
            
            DisplayPlot nMoleculesLiquidHistoryPlot = new DisplayPlot();
            nMoleculesLiquidHistoryPlot.setLabel("# of Liquid Atoms");
            add(nMoleculesLiquidHistoryPlot);
            historyNMoleculesLiquid.addDataSink(nMoleculesLiquidHistoryPlot.getDataSet().makeDataSink());
            DisplayPlot nMoleculesVaporHistoryPlot = new DisplayPlot();
            nMoleculesVaporHistoryPlot.setLabel("# of Vapor Atoms");
            add(nMoleculesVaporHistoryPlot);
            historyNMoleculesVapor.addDataSink(nMoleculesVaporHistoryPlot.getDataSet().makeDataSink());
        }
        
        MeterPressure meterPressureLiquid = new MeterPressure(sim.getSpace());
        meterPressureLiquid.setIntegrator(sim.integratorLiquid);
        fork = new DataFork();
        DataPump pumpPressureLiquid = new DataPump(meterPressureLiquid, fork);
        getController().getDataStreamPumps().add(pumpPressureLiquid);
        AccumulatorAverageCollapsing avgPressureLiquid = new AccumulatorAverageCollapsing();
        avgPressureLiquid.setPushInterval(1);
        fork.addDataSink(avgPressureLiquid);
        AccumulatorHistory historyPressureLiquid = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyPressureLiquid.setTimeDataSource(stepCounter);
        fork.addDataSink(historyPressureLiquid);
        historyPressureLiquid.setPushInterval(1);
        
        MeterPressure meterPressureVapor = new MeterPressure(sim.getSpace());
        meterPressureVapor.setIntegrator(sim.integratorVapor);
        fork = new DataFork();
        DataPump pumpPressureVapor = new DataPump(meterPressureVapor, fork);
        getController().getDataStreamPumps().add(pumpPressureVapor);
        AccumulatorAverageCollapsing avgPressureVapor = new AccumulatorAverageCollapsing();
        avgPressureVapor.setPushInterval(1);
        fork.addDataSink(avgPressureVapor);
        AccumulatorHistory historyPressureVapor = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyPressureVapor.setTimeDataSource(stepCounter);
        fork.addDataSink(historyPressureVapor);
        historyPressureVapor.setPushInterval(1);

        
        DisplayPlot plotHistoryPressure = new DisplayPlot();
        plotHistoryPressure.setLabel("Pressure History");
        plotHistoryPressure.setUnit(new PrefixedUnit(Prefix.MEGA, Pascal.UNIT));

        DisplayTextBoxesCAE displayPressureLiquid = new DisplayTextBoxesCAE();
        displayPressureLiquid.setUnit(new PrefixedUnit(Prefix.MEGA, Pascal.UNIT));
        displayPressureLiquid.setLabel("Liquid pressure (MPa)");
        displayPressureLiquid.setAccumulator(avgPressureLiquid);
        add(displayPressureLiquid);
        historyPressureLiquid.addDataSink(plotHistoryPressure.getDataSet().makeDataSink());
        plotHistoryPressure.setLegend(new DataTag[]{meterPressureLiquid.getTag()}, "Liquid (MPa)");
        IntegratorListenerAction pumpPressureLiquidListener = new IntegratorListenerAction(pumpPressureLiquid);
        sim.integratorLiquid.getEventManager().addListener(pumpPressureLiquidListener);
        pumpPressureLiquidListener.setInterval(500);
        
        DisplayTextBoxesCAE displayPressureVapor = new DisplayTextBoxesCAE();
        displayPressureVapor.setUnit(new PrefixedUnit(Prefix.MEGA, Pascal.UNIT));
        displayPressureVapor.setLabel("Vapor pressure (MPa)");
        displayPressureVapor.setAccumulator(avgPressureVapor);
        add(displayPressureVapor);
        historyPressureVapor.addDataSink(plotHistoryPressure.getDataSet().makeDataSink());
        plotHistoryPressure.setLegend(new DataTag[]{meterPressureVapor.getTag()}, "Vapor (MPa)");
        add(plotHistoryPressure);
        IntegratorListenerAction pumpPressureVaporListener = new IntegratorListenerAction(pumpPressureVapor);
        sim.integratorVapor.getEventManager().addListener(pumpPressureVaporListener);
        pumpPressureVaporListener.setInterval(500);
    }
    
    public static void main(String[] args) {
    	VLESim sim = new VLESim();
        VLE vle = new VLE(sim, sim.getSpace());
        vle.makeAndDisplayFrame();
    }

    public static class Applet extends javax.swing.JApplet {
        public void init() {
            getRootPane().putClientProperty(
                    "defeatSystemEventQueueCheck", Boolean.TRUE);
            VLESim sim = new VLESim();
            VLE vle = new VLE(sim, sim.getSpace());
            getContentPane().add(vle.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}
