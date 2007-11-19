package etomica.modules.vle;

import etomica.action.Action;
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
import etomica.modifier.ModifierGeneral;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pascal;
import etomica.units.Pixel;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.UnitRatio;
import etomica.util.HistoryCollapsingAverage;

public class VLE extends SimulationGraphic {

    private final static String APP_NAME = "Virial / VLE";
    private final static int REPAINT_INTERVAL = 200;

    public VLE(final VLESim sim) {
        super(sim, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        getDisplayBox(sim.boxLiquid).setPixelUnit(new Pixel(8));
        getDisplayBox(sim.boxVapor).setPixelUnit(new Pixel(8));

        DeviceThermoSliderGEMC thermoSlider = new DeviceThermoSliderGEMC(sim.getController());
        thermoSlider.setUnit(Kelvin.UNIT);
        thermoSlider.setPrecision(1);
        thermoSlider.setIntegrators(sim.integratorLiquid, sim.integratorVapor);
        thermoSlider.setIsothermal();
        thermoSlider.setMinimum(200);
        thermoSlider.setMaximum(400);
        thermoSlider.doUpdate();
        thermoSlider.setIsothermalButtonsVisibility(false);
        add(thermoSlider);
        
        DeviceSlider sigmaSlider = new DeviceSlider(sim.getController(), new ModifierGeneral(sim, "sigma"));
        sigmaSlider.setPrecision(2);
        sigmaSlider.setMinimum(2);
        sigmaSlider.setMaximum(5);
        sigmaSlider.setShowValues(true);
        sigmaSlider.setEditValues(true);
        sigmaSlider.setShowBorder(true);
        sigmaSlider.setLabel("sigma (A)");
        add(sigmaSlider);
        sigmaSlider.doUpdate();
        sigmaSlider.setPostAction(new Action() {
            public void actionPerformed() {
                getPaintAction(sim.boxLiquid).actionPerformed();
                getPaintAction(sim.boxVapor).actionPerformed();
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
        add(epsilonSlider);
        epsilonSlider.doUpdate();

        DeviceSlider momentSlider = new DeviceSlider(sim.getController(), new ModifierGeneral(sim, "moment"));
        momentSlider.setUnit(Debye.UNIT);
        momentSlider.setPrecision(2);
        momentSlider.setMinimum(3);
        momentSlider.setMaximum(6);
        momentSlider.setShowValues(true);
        momentSlider.setEditValues(true);
        momentSlider.setShowBorder(true);
        momentSlider.setLabel("moment (Debye)");
        momentSlider.doUpdate();
        add(momentSlider);

        MeterDensity meterDensityLiquid = new MeterDensity(sim.getSpace());
        meterDensityLiquid.setBox(sim.boxLiquid);
        DataFork fork = new DataFork();
        DataPump pumpLiquidDensity = new DataPump(meterDensityLiquid, fork);
        getController().getDataStreamPumps().add(pumpLiquidDensity);
        AccumulatorAverageCollapsing avgLiquidDensity = new AccumulatorAverageCollapsing();
        fork.addDataSink(avgLiquidDensity);
//        AccumulatorHistogram histogramLiquidDensity = new AccumulatorHistogram(new HistogramExpanding(0.01));
//        fork.addDataSink(histogramLiquidDensity);
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integratorGEMC);
        AccumulatorHistory historyDensityLiquid = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyDensityLiquid.setTimeDataSource(stepCounter);
        fork.addDataSink(historyDensityLiquid);
        sim.integratorLiquid.addIntervalAction(pumpLiquidDensity);
        sim.integratorLiquid.setActionInterval(pumpLiquidDensity, 100);
        
        MeterDensity meterDensityVapor = new MeterDensity(sim.getSpace());
        meterDensityVapor.setBox(sim.boxVapor);
        fork = new DataFork();
        DataPump pumpVaporDensity = new DataPump(meterDensityVapor, fork);
        getController().getDataStreamPumps().add(pumpVaporDensity);
        AccumulatorAverageCollapsing avgVaporDensity = new AccumulatorAverageCollapsing();
        fork.addDataSink(avgVaporDensity);
//        AccumulatorHistogram histogramVaporDensity = new AccumulatorHistogram(new HistogramExpanding(0.001));
//        fork.addDataSink(histogramVaporDensity);
        AccumulatorHistory historyDensityVapor = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyDensityVapor.setTimeDataSource(stepCounter);
        fork.addDataSink(historyDensityVapor);
        sim.integratorVapor.addIntervalAction(pumpVaporDensity);
        sim.integratorVapor.setActionInterval(pumpVaporDensity, 100);
        

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
        
        MeterNMolecules meterNMoleculesLiquid = new MeterNMolecules();
        meterNMoleculesLiquid.setBox(sim.boxLiquid);
        AccumulatorHistory historyNMoleculesLiquid = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyNMoleculesLiquid.setTimeDataSource(stepCounter);
        DataPump pumpNMoleculesLiquid = new DataPump(meterNMoleculesLiquid, historyNMoleculesLiquid);
        sim.integratorLiquid.addIntervalAction(pumpNMoleculesLiquid);
        sim.integratorLiquid.setActionInterval(pumpNMoleculesLiquid, 100);
        getController().getDataStreamPumps().add(pumpNMoleculesLiquid);
        MeterNMolecules meterNMoleculesVapor = new MeterNMolecules();
        meterNMoleculesVapor.setBox(sim.boxVapor);
        AccumulatorHistory historyNMoleculesVapor = new AccumulatorHistory(new HistoryCollapsingAverage(100));
        historyNMoleculesVapor.setTimeDataSource(stepCounter);
        DataPump pumpNMoleculesVapor = new DataPump(meterNMoleculesVapor, historyNMoleculesVapor);
        sim.integratorVapor.addIntervalAction(pumpNMoleculesVapor);
        sim.integratorVapor.setActionInterval(pumpNMoleculesVapor, 100);
        getController().getDataStreamPumps().add(pumpNMoleculesVapor);
        
        DisplayPlot nMoleculesLiquidHistoryPlot = new DisplayPlot();
        nMoleculesLiquidHistoryPlot.setLabel("# of Liquid Atoms");
        add(nMoleculesLiquidHistoryPlot);
        historyNMoleculesLiquid.addDataSink(nMoleculesLiquidHistoryPlot.getDataSet().makeDataSink());
        DisplayPlot nMoleculesVaporHistoryPlot = new DisplayPlot();
        nMoleculesVaporHistoryPlot.setLabel("# of Vapor Atoms");
        add(nMoleculesVaporHistoryPlot);
        historyNMoleculesVapor.addDataSink(nMoleculesVaporHistoryPlot.getDataSet().makeDataSink());
        
        MeterPressure meterPressureLiquid = new MeterPressure(sim.getSpace());
        meterPressureLiquid.setIntegrator(sim.integratorLiquid);
        fork = new DataFork();
        DataPump pumpPressureLiquid = new DataPump(meterPressureLiquid, fork);
        getController().getDataStreamPumps().add(pumpPressureLiquid);
        AccumulatorAverageCollapsing avgPressureLiquid = new AccumulatorAverageCollapsing();
        avgPressureLiquid.setPushInterval(10);
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
        avgPressureVapor.setPushInterval(10);
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
        sim.integratorLiquid.addIntervalAction(pumpPressureLiquid);
        sim.integratorLiquid.setActionInterval(pumpPressureLiquid, 500);
        
        DisplayTextBoxesCAE displayPressureVapor = new DisplayTextBoxesCAE();
        displayPressureVapor.setUnit(new PrefixedUnit(Prefix.MEGA, Pascal.UNIT));
        displayPressureVapor.setLabel("Vapor pressure (MPa)");
        displayPressureVapor.setAccumulator(avgPressureVapor);
        add(displayPressureVapor);
        historyPressureVapor.addDataSink(plotHistoryPressure.getDataSet().makeDataSink());
        plotHistoryPressure.setLegend(new DataTag[]{meterPressureVapor.getTag()}, "Vapor (MPa)");
        add(plotHistoryPressure);
        sim.integratorVapor.addIntervalAction(pumpPressureVapor);
        sim.integratorVapor.setActionInterval(pumpPressureVapor, 500);
    }
    
    public static void main(String[] args) {
        VLE vle = new VLE(new VLESim());
        vle.makeAndDisplayFrame();
    }
}
