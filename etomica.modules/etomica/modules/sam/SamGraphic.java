package etomica.modules.sam;

import etomica.api.IAction;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.modifier.ModifierGeneral;
import etomica.units.Bar;
import etomica.units.Pixel;
import etomica.util.HistoryCollapsingAverage;

public class SamGraphic extends SimulationGraphic {
    
    public SamGraphic(final Sam sim) {
        super(sim, SimulationGraphic.TABBED_PANE, "SAM", sim.getSpace());
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
        sim.integrator.setActionInterval(getPaintAction(sim.box), 1);
        getDisplayBox(sim.box).setShowBoundary(false);
        
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).addPlane(new WallPlane(sim.wallPotential));

        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController());
        thermoSlider.setIntegrator(sim.integrator);
        add(thermoSlider);
        
        DisplayTimer timer = new DisplayTimer(sim.integrator);
        add(timer);
        MeterTemperature thermometer = new MeterTemperature(sim, sim.box, 3);
        DisplayTextBox temperatureDisplay = new DisplayTextBox();
        DataPump pump = new DataPump(thermometer, temperatureDisplay);
        sim.integrator.addIntervalAction(pump);
        add(temperatureDisplay);

        ModifierGeneral wallPositionModifier = new ModifierGeneral(sim.wallPotential, "wallPosition");
        DeviceSlider wallPositionSlider = new DeviceSlider(sim.getController(), wallPositionModifier);
        wallPositionSlider.setPrecision(1);
        wallPositionSlider.setShowBorder(true);
        wallPositionSlider.setLabel("Wall position");
        wallPositionSlider.setMinimum(5);
        wallPositionSlider.setMaximum(20);
        wallPositionSlider.setShowValues(true);
        add(wallPositionSlider);

        DeviceSelector primitiveSelect = new DeviceSelector(sim.getController());
        primitiveSelect.addOption("1", new IAction() {
            public void actionPerformed() {
                sim.config.setBasisMolecules(new BasisMonatomic(sim.getSpace()));
                sim.config.setBasisSurface(new BasisMonatomic(sim.getSpace()));
                sim.config.setCellSizeX(4);
                sim.config.setCellSizeZ(4);
                sim.config.initializeCoordinates(sim.box);
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {}
            }
        });
        add(primitiveSelect);
        primitiveSelect.setLabel("Basis size");
        
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
        MeterKineticEnergy meterKE = new MeterKineticEnergy();
        meterKE.setBox(sim.box);
        AccumulatorHistory historyKE = new AccumulatorHistory();
        historyKE.setTimeDataSource(timeCounter);
        pump = new DataPump(meterKE, historyKE);
        getController().getDataStreamPumps().add(pump);
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 10);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.integrator.getPotential());
        meterPE.setBox(sim.box);
        AccumulatorHistory historyPE = new AccumulatorHistory();
        historyPE.setTimeDataSource(timeCounter);
        pump = new DataPump(meterPE, historyPE);
        getController().getDataStreamPumps().add(pump);
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 10);
        MeterEnergy meterEnergy = new MeterEnergy(sim.integrator.getPotential(), sim.box);
        AccumulatorHistory historyE = new AccumulatorHistory();
        historyE.setTimeDataSource(timeCounter);
        pump = new DataPump(meterEnergy, historyE);
        getController().getDataStreamPumps().add(pump);
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 10);

        PotentialCalculationForceSumWall forceSum = new PotentialCalculationForceSumWall(sim.wallPotential);
        sim.integrator.setForceSum(forceSum);
        MeterWallPressure wallPressure = new MeterWallPressure(forceSum);
        wallPressure.setBox(sim.box);
        DataFork fork = new DataFork();
        pump = new DataPump(wallPressure, fork);
        getController().getDataStreamPumps().add(pump);
        AccumulatorAverageCollapsing wallPressureAvg = new AccumulatorAverageCollapsing();
        fork.addDataSink(wallPressureAvg);
        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
        pressureDisplay.setUnit(Bar.UNIT);
        pressureDisplay.setAccumulator(wallPressureAvg);
        sim.integrator.addIntervalAction(pump);
        wallPressureAvg.setPushInterval(10);
        add(pressureDisplay);

        AccumulatorHistory wallPressureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        fork.addDataSink(wallPressureHistory);
        wallPressureHistory.setPushInterval(10);
        DisplayPlot pressurePlot = new DisplayPlot();
        pressurePlot.setUnit(Bar.UNIT);
        wallPressureHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
        pressurePlot.setLabel("Wall Pressure");
        add(pressurePlot);
        
        
        DisplayPlot plot = new DisplayPlot();
        historyKE.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{meterKE.getTag()}, "KE");
        historyPE.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{meterPE.getTag()}, "PE");
        historyE.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{meterEnergy.getTag()}, "E");
        plot.setLabel("Energy");

        add(plot);
    }
    
    public static void main(String[] args) {
        Sam sim = new Sam();
        SamGraphic graphic = new SamGraphic(sim);
        graphic.makeAndDisplayFrame();
//        sim.activityIntegrate.setSleepPeriod(10);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
            getRootPane().putClientProperty(
                            "defeatSystemEventQueueCheck", Boolean.TRUE);
            SamGraphic samGraphic = new SamGraphic(new Sam());

            getContentPane().add(samGraphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}
