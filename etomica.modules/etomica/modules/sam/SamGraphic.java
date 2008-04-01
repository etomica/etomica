package etomica.modules.sam;

import java.awt.Color;

import etomica.api.IAction;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
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
import etomica.graphics.ColorSchemeByType;
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
import etomica.modifier.Modifier;
import etomica.units.Angle;
import etomica.units.Bar;
import etomica.units.Degree;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Pixel;
import etomica.util.HistoryCollapsingAverage;

public class SamGraphic extends SimulationGraphic {
    
    public SamGraphic(final Sam sim) {
        super(sim, SimulationGraphic.TABBED_PANE, "SAM", sim.getSpace());
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
        sim.integrator.setActionInterval(getPaintAction(sim.box), 1);
        getDisplayBox(sim.box).setShowBoundary(false);
        
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).addPlane(new WallPlane(sim.wallPotential));
        
        ((ColorSchemeByType)getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getCH2Type(), new Color(200, 200, 200));
        ((ColorSchemeByType)getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getCH3Type(), new Color(220, 220, 220));
        ((ColorSchemeByType)getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getSulfurType(), new Color(255, 200, 50));
        ((ColorSchemeByType)getDisplayBox(sim.box).getColorScheme()).setColor(sim.speciesSurface.getLeafType(), new Color(218, 165, 32));

        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController());
        thermoSlider.setIntegrator(sim.integrator);
        thermoSlider.setMaximum(500);
        add(thermoSlider);
        
        DisplayTimer timer = new DisplayTimer(sim.integrator);
        add(timer);
        MeterTemperature thermometer = new MeterTemperature(sim, sim.box, 3);
        DisplayTextBox temperatureDisplay = new DisplayTextBox();
        DataPump pump = new DataPump(thermometer, temperatureDisplay);
        sim.integrator.addIntervalAction(pump);
        add(temperatureDisplay);

        Modifier wallPositionModifier = new ModifierWallPosition(sim);
        DeviceSlider wallPositionSlider = new DeviceSlider(sim.getController(), wallPositionModifier);
        wallPositionSlider.setPrecision(1);
        wallPositionSlider.setShowBorder(true);
        wallPositionSlider.setLabel("Wall position");
        wallPositionSlider.setMinimum(15);
        wallPositionSlider.setMaximum(30);
        wallPositionSlider.setShowValues(true);
        wallPositionSlider.setPostAction(getPaintAction(sim.box));
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
        pressureDisplay.setLabel("Wall Pressure (bar)");
        pressureDisplay.setAccumulator(wallPressureAvg);
        sim.integrator.addIntervalAction(pump);
        wallPressureAvg.setPushInterval(10);
        add(pressureDisplay);

        AccumulatorHistory wallPressureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        wallPressureHistory.setTimeDataSource(timeCounter);
        fork.addDataSink(wallPressureHistory);
        wallPressureHistory.setPushInterval(10);
        DisplayPlot pressurePlot = new DisplayPlot();
        pressurePlot.setUnit(Bar.UNIT);
        wallPressureHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
        pressurePlot.setLabel("Pressure");
        pressurePlot.setLegend(new DataTag[]{wallPressure.getTag()}, "Wall Pressure (bar)");
        add(pressurePlot);
        
        MeterTilt meterTilt = new MeterTilt(space, sim.species);
        meterTilt.setBox(sim.box);
        fork = new DataFork();
        pump = new DataPump(meterTilt, fork);
        AccumulatorAverageCollapsing tiltAvg = new AccumulatorAverageCollapsing();
        fork.addDataSink(tiltAvg);
        DisplayTextBoxesCAE tiltDisplay = new DisplayTextBoxesCAE();
        tiltDisplay.setAccumulator(tiltAvg);
        tiltDisplay.setUnit(Degree.UNIT);
        tiltDisplay.setLabel("Tilt angle");
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 10);
        add(tiltDisplay);
        
        AccumulatorHistory tiltHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        tiltHistory.setTimeDataSource(timeCounter);
        fork.addDataSink(tiltHistory);
        DisplayPlot tiltPlot = new DisplayPlot();
        tiltHistory.setDataSink(tiltPlot.getDataSet().makeDataSink());
        tiltPlot.setUnit(Degree.UNIT);
        tiltPlot.setLabel("Tilt");
        add(tiltPlot);
        
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
    
    public static class ModifierWallPosition implements Modifier {
        private final Sam sim;
        protected double surfacePosition;

        public ModifierWallPosition(Sam sim) {
            this.sim = sim;
            update();
        }

        public void update() {
            surfacePosition = ((IAtomPositioned)((IMolecule)sim.box.getMoleculeList(sim.speciesSurface).getAtom(0)).getChildList().getAtom(0)).getPosition().x(1);
        }

        public void setValue(double newValue) {
            double maxAtomPos = -100;
            IAtomSet leafList = sim.box.getLeafList();
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtomPositioned atom = (IAtomPositioned)leafList.getAtom(i);
                double atomPos = atom.getPosition().x(1);
                if (atomPos > maxAtomPos) {
                    maxAtomPos = atomPos;
                }
            }
            if (newValue+surfacePosition < maxAtomPos + 4) {
                throw new RuntimeException("Slow down tiger!");
            }
            sim.wallPotential.setWallPosition(newValue + surfacePosition);
        }

        public double getValue() {
            return sim.wallPotential.getWallPosition() - surfacePosition;
        }

        public Dimension getDimension() {
            return Length.DIMENSION;
        }

        public String getLabel() {
            return "Wall Position";
        }
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
