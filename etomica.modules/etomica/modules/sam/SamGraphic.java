package etomica.modules.sam;

import java.awt.Color;
import java.awt.GridBagLayout;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.action.SimulationRestart;
import etomica.api.IAction;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.units.Bar;
import etomica.units.Degree;
import etomica.units.Dimension;
import etomica.units.Kelvin;
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

        ((SimulationRestart)getController().getReinitButton().getAction()).setConfiguration(sim.config);
        getController().getReinitButton().setPostAction(getPaintAction(sim.box));
        
        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController());
        thermoSlider.setIsothermalButtonsVisibility(false);
        thermoSlider.setUnit(Kelvin.UNIT);
        thermoSlider.setIntegrator(sim.integrator);
        thermoSlider.setMaximum(500);
        add(thermoSlider);
        
        DisplayTimer timer = new DisplayTimer(sim.integrator);
        add(timer);
        MeterTemperature thermometer = new MeterTemperature(sim, sim.box, 3);
        DisplayTextBox temperatureDisplay = new DisplayTextBox();
        temperatureDisplay.setUnit(Kelvin.UNIT);
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
        
        JTabbedPane conformationTabs = new JTabbedPane();
        getPanel().controlPanel.add(conformationTabs, SimulationPanel.getVertGBC());

        JPanel chain1 = new JPanel(new GridBagLayout());
        conformationTabs.add(chain1, "chain 1");
        ModifierGeneral thetaModifier = new ModifierGeneral(sim, "chainTheta");
        DeviceSlider thetaSlider = new DeviceSlider(sim.getController(), thetaModifier);
        thetaSlider.setShowBorder(true);
        thetaSlider.setUnit(Degree.UNIT);
        thetaSlider.setLabel("Theta");
        thetaSlider.setMinimum(0);
        thetaSlider.setMaximum(50);
        thetaSlider.setShowValues(true);
        IAction reinitAction = new IAction() {
            public void actionPerformed() {
                sim.config.initializeCoordinates(sim.box);
                sim.findTetherBonds();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {}
                getPaintAction(sim.box).actionPerformed();
            }
        };
        thetaSlider.setPostAction(reinitAction);
        chain1.add(thetaSlider.graphic(), SimulationPanel.getVertGBC());

        ModifierGeneral psiModifier = new ModifierGeneral(sim, "chainPsi");
        DeviceSlider psiSlider = new DeviceSlider(sim.getController(), psiModifier);
        psiSlider.setShowBorder(true);
        psiSlider.setUnit(Degree.UNIT);
        psiSlider.setLabel("Psi");
        psiSlider.setMinimum(0);
        psiSlider.setMaximum(360);
        psiSlider.setShowValues(true);
        psiSlider.setPostAction(reinitAction);
        chain1.add(psiSlider.graphic(), SimulationPanel.getVertGBC());

        ModifierGeneral phiModifier = new ModifierGeneral(sim, "chainPhi");
        DeviceSlider phiSlider = new DeviceSlider(sim.getController(), phiModifier);
        phiSlider.setShowBorder(true);
        phiSlider.setUnit(Degree.UNIT);
        phiSlider.setLabel("Phi");
        phiSlider.setMinimum(0);
        phiSlider.setMaximum(360);
        phiSlider.setShowValues(true);
        phiSlider.setPostAction(reinitAction);
        chain1.add(phiSlider.graphic(), SimulationPanel.getVertGBC());

        JPanel chain2 = new JPanel(new GridBagLayout());
        conformationTabs.add(chain2, "chain 2");
        ModifierGeneral secondaryThetaModifier = new ModifierGeneral(sim, "secondaryChainTheta");
        DeviceSlider secondaryThetaSlider = new DeviceSlider(sim.getController(), secondaryThetaModifier);
        secondaryThetaSlider.setShowBorder(true);
        secondaryThetaSlider.setUnit(Degree.UNIT);
        secondaryThetaSlider.setLabel("Theta");
        secondaryThetaSlider.setMinimum(0);
        secondaryThetaSlider.setMaximum(50);
        secondaryThetaSlider.setShowValues(true);
        secondaryThetaSlider.setPostAction(reinitAction);
        chain2.add(secondaryThetaSlider.graphic(), SimulationPanel.getVertGBC());

        ModifierGeneral secondaryPsiModifier = new ModifierGeneral(sim, "secondaryChainPsi");
        DeviceSlider secondaryPsiSlider = new DeviceSlider(sim.getController(), secondaryPsiModifier);
        secondaryPsiSlider.setShowBorder(true);
        secondaryPsiSlider.setUnit(Degree.UNIT);
        secondaryPsiSlider.setLabel("Psi");
        secondaryPsiSlider.setMinimum(0);
        secondaryPsiSlider.setMaximum(360);
        secondaryPsiSlider.setShowValues(true);
        secondaryPsiSlider.setPostAction(reinitAction);
        chain2.add(secondaryPsiSlider.graphic(), SimulationPanel.getVertGBC());

        ModifierGeneral secondaryPhiModifier = new ModifierGeneral(sim, "secondaryChainPhi");
        DeviceSlider secondaryPhiSlider = new DeviceSlider(sim.getController(), secondaryPhiModifier);
        secondaryPhiSlider.setShowBorder(true);
        secondaryPhiSlider.setUnit(Degree.UNIT);
        secondaryPhiSlider.setLabel("Phi");
        secondaryPhiSlider.setMinimum(0);
        secondaryPhiSlider.setMaximum(360);
        secondaryPhiSlider.setShowValues(true);
        secondaryPhiSlider.setPostAction(reinitAction);
        chain2.add(secondaryPhiSlider.graphic(), SimulationPanel.getVertGBC());

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
        DataSplitter tiltSplitter = new DataSplitter(); 
        pump = new DataPump(meterTilt, tiltSplitter);
        getController().getDataStreamPumps().add(pump);
        fork = new DataFork();
        tiltSplitter.setDataSink(0, fork);
        AccumulatorAverageCollapsing tiltAvg = new AccumulatorAverageCollapsing();
        tiltAvg.setPushInterval(1);
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
        tiltPlot.setLegend(new DataTag[]{tiltHistory.getTag()}, "net tilt");
        tiltPlot.setUnit(Degree.UNIT);
        tiltPlot.setLabel("Tilt");
        add(tiltPlot);

        AccumulatorHistory tilt2History = new AccumulatorHistory(new HistoryCollapsingAverage());
        tilt2History.setTimeDataSource(timeCounter);
        tiltSplitter.setDataSink(1, tilt2History);
        tilt2History.setDataSink(tiltPlot.getDataSet().makeDataSink());
        tiltPlot.setLegend(new DataTag[]{tilt2History.getTag()}, "avg tilt");
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
