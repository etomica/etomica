/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryScrolling;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2LennardJones;
import etomica.units.Bar;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.dimensions.Angle;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import javax.swing.*;
import java.awt.*;

public class SamGraphic extends SimulationGraphic {
    
    public final ActionMoveWall moveWallAction;
    public final DeviceButton wallButton;
    public final DeviceButton corrugationButton;
    public final DeviceSlider thetaSlider;
    public final DeviceSlider[] phiSlider;
    public final DeviceSlider numCellsSlider;
    
    public SamGraphic(final Sam sim) {
        super(sim, SimulationGraphic.TABBED_PANE, "SAM");
        DisplayBox displayBox = getDisplayBox(sim.box);
        displayBox.setPixelUnit(new Pixel(8));
        setPaintInterval(sim.box, 1);
        displayBox.setShowBoundary(false);
        
        final IAction resetAction = getController().getSimRestart().getDataResetAction();
        
        ((DisplayBoxCanvasG3DSys)displayBox.canvas).addPlane(new WallPlane(space, sim.wallPotential));
        
        ((ColorSchemeByType)displayBox.getColorScheme()).setColor(sim.species.getCH2Type(), new Color(190, 190, 190));
        ((ColorSchemeByType)displayBox.getColorScheme()).setColor(sim.species.getCH3Type(), new Color(230, 230, 230));
        ((ColorSchemeByType)displayBox.getColorScheme()).setColor(sim.species.getSulfurType(), new Color(255, 200, 50));
        ((ColorSchemeByType)displayBox.getColorScheme()).setColor(sim.speciesSurface.getLeafType(), new Color(218, 165, 32));
        
        DiameterHashByType diameterHash = (DiameterHashByType)displayBox.getDiameterHash();
        diameterHash.setDiameter(sim.speciesSurface.getLeafType(), 3.0);
        diameterHash.setDiameter(sim.species.getCH2Type(), sim.p2CH2.getSigma());
        diameterHash.setDiameter(sim.species.getCH3Type(), sim.p2CH3.getSigma());
        diameterHash.setDiameter(sim.species.getSulfurType(), sim.p2S.getSigma());


        getController().getReinitButton().setPostAction(getPaintAction(sim.box));
        
        final IAction moveWallToggle = new IAction() {
            public void actionPerformed() {
                if (moveWallAction.enabled) {
                    moveWallAction.enabled = false;
                    wallButton.setLabel("Lower wall");
                    thetaSlider.setEnabled(true);
                    for (int i=0; i<4; i++) {
                        phiSlider[i].setEnabled(true);
                    }
                    numCellsSlider.setEnabled(true);
                    corrugationButton.getButton().setEnabled(true);
                }
                else {
                    moveWallAction.enabled = true;
                    thetaSlider.setEnabled(false);
                    for (int i=0; i<4; i++) {
                        phiSlider[i].setEnabled(false);
                    }
                    numCellsSlider.setEnabled(false);
                    corrugationButton.getButton().setEnabled(false);
                    wallButton.setLabel("Stop wall");
                }
            }
        };
        wallButton = new DeviceButton(sim.getController(), moveWallToggle);
        wallButton.setLabel("Lower wall");
        
        JTabbedPane moreTabs = new JTabbedPane();
        getPanel().controlPanel.add(moreTabs, SimulationPanel.getVertGBC());
        JPanel wallTab = new JPanel(new GridBagLayout());
        wallTab.add(wallButton.graphic(), SimulationPanel.getVertGBC());
        moreTabs.add(wallTab, "wall");

        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        thermoSlider.setIsothermalButtonsVisibility(false);
        thermoSlider.setUnit(Kelvin.UNIT);
        thermoSlider.setMaximum(500);
        thermoSlider.setSliderPostAction(resetAction);
        JPanel stateTab = new JPanel(new GridBagLayout());
        moreTabs.add(stateTab, "state");
        stateTab.add(thermoSlider.graphic(), SimulationPanel.getVertGBC());
        
        DisplayTimer timer = new DisplayTimer(sim.integrator);
        add(timer);
        MeterTemperature thermometer = new MeterTemperature(sim, sim.box, 3);
        DisplayTextBox temperatureDisplay = new DisplayTextBox();
        temperatureDisplay.setUnit(Kelvin.UNIT);
        DataPump pump = new DataPump(thermometer, temperatureDisplay);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        add(temperatureDisplay);

        Modifier wallPositionModifier = new ModifierWallPosition(sim);
        final DeviceSlider wallPositionSlider = new DeviceSlider(sim.getController(), wallPositionModifier);
        wallPositionSlider.setPrecision(1);
        wallPositionSlider.setShowBorder(true);
        wallPositionSlider.setLabel("Wall position");
        wallPositionSlider.setMinimum(15);
        double surfacePosition = sim.box.getMoleculeList(sim.speciesSurface).getMolecule(0).getChildList().getAtom(0).getPosition().getX(1);
        double sliderValue = sim.wallPotential.getWallPosition()-surfacePosition;
        if (sliderValue < 30) {
            wallPositionSlider.setMaximum(30);
        }
        else {
            wallPositionSlider.setMaximum(40);
            wallPositionSlider.setNMajor(5);
        }
        wallPositionSlider.setShowValues(true);
        wallPositionSlider.setPostAction(getPaintAction(sim.box));
        wallTab.add(wallPositionSlider.graphic(), SimulationPanel.getVertGBC());
        
        JTabbedPane conformationTabs = new JTabbedPane();
        getPanel().controlPanel.add(conformationTabs, SimulationPanel.getVertGBC());

        JPanel chainThetaPsi = new JPanel(new GridBagLayout());
        conformationTabs.add(chainThetaPsi, "theta, psi");
        ModifierGeneral thetaModifier = new ModifierGeneral(sim, "chainTheta");
        thetaSlider = new DeviceSlider(sim.getController(), thetaModifier);
        thetaSlider.setShowBorder(true);
        thetaSlider.setUnit(Degree.UNIT);
        thetaSlider.setLabel("Theta");
        thetaSlider.setMinimum(0);
        thetaSlider.setMaximum(50);
        thetaSlider.setShowValues(true);
        IAction reinitAction = new IAction() {
            public void actionPerformed() {
                sim.wallPotential.setWallPosition(20);
                wallPositionSlider.doUpdate();
                sim.config.initializeCoordinates(sim.box);
                sim.findTetherBonds();
                sim.integrator.resetStepCount();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {}
                getDisplayBox(sim.box).repaint();
                resetAction.actionPerformed();
            }
        };
        // SimulationRestart isn't enough, we need to re-findtetherbonds
        getController().getReinitButton().setAction(reinitAction);
        thetaSlider.setPostAction(reinitAction);
        chainThetaPsi.add(thetaSlider.graphic(), SimulationPanel.getVertGBC());

        ModifierGeneral psiModifier = new ModifierGeneral(sim, "chainPsi");
        DeviceSlider psiSlider = new DeviceSlider(sim.getController(), psiModifier);
        psiSlider.setShowBorder(true);
        psiSlider.setUnit(Degree.UNIT);
        psiSlider.setLabel("Psi");
        psiSlider.setMinimum(0);
        psiSlider.setMaximum(360);
        psiSlider.setShowValues(true);
        psiSlider.setPostAction(reinitAction);
        chainThetaPsi.add(psiSlider.graphic(), SimulationPanel.getVertGBC());

        JPanel chainPhi = null;
        phiSlider = new DeviceSlider[4];
        for (int i=0; i<4; i++) {
            if (i==0) {
                chainPhi = new JPanel(new GridBagLayout());
                conformationTabs.add(chainPhi, "phi 1,2");
            }
            else if (i==2) {
                chainPhi = new JPanel(new GridBagLayout());
                conformationTabs.add(chainPhi, "phi 3,4");
            }

            final int iChain = i;
            Modifier phiModifier = new Modifier() {
                public Dimension getDimension() {
                    return Angle.DIMENSION;
                }

                public String getLabel() {
                    return "Phi "+iChain;
                }

                public double getValue() {
                    return sim.getChainPhi(iChain);
                }

                public void setValue(double newValue) {
                    sim.setChainPhi(iChain, newValue);
                }
            };
            phiSlider[i] = new DeviceSlider(sim.getController(), phiModifier);
            phiSlider[i].setShowBorder(true);
            phiSlider[i].setUnit(Degree.UNIT);
            phiSlider[i].setLabel("Phi "+(i+1));
            phiSlider[i].setMinimum(0);
            phiSlider[i].setMaximum(360);
            phiSlider[i].setShowValues(true);
            phiSlider[i].setPostAction(reinitAction);
            chainPhi.add(phiSlider[i].graphic(), SimulationPanel.getVertGBC());
        }

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
        MeterKineticEnergy meterKE = new MeterKineticEnergy(sim.box);
        AccumulatorHistory historyKE = new AccumulatorHistory();
        historyKE.setTimeDataSource(timeCounter);
        pump = new DataPump(meterKE, historyKE);
        getController().getDataStreamPumps().add(pump);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(10);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());
        meterPE.setBox(sim.box);
        AccumulatorHistory historyPE = new AccumulatorHistory();
        historyPE.setTimeDataSource(timeCounter);
        pump = new DataPump(meterPE, historyPE);
        getController().getDataStreamPumps().add(pump);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(10);
        MeterEnergy meterEnergy = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        AccumulatorHistory historyE = new AccumulatorHistory();
        historyE.setTimeDataSource(timeCounter);
        pump = new DataPump(meterEnergy, historyE);
        getController().getDataStreamPumps().add(pump);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(10);

        Modifier nCellsModifier = new Modifier() {
            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            public String getLabel() {
                return "num cells";
            }

            public double getValue() {
                return sim.getNumXCells()/2;
            }

            public void setValue(double newValue) {
                if (newValue > 2 && newValue < 4) return;
                sim.setNumXCells(2*(int)newValue);
                sim.setNumZCells((int)newValue);
            }
        };
        numCellsSlider = new DeviceSlider(sim.getController(), nCellsModifier);
        numCellsSlider.setShowBorder(true);
        numCellsSlider.setLabel("# of cells");
        numCellsSlider.setMinimum(1);
        numCellsSlider.setMaximum(4);
        numCellsSlider.doUpdate();
        numCellsSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                getDisplayBox(sim.box).repaint();
                resetAction.actionPerformed();
            }
        });

        getPanel().controlPanel.add(numCellsSlider.graphic(), SimulationPanel.getVertGBC());
        
        MeterWallPressure wallPressure = new MeterWallPressure(sim.forceSum);
        wallPressure.setBox(sim.box);
        DataFork fork = new DataFork();
        pump = new DataPump(wallPressure, fork);
        getController().getDataStreamPumps().add(pump);
        AccumulatorAverageCollapsing wallPressureAvg = new AccumulatorAverageCollapsing();
        fork.addDataSink(wallPressureAvg);
        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
        pressureDisplay.setUnit(Bar.UNIT);
        pressureDisplay.setLabel("Wall Stress (bar)");
        pressureDisplay.setAccumulator(wallPressureAvg);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        wallPressureAvg.setPushInterval(10);
        add(pressureDisplay);

        AccumulatorHistory wallPressureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        wallPressureHistory.setTimeDataSource(timeCounter);
        fork.addDataSink(wallPressureHistory);
        wallPressureHistory.setPushInterval(10);
        DisplayPlot pressurePlot = new DisplayPlot();
        pressurePlot.setUnit(Bar.UNIT);
        wallPressureHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
        pressurePlot.setLabel("Stress");
        pressurePlot.setLegend(new DataTag[]{wallPressure.getTag()}, "Wall Stress (bar)");
        add(pressurePlot);
        
        final MeterTilt meterTilt = new MeterTilt(space, sim.species);
        meterTilt.setBox(sim.box);
        DataSplitter tiltSplitter = new DataSplitter(); 
        pump = new DataPump(meterTilt, tiltSplitter);
        getController().getDataStreamPumps().add(pump);
        DataFork tiltFork = new DataFork();
        tiltSplitter.setDataSink(0, tiltFork);
        AccumulatorAverageCollapsing tiltAvg = new AccumulatorAverageCollapsing();
        tiltAvg.setPushInterval(1);
        tiltFork.addDataSink(tiltAvg);
        DisplayTextBoxesCAE tiltDisplay = new DisplayTextBoxesCAE();
        tiltDisplay.setAccumulator(tiltAvg);
        tiltDisplay.setUnit(Degree.UNIT);
        tiltDisplay.setLabel("Tilt angle");
        pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(10);
        add(tiltDisplay);
        
        AccumulatorHistory tiltHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        tiltHistory.setTimeDataSource(timeCounter);
        tiltFork.addDataSink(tiltHistory);
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
        
        AccumulatorHistogram stressStrainHistogram = new AccumulatorHistogram(new HistogramNotSoSimple(100, new DoubleRange(0, 40)));
        DataPipe stressStrainPipe = new DataPipeStressStrain(sim);
        fork.addDataSink(stressStrainPipe);
        stressStrainPipe.setDataSink(stressStrainHistogram);
        final DisplayPlot stressStrainPlot = new DisplayPlot();
        stressStrainHistogram.setDataSink(stressStrainPlot.getDataSet().makeDataSink());
        stressStrainPlot.setLabel("Stress vs. Strain");
        stressStrainPlot.setXLabel("Wall Position");
        stressStrainPlot.setDoLegend(false);
        add(stressStrainPlot);
        
        AccumulatorHistory energyTilt = new AccumulatorHistory(new HistoryScrolling(1));
        energyTilt.setTimeDataSource(new DataSourceScalar("Tilt Angle", Angle.DIMENSION) {
            public double getDataAsScalar() {
                return meterTilt.getData().getValue(0);
            }
        });
        DataPump energyTiltPump = new DataPump(meterPE, energyTilt);
        IntegratorListenerAction energyTiltPumpListener = new IntegratorListenerAction(energyTiltPump);
        sim.integrator.getEventManager().addListener(energyTiltPumpListener);
        energyTiltPumpListener.setInterval(10);
        final DisplayPlot energyTiltPlot = new DisplayPlot();
        energyTiltPlot.setXUnit(Degree.UNIT);
        energyTilt.setDataSink(energyTiltPlot.getDataSet().makeDataSink());
        energyTiltPlot.setLabel("Energy vs. Tilt");
        energyTiltPlot.setDoLegend(false);
        energyTiltPlot.setDoClear(false);
        energyTiltPlot.setDoDrawLines(new DataTag[]{energyTilt.getTag()}, false);
        add(energyTiltPlot);
        getController().getResetAveragesButton().setPostAction(new IAction() {
            public void actionPerformed() {
                energyTiltPlot.getPlot().clear(false);
                stressStrainPlot.getPlot().clear(false);
            }
        });

        DisplayPlot plot = new DisplayPlot();
        historyKE.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{meterKE.getTag()}, "KE");
        historyPE.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{meterPE.getTag()}, "PE");
        historyE.setDataSink(plot.getDataSet().makeDataSink());
        plot.setLegend(new DataTag[]{meterEnergy.getTag()}, "E");
        plot.setLabel("Energy");
        
        corrugationButton = new DeviceButton(sim.getController());
        IAction corrugationAction = new IAction() {
            public void actionPerformed() {
                sinusoidalEnabled = !sinusoidalEnabled;
                sim.p1SurfaceBond.setB(sinusoidalEnabled ? sim.sinusoidalB : 0);
                ((P2LennardJones)sim.p2SulfurSurfaceLJ.getWrappedPotential()).setEpsilon(sinusoidalEnabled ? 0 : 
                    ((P2LennardJones)sim.p2CH2Surface.getWrappedPotential()).getEpsilon());
                sim.p2SurfaceBond.setSpringConstant(sinusoidalEnabled ? 0 : sim.harmonicStrength);
                sim.integrator.setSulfurType(sinusoidalEnabled ? sim.species.getSulfurType() : null);

                sim.wallPotential.setWallPosition(20);
                wallPositionSlider.doUpdate();
                sim.config.initializeCoordinates(sim.box);
                sim.findTetherBonds();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {}
                getDisplayBox(sim.box).repaint();
                resetAction.actionPerformed();
                corrugationButton.setLabel(sinusoidalEnabled ? "Use harmonic corrugation" : "Use sinusoidal corrugation");
            }
            boolean sinusoidalEnabled = false;
        };
        corrugationButton.setAction(corrugationAction);
        stateTab.add(corrugationButton.graphic(), SimulationPanel.getVertGBC());
        corrugationButton.setLabel("Use sinusoidal corrugation");

        add(plot);
        
        moveWallAction = new ActionMoveWall(wallPositionSlider, sim, moveWallToggle);
        IntegratorListenerAction moveWallActionListener = new IntegratorListenerAction(moveWallAction);
        sim.integrator.getEventManager().addListener(moveWallActionListener);
        moveWallActionListener.setInterval(200);
    }
    
    public static void main(String[] args) {
        Sam sim = new Sam();
        SamGraphic graphic = new SamGraphic(sim);
        graphic.makeAndDisplayFrame();
    }

    public static class DataPipeStressStrain extends DataProcessor {
        protected final Sam sim;
        protected DataDoubleArray data;
        protected DataInfoDoubleArray dataInfo;

        public DataPipeStressStrain(Sam sim) {
            super();
            this.sim = sim;
        }

        protected IData processData(IData inputData) {
            double[] xy = data.getData();
            double surfacePosition = sim.box.getMoleculeList(sim.speciesSurface).getMolecule(0).getChildList().getAtom(0).getPosition().getX(1);
            xy[0] = sim.wallPotential.wallPosition-surfacePosition;
            xy[1] = inputData.getValue(0);
            return data;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            data = new DataDoubleArray(2);
            dataInfo = new DataInfoDoubleArray("Strain and Stress", Null.DIMENSION, new int[]{2});
            return dataInfo;
        }
    }

    public class ActionMoveWall implements IAction {
        private final DeviceSlider wallPositionSlider;
        private final Sam sim;
        private final IAction moveWallToggle;
        public boolean enabled;

        protected ActionMoveWall(DeviceSlider wallPositionSlider, Sam sim,
                IAction moveWallToggle) {
            this.wallPositionSlider = wallPositionSlider;
            this.sim = sim;
            this.moveWallToggle = moveWallToggle;
            enabled = false;
        }

        public void actionPerformed() {
            if (!enabled) {
                return;
            }
            double position = sim.wallPotential.getWallPosition()-0.1;
            sim.wallPotential.setWallPosition(position);
            wallPositionSlider.doUpdate();
            double surfacePosition = sim.box.getMoleculeList(sim.speciesSurface).getMolecule(0).getChildList().getAtom(0).getPosition().getX(1);
            if (position-surfacePosition < 15.01) {
                // turn ourselves off
                moveWallToggle.actionPerformed();
            }
        }
    }

    public static class ModifierWallPosition implements Modifier {
        private final Sam sim;
        protected double surfacePosition;

        public ModifierWallPosition(Sam sim) {
            this.sim = sim;
            update();
        }

        public void update() {
            surfacePosition = sim.box.getMoleculeList(sim.speciesSurface).getMolecule(0).getChildList().getAtom(0).getPosition().getX(1);
        }

        public void setValue(double newValue) {
            double maxAtomPos = -100;
            IAtomList leafList = sim.box.getLeafList();
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtom atom = leafList.getAtom(i);
                double atomPos = atom.getPosition().getX(1);
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
