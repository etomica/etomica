/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.space.Vector;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataProcessorFunction;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.geometry.Plane;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.CriterionPositionWall;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.Pixel;
import etomica.units.dimensions.Quantity;
import etomica.math.function.Function;
import etomica.data.history.HistoryCollapsingDiscard;

/**
 * Colloid module app.  Design by Alberto Striolo.
 * 
 * @author Andrew Schultz
 */
public class ColloidGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Colloid";
    private final static int REPAINT_INTERVAL = 10;
    protected ColloidSim sim;
    protected final MeterProfileByVolume densityProfileMeter, colloidDensityProfileMeter;
    protected final MeterEnd2End meterE2E;
    
    public ColloidGraphic(final ColloidSim simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

    	getController().getSimRestart().setConfiguration(sim.configuration);
    	getController().getSimRestart().setIgnoreOverlap(true);
    	getController().getReinitButton().setPostAction(getPaintAction(sim.box));

    	((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), sim.p2mm.getCoreDiameter());
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesColloid.getLeafType(), 7.5);

        final Plane planeBottom = new Plane(space,0,1,0,sim.box.getBoundary().getBoxSize().getX(1)*0.5);
        final Plane planeTop = new Plane(space,0,1,0,-sim.box.getBoundary().getBoxSize().getX(1)*0.5);
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).addPlane(planeBottom);
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).addPlane(planeTop);
    	
    	getDisplayBox(sim.box).setPixelUnit(new Pixel(2));

        sim.activityIntegrate.setSleepPeriod(0);
       
        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");
        
        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));

        DeviceSelector graftSelector = new DeviceSelector(sim.getController());
        graftSelector.addOption("1", new GraftAction(1));
        graftSelector.addOption("2", new GraftAction(2));
        graftSelector.addOption("4", new GraftAction(4));
        graftSelector.addOption("6", new GraftAction(6));
        graftSelector.addOption("8", new GraftAction(8));
        graftSelector.addOption("12", new GraftAction(12));
        graftSelector.setSelected(5);
        graftSelector.setLabel("# of grafted chains");
        add(graftSelector);
        
        DeviceSlider chainLengthSlider = new DeviceSlider(sim.getController());
        chainLengthSlider.setModifier(new Modifier() {
            public void setValue(double newValue) {
                if (newValue > 2*(sim.box.getBoundary().getBoxSize().getX(1) - sim.getColloidSigma())) {
                    throw new IllegalArgumentException("too large for box size");
                }
                sim.setChainLength((int)newValue);
                meterE2E.setChainLength((int)newValue);
                int n = sim.box.getLeafList().getAtomCount();
                sim.integrator.setThermostatInterval((1000+(n-1))/n);
            }
            public double getValue() {
                return sim.getChainLength();
            }
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }
            public String getLabel() {
                return "chain length";
            }
        });
        chainLengthSlider.setLabel("Chain length");
        chainLengthSlider.setShowBorder(true);
        chainLengthSlider.setMaximum(50);
        chainLengthSlider.setNMajor(5);
        chainLengthSlider.setPostAction(getPaintAction(sim.box));
        chainLengthSlider.setShowValues(true);
        add(chainLengthSlider);
        
        MeterTemperature meterTemperature = new MeterTemperature(sim.box, 3);
        DisplayTextBox displayTemperature = new DisplayTextBox();
        DataPumpListener tempPump = new DataPumpListener(meterTemperature, displayTemperature);
        sim.getIntegrator().getEventManager().addListener(tempPump);
        add(displayTemperature);
        
        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        thermoSlider.setMaximum(5);
        thermoSlider.setShowValues(true);
        thermoSlider.setPrecision(1);
        add(thermoSlider);
        
        DeviceSlider boxSizeSlider = new DeviceSlider(sim.getController());
        boxSizeSlider.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                if (newValue < sim.getColloidSigma() + 1 || newValue < 12) {
                    throw new IllegalArgumentException("box size too small");
                }
                Vector v = space.makeVector();
                v.E(sim.box.getBoundary().getBoxSize());
                v.setX(1, newValue);
                sim.box.setNMolecules(sim.species, 0);
                sim.box.setNMolecules(sim.speciesColloid, 0);
                sim.box.getBoundary().setBoxSize(v);
                sim.configuration.initializeCoordinates(sim.box);
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {}
                planeBottom.setDistanceToOrigin(0.5*newValue);
                planeTop.setDistanceToOrigin(-0.5*newValue);

                densityProfileMeter.reset();
                colloidDensityProfileMeter.reset();

                getDisplayBox(sim.box).repaint();
            }
            
            public double getValue() {
                return sim.box.getBoundary().getBoxSize().getX(1);
            }
            
            public String getLabel() {
                return "Box size";
            }
            
            public Dimension getDimension() {
                return Length.DIMENSION;
            }
        });
        boxSizeSlider.setMaximum(150);
        boxSizeSlider.setShowValues(true);
        boxSizeSlider.setLabel("Box Size");
        boxSizeSlider.setShowBorder(true);
        add(boxSizeSlider);

//        JTabbedPane potentialTabs = new JTabbedPane();
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
//        JPanel monomerPanel = new JPanel(new GridLayout(0,2));
//        DeviceBox monomerRangeBox = new DeviceBox();
//        monomerRangeBox.setController(sim.getController());
//        monomerRangeBox.setLabel("Wall range");
//        monomerRangeBox.setModifier(new WallRangeModifier(sim.p1WallMonomer, sim.criterionWallMonomer));
//        monomerPanel.add(monomerRangeBox.graphic());
        
//        DeviceBox monomerEpsilonBox = new DeviceBox();
//        monomerEpsilonBox.setController(sim.getController());
//        monomerEpsilonBox.setLabel("W-W epsilon");
//        monomerEpsilonBox.setModifier(new Modifier() {
//            
//            public void setValue(double newValue) {
//                sim.epsWallWall = newValue;
//                double epsMW = Math.sqrt(newValue*sim.p2mm.getEpsilon());
//                sim.p1WallMonomer.setEpsilon(epsMW);
//            }
//            
//            public double getValue() { return sim.epsWallWall; }
//            public String getLabel() { return "W-W epsilon"; }
//            public Dimension getDimension() { return Energy.DIMENSION; }
//        });
//
//        monomerPanel.add(monomerEpsilonBox.graphic());

        DeviceBox mmEpsilonBox = new DeviceBox();
        mmEpsilonBox.setController(sim.getController());
        mmEpsilonBox.setLabel("monomer epsilon");
        mmEpsilonBox.setModifier(new ModifierGeneral(sim.p2mm, "epsilon"));
        mmEpsilonBox.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                sim.p2mm.setEpsilon(newValue);
                double epsMW = Math.sqrt(newValue*sim.epsWallWall);
                sim.p1WallMonomer.setEpsilon(epsMW);
            }
            
            public double getValue() { return sim.p2mm.getEpsilon(); }
            public String getLabel() { return "monomer epsilon"; }
            public Dimension getDimension() { return Energy.DIMENSION; }
        });
        gbc.gridx = 0;
        gbc.gridy = 0;
        potentialPanel.add(mmEpsilonBox.graphic(), gbc);

//        potentialPanel.add(monomerPanel, "monomer");

//        JPanel colloidPanel = new JPanel(new GridLayout(0,2));
//        DeviceBox colloidRangeBox = new DeviceBox();
//        colloidRangeBox.setController(sim.getController());
//        colloidRangeBox.setLabel("Wall range");
//        colloidRangeBox.setModifier(new WallRangeModifier(sim.p1WallColloid, null));
//        colloidPanel.add(colloidRangeBox.graphic());

        DeviceBox colloidSigmaBox = new DeviceBox();
        colloidSigmaBox.setController(sim.getController());
        colloidSigmaBox.setLabel("colloid sigma");
        colloidSigmaBox.setModifier(new ModifierGeneral(sim, "colloidSigma"));
        colloidSigmaBox.setPostAction(new IAction() {
            public void actionPerformed() {
                ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesColloid.getLeafType(), sim.getColloidSigma());
                getPaintAction(sim.box).actionPerformed();
            }
        });
        gbc.gridx = 1;
        potentialPanel.add(colloidSigmaBox.graphic(), gbc);

        DeviceBox colloidEpsilonBox = new DeviceBox();
        colloidEpsilonBox.setController(sim.getController());
        colloidEpsilonBox.setLabel("colloid-wall epsilon");
        colloidEpsilonBox.setModifier(new ModifierGeneral(sim.p1WallColloid, "epsilon"));
        gbc.gridx = 0;
        gbc.gridy = 1;
        gbc.gridwidth = 2;
        potentialPanel.add(colloidEpsilonBox.graphic(), gbc);


//        potentialPanel.add(colloidPanel, "colloid");
        
        getPanel().controlPanel.add(potentialPanel,SimulationPanel.getVertGBC());


        densityProfileMeter = new MeterProfileByVolume(space);
        densityProfileMeter.setProfileDim(1);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(100);
        DataPumpListener profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, 10);
        sim.integrator.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        colloidDensityProfileMeter = new MeterProfileByVolume(space);
        colloidDensityProfileMeter.setProfileDim(1);
        colloidDensityProfileMeter.setBox(sim.box);
        meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.speciesColloid);
        colloidDensityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed colloidDensityProfileAvg = new AccumulatorAverageFixed(10);
        colloidDensityProfileAvg.setPushInterval(100);
        profilePump = new DataPumpListener(colloidDensityProfileMeter, colloidDensityProfileAvg, 10);
        sim.integrator.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag()}, "monomer");
        profilePlot.getPlot().setTitle("Monomer");
        profilePlot.setXLabel("y");
        profilePlot.setDoLegend(false);
        profilePlot.setLabel("Monomer Density");

        DisplayPlot colloidProfilePlot = new DisplayPlot();
        colloidDensityProfileAvg.addDataSink(colloidProfilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{colloidDensityProfileAvg.AVERAGE});
        colloidProfilePlot.setLegend(new DataTag[]{colloidDensityProfileAvg.getTag()}, "colloid");
        colloidProfilePlot.getPlot().setTitle("Colloid");
        colloidProfilePlot.setXLabel("y");
        colloidProfilePlot.setDoLegend(false);
        colloidProfilePlot.setLabel("Colloid Density");

        JPanel plotPanel = new JPanel(new GridLayout(2,1));
        plotPanel.add(profilePlot.graphic());
        plotPanel.add(colloidProfilePlot.graphic());

        java.awt.Dimension d = profilePlot.getPlot().getPreferredSize();
        d.width -= 100;
        d.height = 210;
        profilePlot.getPlot().setSize(d);
        colloidProfilePlot.getPlot().setSize(d);

        addAsTab(plotPanel, "Density Profiles", true);

        meterE2E = new MeterEnd2End(space);
        meterE2E.setBox(sim.box);
        meterE2E.setChainLength(sim.getChainLength());
        AccumulatorAverageCollapsing avgE2E = new AccumulatorAverageCollapsing();
        avgE2E.setPushInterval(10);
        DataPumpListener pumpE2E = new DataPumpListener(meterE2E, avgE2E, 10);
        sim.integrator.getEventManager().addListener(pumpE2E);
        dataStreamPumps.add(pumpE2E);
//        DisplayTextBoxesCAE displayE2E = new DisplayTextBoxesCAE();
//        displayE2E.setAccumulator(avgE2E);
//        displayE2E.setLabel("end to end distance");
//        add(displayE2E);
        
        DataProcessorFunction sqrtE2E = new DataProcessorFunction(Function.Sqrt.INSTANCE);
        avgE2E.addDataSink(sqrtE2E, new StatType[]{avgE2E.AVERAGE});
        AccumulatorHistory historyE2E = new AccumulatorHistory(new HistoryCollapsingDiscard());
        sqrtE2E.setDataSink(historyE2E);
        DisplayPlot runningAvgE2E = new DisplayPlot();
        historyE2E.setDataSink(runningAvgE2E.getDataSet().makeDataSink());
        runningAvgE2E.setLabel("End-to-End Distance");
        add(runningAvgE2E);
    }

    public static class WallRangeModifier implements Modifier {
        public WallRangeModifier(P1Wall p1, CriterionPositionWall criterion) {
            this.p1 = p1;
            this.criterion = criterion;
        }

        public Dimension getDimension() {
            return Length.DIMENSION;
        }

        public String getLabel() {
            return "Wall range";
        }

        public double getValue() {
            return p1.getRange();
        }

        public void setValue(double newValue) {
            if (newValue == getValue()) return;
            p1.setRange(newValue);
            if (criterion != null) {
                criterion.setNeighborRange(1.5*newValue);
            }
        }
        
        protected final P1Wall p1;
        protected final CriterionPositionWall criterion;
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        ColloidGraphic swmdGraphic = new ColloidGraphic(new ColloidSim(space), space);
		SimulationGraphic.makeAndDisplayFrame
		        (swmdGraphic.getPanel(), APP_NAME);
    }
    
    public class GraftAction implements IAction {
        public GraftAction(int numGraft) {
            nGraft = numGraft;
        }
        public void actionPerformed() {
            sim.setNumGraft(nGraft);
            int n = sim.box.getLeafList().getAtomCount();
            sim.integrator.setThermostatInterval((1000+(n-1))/n);
            getDisplayBox(sim.box).repaint();
        }
        protected final int nGraft;
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            int dim = 3;
            Space sp = Space.getInstance(dim);
            ColloidGraphic swmdGraphic = new ColloidGraphic(new ColloidSim(sp), sp);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
