/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Quantity;
import etomica.util.Constants.CompassDirection;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

/**
 * Module for chain reaction (polymerization) using ChainEquilibriumSim as the
 * simulation class.  Original module by William Scharmach and Matt Moynihan.
 * Later revamped based on module redesign by William M. Chirdon.
 * 
 * @author William Scharmach
 * @author Matt Moynihan
 * @author Andrew Schultz
 */
public class ChainEquilibriumGraphic extends SimulationGraphic {

	private static final String APP_NAME = "Stepwise Growth Polymerization";
	private static final int REPAINT_INTERVAL = 2;

    protected ChainEquilibriumSim sim;

    public ChainEquilibriumGraphic(ChainEquilibriumSim simulation, Space _space) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        getPanel().toolbar.addAuthor("Dr. William Chirdon");
        this.sim = simulation;
        
        int dataInterval = (int) (.04 / sim.integratorHard.getTimeStep());
        
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        
        getDisplayBox(sim.box).setPixelUnit(new Pixel(7));

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******	
        int eMin = 0, eMax = 40;

        // **** Stuff that Modifies the Simulation

        final IAction resetAction = getController().getSimRestart().getDataResetAction();
        
        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.controller1, sim.activityIntegrate);

        // Sliders on Well depth page
        final DeviceSlider ABSlider = sliders(eMin, eMax, "Diol-Carboxylic Acid", sim.ABbonded);
//        final DeviceSlider ACSlider = sliders(eMin, eMax, "Diol-Crosslinker", sim.ACbonded);
        ABSlider.setPostAction(resetAction);
//        ACSlider.setPostAction(resetAction);
        
        DeviceBox solventThermoFrac = new DeviceBox();
        solventThermoFrac.setController(sim.getController());
        solventThermoFrac.setModifier(new ModifierGeneral(sim.ABbonded, "solventThermoFrac"));
        solventThermoFrac.setLabel("fraction heat transfer to solvent");
        DisplayTextBox tBox = new DisplayTextBox();

        DisplayTimer displayTimer = new DisplayTimer(sim.integratorHard);
        add(displayTimer);
        
        DataSourceCountTime timer = new DataSourceCountTime(sim.integratorHard);

        DataFork tFork = new DataFork();
        final DataPump tPump = new DataPump (sim.thermometer, tFork);
        tFork.addDataSink(tBox);
        add(tBox);
        dataStreamPumps.add(tPump);
        AccumulatorHistory tHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        tHistory.setTimeDataSource(timer);
        tFork.addDataSink(tHistory);
        DisplayPlot tPlot = new DisplayPlot();
        tHistory.addDataSink(tPlot.getDataSet().makeDataSink());
        tPlot.setUnit(Kelvin.UNIT);
        tPlot.setLabel("Temperature");
        tPlot.getPlot().setYLabel("Temperature (K)");
        tPlot.setDoLegend(false);
        add(tPlot);
        IntegratorListenerAction tPumpListener = new IntegratorListenerAction(tPump);
        sim.integratorHard.getEventManager().addListener(tPumpListener);
        tPumpListener.setInterval(dataInterval);

        final MeterChainLength molecularCount = new MeterChainLength(sim.agentManager);
        molecularCount.setBox(sim.box);
        AccumulatorAverage accumulator = new AccumulatorAverageFixed(10);
        accumulator.setPushInterval(10);
        DataFork mwFork = new DataFork();
        final DataPump mwPump = new DataPump(molecularCount,mwFork);
        mwFork.addDataSink(accumulator);
        dataStreamPumps.add(mwPump);
        IntegratorListenerAction mwPumpListener = new IntegratorListenerAction(mwPump);
        sim.integratorHard.getEventManager().addListener(mwPumpListener);
        mwPumpListener.setInterval(dataInterval);
        
        MolecularWeightAvg molecularWeightAvg = new MolecularWeightAvg();
        mwFork.addDataSink(molecularWeightAvg);
        DataFork mwAvgFork = new DataFork();
        molecularWeightAvg.setDataSink(mwAvgFork);
        AccumulatorAverageCollapsing mwAvg = new AccumulatorAverageCollapsing();
        mwAvg.setPushInterval(1);
        mwAvgFork.addDataSink(mwAvg);
        final AccumulatorHistory mwHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        mwHistory.setTimeDataSource(timer);
        mwAvgFork.addDataSink(mwHistory);

        MolecularWeightAvg2 molecularWeightAvg2 = new MolecularWeightAvg2();
        mwFork.addDataSink(molecularWeightAvg2);
        DataFork mwAvg2Fork = new DataFork();
        molecularWeightAvg2.setDataSink(mwAvg2Fork);
        AccumulatorAverageCollapsing mwAvg2 = new AccumulatorAverageCollapsing();
        mwAvg2.setPushInterval(1);
        mwAvg2Fork.addDataSink(mwAvg2);
        final AccumulatorHistory mw2History = new AccumulatorHistory(new HistoryCollapsingAverage());
        mw2History.setTimeDataSource(timer);
        mwAvg2Fork.addDataSink(mw2History);

        MeterConversion reactionConversionDiol = new MeterConversion(sim.box, sim.agentManager);
        reactionConversionDiol.setSpecies(new ISpecies[]{sim.speciesA});
        final HistoryCollapsingAverage conversionHistoryDiol = new HistoryCollapsingAverage();
        AccumulatorHistory conversionHistoryAccDiol = new AccumulatorHistory(conversionHistoryDiol);
        conversionHistoryAccDiol.setTimeDataSource(timer);
        final DataPump conversionPumpDiol = new DataPump(reactionConversionDiol, conversionHistoryAccDiol);
        IntegratorListenerAction conversionPumpDiolListener = new IntegratorListenerAction(conversionPumpDiol);
        sim.integratorHard.getEventManager().addListener(conversionPumpDiolListener);
        conversionPumpDiolListener.setInterval(dataInterval);
        dataStreamPumps.add(conversionPumpDiol);

        MeterConversion reactionConversionAcid = new MeterConversion(sim.box, sim.agentManager);
        reactionConversionAcid.setSpecies(new ISpecies[]{sim.speciesB});
        final HistoryCollapsingAverage conversionHistoryAcid = new HistoryCollapsingAverage();
        AccumulatorHistory conversionHistoryAccAcid = new AccumulatorHistory(conversionHistoryAcid);
        conversionHistoryAccAcid.setTimeDataSource(timer);
        final DataPump conversionPumpAcid = new DataPump(reactionConversionAcid, conversionHistoryAccAcid);
        IntegratorListenerAction conversionPumpAcidListener = new IntegratorListenerAction(conversionPumpAcid);
        sim.integratorHard.getEventManager().addListener(conversionPumpAcidListener);
        conversionPumpAcidListener.setInterval(dataInterval);
        dataStreamPumps.add(conversionPumpAcid);

        final IAction resetData = new IAction() {
            public void actionPerformed() {
                sim.integratorHard.resetStepCount();
                molecularCount.reset();
                conversionPumpDiol.actionPerformed();
                conversionPumpAcid.actionPerformed();
                mwPump.actionPerformed();
                tPump.actionPerformed();
            }
        };

        getController().getResetAveragesButton().setLabel("Reset");
        getController().getResetAveragesButton().setPostAction(resetData);

        DisplayPlot compositionPlot = new DisplayPlot();
        accumulator.addDataSink(compositionPlot.getDataSet().makeDataSink(),new AccumulatorAverage.StatType[]{accumulator.AVERAGE});
        compositionPlot.setDoLegend(false);

        DisplayPlot mwPlot = new DisplayPlot();
        mwPlot.setLabel("Molecular Weight");
        mwHistory.addDataSink(mwPlot.getDataSet().makeDataSink());
        mwPlot.setLegend(new DataTag[]{mwHistory.getTag()}, "Number Avg");
        mw2History.addDataSink(mwPlot.getDataSet().makeDataSink());
        mwPlot.setLegend(new DataTag[]{mw2History.getTag()}, "Weight Avg");

        DisplayPlot conversionPlot = new DisplayPlot();
        conversionHistoryAccDiol.addDataSink(conversionPlot.getDataSet().makeDataSink());
        conversionPlot.setLegend(new DataTag[]{reactionConversionDiol.getTag()}, "diol conversion");
        conversionHistoryAccAcid.addDataSink(conversionPlot.getDataSet().makeDataSink());
        conversionPlot.setLegend(new DataTag[]{reactionConversionAcid.getTag()}, "acid conversion");

        DisplayTextBoxesCAE mwBox = new DisplayTextBoxesCAE();
        mwBox.setAccumulator(mwAvg);
        mwBox.setLabel("Number Avg Molecular Weight");
        add(mwBox);

        DisplayTextBoxesCAE mw2Box = new DisplayTextBoxesCAE();
        mw2Box.setAccumulator(mwAvg2);
        mw2Box.setLabel("Weight Avg Molecular Weight");
        add(mw2Box);

        ((SimulationRestart)getController().getReinitButton().getAction()).setConfiguration(sim.config);
        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                sim.resetBonds();
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        });

        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.controller1, sim.integratorHard);
        temperatureSelect.setUnit(Kelvin.UNIT);
        temperatureSelect.setMaximum(1200);
        temperatureSelect.setIsothermal();
        temperatureSelect.setSliderPostAction(resetAction);
        temperatureSelect.setRadioGroupPostAction(resetAction);
        
        ColorSchemeStepWise colorScheme = new ColorSchemeStepWise(sim, sim.agentManager);
        colorScheme.setColor(sim.speciesA.getLeafType(), 1, new Color(255, 120, 120));
        colorScheme.setColor(sim.speciesB.getLeafType(), 1, new Color(120, 120, 255));
        colorScheme.setColor(sim.speciesA.getLeafType(), 2, new Color(200, 0, 0));
        colorScheme.setColor(sim.speciesB.getLeafType(), 2, new Color(0, 0, 200));
        colorScheme.setColor(sim.speciesB.getLeafType(), 3, Color.GREEN);
        getDisplayBox(sim.box).setColorScheme(colorScheme);

        IAction reconfig = new IAction() {
            public void actionPerformed() {
                getController().getSimRestart().actionPerformed();
                sim.resetBonds();
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        };
        DeviceBox nMonoOlBox = new DeviceBox();
        nMonoOlBox.setInteger(true);
        nMonoOlBox.setController(sim.getController());
        nMonoOlBox.setLabel("Mono-ol (light red)");
        nMonoOlBox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNMonoOl();}
            public void setValue(double newValue) {sim.setNMonoOl((int)newValue);}
        });
        nMonoOlBox.setPostAction(reconfig);
        DeviceBox nMonoAcidBox = new DeviceBox();
        nMonoAcidBox.setInteger(true);
        nMonoAcidBox.setController(sim.getController());
        nMonoAcidBox.setLabel("Mono-acid (light blue)");
        nMonoAcidBox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNMonoAcid();}
            public void setValue(double newValue) {sim.setNMonoAcid((int)newValue);}
        });
        nMonoAcidBox.setPostAction(reconfig);
        DeviceBox nDiolBox = new DeviceBox();
        nDiolBox.setInteger(true);
        nDiolBox.setController(sim.getController());
        nDiolBox.setLabel("Di-ol (red)");
        nDiolBox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNDiol();}
            public void setValue(double newValue) {sim.setNDiol((int)newValue);}
        });
        nDiolBox.setPostAction(reconfig);
        DeviceBox nDiAcidBox = new DeviceBox();
        nDiAcidBox.setInteger(true);
        nDiAcidBox.setController(sim.getController());
        nDiAcidBox.setLabel("Di-acid (blue)");
        nDiAcidBox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNDiAcid();}
            public void setValue(double newValue) {sim.setNDiAcid((int)newValue);}
        });
        nDiAcidBox.setPostAction(reconfig);
        DeviceBox nCrossLinkerBox = new DeviceBox();
        nCrossLinkerBox.setInteger(true);
        nCrossLinkerBox.setController(sim.getController());
        nCrossLinkerBox.setLabel("Crosslinker (green)");
        nCrossLinkerBox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNCrossLinkersAcid();}
            public void setValue(double newValue) {sim.setNCrossLinkersAcid((int)newValue);}
        });
        nCrossLinkerBox.setPostAction(reconfig);

        tBox.setUnit(Kelvin.UNIT);
        tBox.setLabel("Measured Temperature");
        tBox.setLabelPosition(CompassDirection.NORTH);

        compositionPlot.setLabel("Composition");
        conversionPlot.setLabel("Conversion");

        add(compositionPlot);
        add(mwPlot);
        JPanel conversionPanel = new JPanel(new GridBagLayout());
        conversionPanel.add(conversionPlot.graphic(), vertGBC);
        
        DeviceBox conversionHistoryLength = new DeviceBox();
        conversionHistoryLength.setInteger(true);
        conversionHistoryLength.setController(sim.getController());
        conversionHistoryLength.setModifier(new Modifier() {

            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "history length";
            }

            public double getValue() {
                return conversionHistoryDiol.getHistoryLength();
            }

            public void setValue(double newValue) {
                conversionHistoryDiol.setHistoryLength((int)newValue);
                conversionHistoryAcid.setHistoryLength((int)newValue);
            }
        });
        conversionPanel.add(conversionHistoryLength.graphic(),vertGBC);

        final DeviceButton atomFilterButton;
        if (space.D() == 3) {
            final AtomTestChainLength atomFilter = new AtomTestChainLength(sim.agentManager);
            atomFilter.setBox(sim.box);
            atomFilterButton = new DeviceButton(sim.getController());
            atomFilterButton.setAction(new IAction() {
                public void actionPerformed() {
                    DisplayBox displayBox = getDisplayBox(sim.box);
                    if (displayBox.setAtomTestDoDisplay() == null) {
                        displayBox.setAtomTestDoDisplay(atomFilter);
                        atomFilterButton.setLabel("Show all");
                    }
                    else {
                        displayBox.setAtomTestDoDisplay(null);
                        atomFilterButton.setLabel("Show only longest chain");
                    }
                    displayBox.repaint();
                }
            });
            atomFilterButton.setLabel("Show only longest chain");
        }
        else {
            atomFilterButton = null;
        }

        getPanel().tabbedPane.add("Conversion" , conversionPanel);

        JPanel speciesEditors = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel epsilonSliders = new JPanel(new java.awt.GridBagLayout());
        JPanel controls = new JPanel(new java.awt.GridBagLayout());

        speciesEditors.add(nMonoOlBox.graphic());
        speciesEditors.add(nMonoAcidBox.graphic());
        speciesEditors.add(nDiolBox.graphic());
        speciesEditors.add(nDiAcidBox.graphic());
        speciesEditors.add(nCrossLinkerBox.graphic());

        epsilonSliders.add(ABSlider.graphic(null), vertGBC);
//        epsilonSliders.add(ACSlider.graphic(null), vertGBC);
        epsilonSliders.add(solventThermoFrac.graphic(), vertGBC);

        final JTabbedPane sliderPanel = new JTabbedPane();
        //panel for all the controls
        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);
        sliderPanel.add(controls, "Controls");
        sliderPanel.add(epsilonSliders, "Reaction Energy (kJ/mol)");
        sliderPanel.add(speciesEditors, "Number of Molecules");
        controls.add(delaySlider.graphic(), vertGBC);
        if (space.D() == 3) {
            controls.add(atomFilterButton.graphic(), vertGBC);
        }

        //set the number of significant figures displayed on the table.
        javax.swing.table.DefaultTableCellRenderer numberRenderer = new javax.swing.table.DefaultTableCellRenderer() {
            java.text.NumberFormat formatter;
            {
                formatter = java.text.NumberFormat.getInstance();
                formatter.setMaximumFractionDigits(6);
            }

            public void setValue(Object value) {
                setText((value == null) ? "" : formatter.format(value));
            }
        };

        numberRenderer.setHorizontalAlignment(SwingConstants.RIGHT);
    }

    public static void main(String[] args) {
        int D = 2;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-dim") && i + 1 < args.length) {
                i++;
                D = Integer.parseInt(args[i]);
            }
        }
        ChainEquilibriumSim sim = new ChainEquilibriumSim(Space.getInstance(D));
        ChainEquilibriumGraphic graphic = new ChainEquilibriumGraphic(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }

    public DeviceSlider sliders(int eMin, int eMax, String s, P2SquareWellBonded p){

        DeviceSlider AASlider = new DeviceSlider(sim.getController(), new ModifierGeneral(p, "epsilon"));
        AASlider.setUnit(new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT));
        AASlider.doUpdate();
        AASlider.setShowBorder(true);
        AASlider.setLabel(s);
        AASlider.setMinimum(eMin);
        AASlider.setMaximum(eMax);
        AASlider.setNMajor(4);
//        AASlider.getSlider().setSnapToTicks(true);

        return AASlider;
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
			getRootPane().putClientProperty("defeatSystemEventQueueCheck", Boolean.TRUE);
            int D = 2;
            String dimStr = getParameter("dim");
            if (dimStr != null) {
                D = Integer.parseInt(dimStr);
            }
	        ChainEquilibriumSim sim = new ChainEquilibriumSim(Space.getInstance(D));
			getContentPane().add(new ChainEquilibriumGraphic(sim, sim.getSpace()).getPanel());
        }
    }
}
