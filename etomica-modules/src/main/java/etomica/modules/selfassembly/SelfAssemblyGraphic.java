/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.selfassembly;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.Atom;
import etomica.atom.AtomTest;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2SquareWell;
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
 * @author Andrew Schultz
 */
public class SelfAssemblyGraphic extends SimulationGraphic {

	private static final String APP_NAME = "Self-assembly";
	private static final int REPAINT_INTERVAL = 2;
    boolean displayA = true;
    boolean displayB1 = true;
    boolean displayB2 = true;


    protected SelfAssemblySim sim;

    public SelfAssemblyGraphic(SelfAssemblySim simulation, Space _space) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        
        int dataInterval = (int) (.04 / sim.integratorHard.getTimeStep());
        
        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        
        getDisplayBox(sim.box).setPixelUnit(new Pixel(7));

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******	
        int eMin = 0, eMax = 40;

        // **** Stuff that Modifies the Simulation

        final IAction resetAction = getController().getSimRestart().getDataResetAction();
        
        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.controller1);

        // Sliders on Well depth page
//        final DeviceSlider ABSlider = sliders(eMin, eMax, "Diol-Carboxylic Acid", sim.ABbonded);
//        final DeviceSlider ACSlider = sliders(eMin, eMax, "Diol-Crosslinker", sim.ACbonded);
//        ABSlider.setPostAction(resetAction);
//        ACSlider.setPostAction(resetAction);
        
//        DeviceBox solventThermoFrac = new DeviceBox();
//        solventThermoFrac.setController(sim.getController());
//        solventThermoFrac.setModifier(new ModifierGeneral(sim.ABbonded, "solventThermoFrac"));
//        solventThermoFrac.setLabel("fraction heat transfer to solvent");
//        DisplayTextBox tBox = new DisplayTextBox();

        DisplayTimer displayTimer = new DisplayTimer(sim.integratorHard);
        add(displayTimer);
        
        DataSourceCountTime timer = new DataSourceCountTime(sim.integratorHard);

//        DataFork tFork = new DataFork();
//        final DataPump tPump = new DataPump (sim.thermometer, tFork);
//        tFork.addDataSink(tBox);
//        add(tBox);
//        dataStreamPumps.add(tPump);
//        AccumulatorHistory tHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        tHistory.setTimeDataSource(timer);
//        tFork.addDataSink(tHistory);
//        DisplayPlotXChart tPlot = new DisplayPlotXChart();
//        tHistory.addDataSink(tPlot.getDataSet().makeDataSink());
//        tPlot.setUnit(Kelvin.UNIT);
//        tPlot.setLabel("Temperature");
//        tPlot.getPlot().setYLabel("Temperature (K)");
//        tPlot.setDoLegend(false);
//        add(tPlot);
//        IntegratorListenerAction tPumpListener = new IntegratorListenerAction(tPump);
//        sim.integratorHard.getEventManager().addListener(tPumpListener);
//        tPumpListener.setInterval(dataInterval);


        final IAction resetData = new IAction() {
            public void actionPerformed() {
                sim.integratorHard.resetStepCount();
//                molecularCount.reset();
//                conversionPumpDiol.actionPerformed();
//                conversionPumpAcid.actionPerformed();
//                mwPump.actionPerformed();
//                tPump.actionPerformed();
            }
        };

        getController().getResetAveragesButton().setLabel("Reset");
        getController().getResetAveragesButton().setPostAction(resetData);

        ((SimulationRestart)getController().getReinitButton().getAction()).setConfiguration(sim.config);
        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        });

        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.controller1, sim.integratorHard);
        temperatureSelect.setUnit(Kelvin.UNIT);
        temperatureSelect.setMaximum(2000);
        temperatureSelect.setIsothermal();
        temperatureSelect.setSliderPostAction(resetAction);
        temperatureSelect.setRadioGroupPostAction(resetAction);
        
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.getTypeA(),new Color(200, 0, 0));
        colorScheme.setColor(sim.getTypeB1(), new Color(0, 0, 200));
        colorScheme.setColor(sim.getTypeB2(), new Color(0, 200, 0));
        getDisplayBox(sim.box).setColorScheme(colorScheme);

        getController().getSimRestart().setIgnoreOverlap(true);
        IAction reconfig = new IAction() {
            public void actionPerformed() {
                getController().getSimRestart().actionPerformed();
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        };

// ******** Set number of molecules, atoms per copolymer segment  *****************
        DeviceBox ABox = new DeviceBox();
        ABox.setInteger(true);
        ABox.setController(sim.getController());
        ABox.setLabel("Solvent");
        ABox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "nA";}
            public double getValue() {return sim.getNA();}
            public void setValue(double newValue) {sim.setNA((int)newValue);}
        });
        ABox.setPostAction(reconfig);
        DeviceBox BBox = new DeviceBox();
        BBox.setInteger(true);
        BBox.setController(sim.getController());
        BBox.setLabel("Polymer");
        BBox.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNB();}
            public void setValue(double newValue) {sim.setNB((int)newValue);}
        });
        BBox.setPostAction(reconfig);

        DeviceBox nB1 = new DeviceBox();
        nB1.setInteger(true);
        nB1.setController(sim.getController());
        nB1.setLabel("B1 segments (blue)");
        nB1.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNB1();}
            public void setValue(double newValue) {sim.setNB1((int)newValue);}
        });
        nB1.setPostAction(reconfig);

        DeviceBox nB2 = new DeviceBox();
        nB2.setInteger(true);
        nB2.setController(sim.getController());
        nB2.setLabel("B2 segments (green)");
        nB2.setModifier(new Modifier() {
            public Dimension getDimension() {return Quantity.DIMENSION;}
            public String getLabel() {return "n";}
            public double getValue() {return sim.getNB2();}
            public void setValue(double newValue) {sim.setNB2((int)newValue);}
        });
        nB2.setPostAction(reconfig);

//        tBox.setUnit(Kelvin.UNIT);
//        tBox.setLabel("Measured Temperature");
//        tBox.setLabelPosition(CompassDirection.NORTH);

        // ******** Turn off/on display of some atoms (3D only) **********
        AtomTest displayTest = new AtomTest() {
            public boolean test(IAtom a) {
                AtomType atomType = a.getType();
                if(atomType == sim.typeA) return displayA;
                else if(atomType == sim.typeB1) return displayB1;
                else return displayB2;
            }
        };

        final DeviceButton atomFilterButtonA, atomFilterButtonB1, atomFilterButtonB2;
        if (space.D() == 3) {
            getDisplayBox(sim.box).setAtomTestDoDisplay(displayTest);
            atomFilterButtonA = new DeviceButton(sim.getController());
            atomFilterButtonA.setAction(new IAction() {
                public void actionPerformed() {
                    displayA = !displayA;
                    if(displayA) atomFilterButtonA.setLabel("Hide A");
                    else atomFilterButtonA.setLabel("Show A");
                    getDisplayBox(sim.box).repaint();
                }
            });
            atomFilterButtonA.setLabel("Hide A");

            atomFilterButtonB1 = new DeviceButton(sim.getController());
            atomFilterButtonB1.setAction(new IAction() {
                public void actionPerformed() {
                    displayB1 = !displayB1;
                    if(displayB1) atomFilterButtonB1.setLabel("Hide B1");
                    else atomFilterButtonB1.setLabel("Show B1");
                    getDisplayBox(sim.box).repaint();
                }
            });
            atomFilterButtonB1.setLabel("Hide B1");

            atomFilterButtonB2 = new DeviceButton(sim.getController());
            atomFilterButtonB2.setAction(new IAction() {
                public void actionPerformed() {
                    displayB2 = !displayB2;
                    if(displayB2) atomFilterButtonB2.setLabel("Hide B2");
                    else atomFilterButtonB2.setLabel("Show B2");
                    getDisplayBox(sim.box).repaint();
                }
            });
            atomFilterButtonB2.setLabel("Hide B2");
        } else {
            atomFilterButtonA = null;
            atomFilterButtonB1 = null;
            atomFilterButtonB2 = null;
        }

        // *******  Edit potential parameters; sigma, epsilon, lambda
        JPanel speciesEditors = new JPanel(new GridLayout(0, 1));
        JPanel controls = new JPanel(new GridBagLayout());

        speciesEditors.add(ABox.graphic());
        speciesEditors.add(BBox.graphic());
        speciesEditors.add(nB1.graphic());
        speciesEditors.add(nB2.graphic());

//        epsilonSliders.add(ABSlider.graphic(), vertGBC);
//        epsilonSliders.add(solventThermoFrac.graphic(), vertGBC);

        JPanel AASliders = potentialEditors(sim.p2AA);
        JPanel AB1Sliders = potentialEditors(sim.p2AB1);
        JPanel AB2Sliders = potentialEditors(sim.p2AB2);
        JPanel B1B1Sliders = potentialEditors(sim.p2B1B1);
        JPanel B1B2Sliders = potentialEditors(sim.p2B1B2);
        JPanel B2B2Sliders = potentialEditors(sim.p2B2B2);


        // ***********  Assemble controls
        final JTabbedPane sliderPanel = new JTabbedPane();
        //panel for all the controls
        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);
        sliderPanel.add(controls, "Controls");
        sliderPanel.add(speciesEditors, "Number of Molecules");
        sliderPanel.add(AASliders,"A-A");
        sliderPanel.add(AB1Sliders,"A-B1");
        sliderPanel.add(AB2Sliders,"A-B2");
        sliderPanel.add(B1B1Sliders,"B1-B1");
        sliderPanel.add(B1B2Sliders,"B1-B2");
        sliderPanel.add(B2B2Sliders,"B2-B2");
        controls.add(delaySlider.graphic(), vertGBC);
        if (space.D() == 3) {
            controls.add(atomFilterButtonA.graphic(), vertGBC);
            controls.add(atomFilterButtonB1.graphic(), vertGBC);
            controls.add(atomFilterButtonB2.graphic(), vertGBC);
        }

    }

    public static void main(String[] args) {
        int D = 2;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-dim") && i + 1 < args.length) {
                i++;
                D = Integer.parseInt(args[i]);
            }
        }
        SelfAssemblySim sim = new SelfAssemblySim(Space.getInstance(D));
        SelfAssemblyGraphic graphic = new SelfAssemblyGraphic(sim, sim.getSpace());
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }

    public JPanel potentialEditors(P2SquareWell p) {
        JPanel panel = new JPanel(new GridBagLayout());

        DeviceSlider epsSlider = slider(0.0,500.0,"epsilon", 0, p);
        epsSlider.setUnit(Kelvin.UNIT);

        DeviceSlider sigSlider = slider(0.0,2.0,"coreDiameter", 2, p);

        DeviceSlider lamSlider = slider(1.0,2.0,"lambda", 2, p);
        lamSlider.setNMajor(5);

        panel.add(sigSlider.graphic(),SimulationPanel.getVertGBC());
        panel.add(epsSlider.graphic(),SimulationPanel.getVertGBC());
        panel.add(lamSlider.graphic(),SimulationPanel.getVertGBC());
        return panel;
    }

    public DeviceSlider slider(double eMin, double eMax, String s, int precision, P2SquareWell p){

        DeviceSlider AASlider = new DeviceSlider(sim.getController(), p, s, precision);
//        AASlider.setPrecision(precision);
        double value = AASlider.getValue();
        System.out.println(s+" "+value);
        AASlider.doUpdate();
        AASlider.setShowBorder(true);
        AASlider.setLabel(s);
        AASlider.setMinimum(eMin);
        AASlider.setMaximum(eMax);
        AASlider.setNMajor(4);
        AASlider.setValue(value);

        return AASlider;
    }

}
