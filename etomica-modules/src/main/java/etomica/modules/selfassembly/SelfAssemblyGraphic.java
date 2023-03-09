/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.selfassembly;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTest;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByElementType;
import etomica.atom.IAtom;
import etomica.graphics.*;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierSQWCore;
import etomica.modifier.ModifierSQWEpsilon;
import etomica.modifier.ModifierSQWLambda;
import etomica.potential.P2HardGeneric;
import etomica.space.Space;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Quantity;

import javax.swing.*;
import java.awt.*;

/**
 *
 */
public class SelfAssemblyGraphic extends SimulationGraphic {

    private static final String APP_NAME = "Self-assembly";
    private static final int REPAINT_INTERVAL = 2;
    boolean displayA = true;
    boolean displayB1 = true;
    boolean displayB2 = true;


    protected SelfAssemblySim sim;

    public SelfAssemblyGraphic(SelfAssemblySim simulation) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;

        getDisplayBox(sim.box).setPixelUnit(new Pixel(7));

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******

        // **** Stuff that Modifies the Simulation

        final IAction resetAction = getController().getSimRestart().getDataResetAction();

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController());


        DisplayTimer displayTimer = new DisplayTimer(sim.integratorHard);
        add(displayTimer);


        final IAction resetData = new IAction() {
            public void actionPerformed() {
                sim.integratorHard.resetStepCount();
            }
        };

        getController().getResetAveragesButton().setLabel("Reset");
        getController().getResetAveragesButton().setPostAction(resetData);

        ((SimulationRestart) getController().getReinitButton().getAction()).setConfiguration(sim.config);
        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                getDisplayBox(sim.box).repaint();
                resetData.actionPerformed();
            }
        });

        DeviceThermoSlider temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integratorHard);
        temperatureSelect.setUnit(Kelvin.UNIT);
        temperatureSelect.setMaximum(2000);
        temperatureSelect.setIsothermal();
        temperatureSelect.setSliderPostAction(resetAction);
        temperatureSelect.setRadioGroupPostAction(resetAction);

        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.getTypeA(), new Color(200, 0, 0));
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
        DeviceBox ABox = new DeviceBox(sim.getController());
        ABox.setInteger(true);
        ABox.setLabel("Solvent");
        ABox.setModifier(new Modifier() {
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "nA";
            }

            public double getValue() {
                return sim.getNA();
            }

            public void setValue(double newValue) {
                sim.setNA((int) newValue);
            }
        });
        ABox.setPostAction(reconfig);
        DeviceBox BBox = new DeviceBox(sim.getController());
        BBox.setInteger(true);
        BBox.setLabel("Polymer");
        BBox.setModifier(new Modifier() {
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "n";
            }

            public double getValue() {
                return sim.getNB();
            }

            public void setValue(double newValue) {
                sim.setNB((int) newValue);
            }
        });
        BBox.setPostAction(reconfig);

        DeviceBox nB1 = new DeviceBox(sim.getController());
        nB1.setInteger(true);
        nB1.setLabel("B1 segments (blue)");
        nB1.setModifier(new Modifier() {
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "n";
            }

            public double getValue() {
                return sim.getNB1();
            }

            public void setValue(double newValue) {
                sim.setNB1((int) newValue);
            }
        });
        nB1.setPostAction(reconfig);

        DeviceBox nB2 = new DeviceBox(sim.getController());
        nB2.setInteger(true);
        nB2.setLabel("B2 segments (green)");
        nB2.setModifier(new Modifier() {
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "n";
            }

            public double getValue() {
                return sim.getNB2();
            }

            public void setValue(double newValue) {
                sim.setNB2((int) newValue);
            }
        });
        nB2.setPostAction(reconfig);

//        tBox.setUnit(Kelvin.UNIT);
//        tBox.setLabel("Measured Temperature");
//        tBox.setLabelPosition(CompassDirection.NORTH);

        // ******** Turn off/on display of some atoms (3D only) **********
        AtomTest displayTest = new AtomTest() {
            public boolean test(IAtom a) {
                AtomType atomType = a.getType();
                if (atomType == sim.typeA) return displayA;
                else if (atomType == sim.typeB1) return displayB1;
                else return displayB2;
            }
        };

        final DeviceButton atomFilterButtonA, atomFilterButtonB1, atomFilterButtonB2;
        boolean showAtomFilterButtons = space.D() == 3 || space.D() == 2;
        if (showAtomFilterButtons) {
            getDisplayBox(sim.box).setAtomTestDoDisplay(displayTest);
            atomFilterButtonA = new DeviceButton(sim.getController());
            atomFilterButtonA.setAction(new IAction() {
                public void actionPerformed() {
                    displayA = !displayA;
                    if (displayA) atomFilterButtonA.setLabel("Hide A");
                    else atomFilterButtonA.setLabel("Show A");
                    getDisplayBox(sim.box).repaint();
                }
            });
            atomFilterButtonA.setLabel("Hide A");

            atomFilterButtonB1 = new DeviceButton(sim.getController());
            atomFilterButtonB1.setAction(new IAction() {
                public void actionPerformed() {
                    displayB1 = !displayB1;
                    if (displayB1) atomFilterButtonB1.setLabel("Hide B1");
                    else atomFilterButtonB1.setLabel("Show B1");
                    getDisplayBox(sim.box).repaint();
                }
            });
            atomFilterButtonB1.setLabel("Hide B1");

            atomFilterButtonB2 = new DeviceButton(sim.getController());
            atomFilterButtonB2.setAction(new IAction() {
                public void actionPerformed() {
                    displayB2 = !displayB2;
                    if (displayB2) atomFilterButtonB2.setLabel("Hide B2");
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

        PotentialEditors AASliders = new PotentialEditors(sim.p2AA);
        PotentialEditors AB1Sliders = new PotentialEditors(sim.p2AB1);
        PotentialEditors AB2Sliders = new PotentialEditors(sim.p2AB2);
        PotentialEditors B1B2Sliders = new PotentialEditors(sim.p2B1B2);
        PotentialEditors B1B1Sliders = new PotentialEditors(sim.p2B1B1);
        PotentialEditors B2B2Sliders = new PotentialEditors(sim.p2B2B2);

        AASliders.sigSlider.setPostAction((IAction) () -> {
            ((DiameterHashByElementType) getDisplayBox(sim.box).getDiameterHash()).setDiameter("A", sim.p2AA.getCollisionDiameter(0));
            getDisplayBox(sim.box).repaint();
//            double newRange = Math.min(sim.box().getBoundary().getBoxSize().getX(0)*0.5, 2*sim.p2AA.getRange());
//            ((PotentialMasterList)sim.potentialMaster).setRange(newRange);
            setNbrRange();
            sim.integratorHard.reset();
//                reconfig.actionPerformed();
            sim.sigAA = sim.p2AA.getCollisionDiameter(0);
//                System.out.println(sim.sigAA);
        });

        B1B1Sliders.sigSlider.setPostAction(() -> {
            ((DiameterHashByElementType) getDisplayBox(sim.box).getDiameterHash()).setDiameter("B1", sim.p2B1B1.getCollisionDiameter(0));
            getDisplayBox(sim.box).repaint();
            setNbrRange();
            sim.integratorHard.reset();
            sim.sigB1B1 = sim.p2B1B1.getCollisionDiameter(0);
        });

        B2B2Sliders.sigSlider.setPostAction(() -> {
            ((DiameterHashByElementType) getDisplayBox(sim.box).getDiameterHash()).setDiameter("B2", sim.p2B2B2.getCollisionDiameter(0));
            getDisplayBox(sim.box).repaint();
            setNbrRange();
            sim.integratorHard.reset();
            sim.sigB2B2 = sim.p2B2B2.getCollisionDiameter(0);
        });

        AB1Sliders.sigSlider.setPostAction(() -> {
            setNbrRange();
            sim.integratorHard.reset();
            sim.sigAB1 = sim.p2AB1.getCollisionDiameter(0);
        });
        AB2Sliders.sigSlider.setPostAction(() -> {
            setNbrRange();
            sim.integratorHard.reset();
            sim.sigAB2 = sim.p2AB2.getCollisionDiameter(0);
        });
        B1B2Sliders.sigSlider.setPostAction(() -> {
            setNbrRange();
            sim.integratorHard.reset();
            sim.sigB1B2 = sim.p2B1B2.getCollisionDiameter(0);
        });


        final DeviceButton printParamsButton = new DeviceButton(sim.getController());
        printParamsButton.setAction(new IAction() {
            public void actionPerformed() {
                System.out.println("nA: " + sim.nA);
                System.out.println("nB: " + sim.nB);
                System.out.println("nB1: " + sim.nB1);
                System.out.println("nB2: " + sim.nB2);
                printPotentialParameters("A-A", sim.p2AA);
                printPotentialParameters("A-B1", sim.p2AB1);
                printPotentialParameters("A-B2", sim.p2AB2);
                printPotentialParameters("B1-B1", sim.p2B1B1);
                printPotentialParameters("B1-B2", sim.p2B1B2);
                printPotentialParameters("B2-B2", sim.p2B2B2);
                System.out.println("Temperature: " + temperatureSelect.getTemperature());
            }
        });
        printParamsButton.setLabel("Print parameters");


        // ***********  Assemble controls
        final JTabbedPane sliderPanel = new JTabbedPane();
        //panel for all the controls
        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);
        sliderPanel.add(controls, "Controls");
        sliderPanel.add(speciesEditors, "Number of Molecules");
        sliderPanel.add(AASliders.panel, "A-A");
        sliderPanel.add(AB1Sliders.panel, "A-B1");
        sliderPanel.add(AB2Sliders.panel, "A-B2");
        sliderPanel.add(B1B1Sliders.panel, "B1-B1");
        sliderPanel.add(B1B2Sliders.panel, "B1-B2");
        sliderPanel.add(B2B2Sliders.panel, "B2-B2");
        controls.add(delaySlider.graphic(), vertGBC);
        if (showAtomFilterButtons) {
            controls.add(atomFilterButtonA.graphic(), vertGBC);
            controls.add(atomFilterButtonB1.graphic(), vertGBC);
            controls.add(atomFilterButtonB2.graphic(), vertGBC);
        }
        controls.add(printParamsButton.graphic(), vertGBC);

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
        SelfAssemblyGraphic graphic = new SelfAssemblyGraphic(sim);
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }

    private class PotentialEditors {
        final JPanel panel = new JPanel(new GridBagLayout());
        DeviceSlider epsSlider;
        DeviceSlider sigSlider;
        DeviceSlider lamSlider;

        PotentialEditors(P2HardGeneric p) {
            epsSlider = slider(0.0, 500.0, "epsilon", 0, p);
            epsSlider.setUnit(Kelvin.UNIT);

            double mySig = (p == sim.p2AA) ? 20.0 : 2.0;
            sigSlider = slider(0.0, mySig, "coreDiameter", 2, p);

            lamSlider = slider(1.0, 2.0, "lambda", 2, p);
            lamSlider.setNMajor(5);

            panel.add(sigSlider.graphic(), SimulationPanel.getVertGBC());
            panel.add(epsSlider.graphic(), SimulationPanel.getVertGBC());
            panel.add(lamSlider.graphic(), SimulationPanel.getVertGBC());
        }
    }

    public DeviceSlider slider(double eMin, double eMax, String s, int precision, P2HardGeneric p) {

        Modifier m = null;
        switch (s) {
            case "epsilon":
                m = new ModifierSQWEpsilon(p);
                break;
            case "coreDiameter":
                m = new ModifierSQWCore(p);
                break;
            case "lambda":
                m = new ModifierSQWLambda(p);
                break;
            default:
                throw new RuntimeException("unknown prop " + s);
        }
        DeviceSlider mySlider = new DeviceSlider(sim.getController(), m, precision);
        mySlider.doUpdate();
        mySlider.setShowBorder(true);
        mySlider.setLabel(s);
        mySlider.setMinimum(eMin);
        mySlider.setMaximum(eMax);
        mySlider.setNMajor(4);
        mySlider.setShowValues(true);
        mySlider.setEditValues(true);
        mySlider.setPostAction(new IAction() {
            @Override
            public void actionPerformed() {
                sim.integratorHard.reset();
            }
        });

        return mySlider;
    }

    private void printPotentialParameters(String label, P2HardGeneric p) {
        System.out.println(label + " parameters");
        System.out.println("  sigma: " + p.getCollisionDiameter(0));
        System.out.println("epsilon: " + (-p.getEnergyForState(1)));
        System.out.println(" lambda: " + (p.getCollisionDiameter(1) / p.getCollisionDiameter(0)) + "\n");
    }

    private void setNbrRange() {
        double range = Math.max(2.0 * sim.p2AA.getRange(), 2.0 * sim.p2AB1.getRange());
        range = Math.max(range, 2.0 * sim.p2AB2.getRange());
        range = Math.max(range, 2.0 * sim.p2B1B1.getRange());
        range = Math.max(range, 2.0 * sim.p2B1B2.getRange());
        range = Math.max(range, 2.0 * sim.p2B2B2.getRange());
        range = Math.min(sim.box().getBoundary().getBoxSize().getX(0) * 0.5, range);
        sim.neighborManager.setNeighborRange(range);
    }

}
