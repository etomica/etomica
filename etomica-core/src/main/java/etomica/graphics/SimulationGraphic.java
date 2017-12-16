/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.Controller;
import etomica.box.Box;
import etomica.graphics.DisplayPlot.PopupListener;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.IntegratorListenerAction;
import etomica.simulation.Simulation;
import etomica.simulation.SimulationContainer;
import etomica.simulation.prototypes.HSMD2D;
import etomica.space.Space;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * General class for graphical presentation of the elements of a molecular simulation.
 *
 * @author David Kofke
 */

public class SimulationGraphic implements SimulationContainer {

    public static final int GRAPHIC_ONLY = 0;
    public static final int TABBED_PANE = 1;

    //anonymous class to handle window closing
    public static final java.awt.event.WindowAdapter WINDOW_CLOSER = new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
            System.exit(0);
        }
    };
    private static int DEFAULT_UPDATE_INTERVAL = 100;

    static {
        try {
            javax.swing.UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
        } catch (Exception e) {
        }
    }

    protected final Simulation simulation;
    protected final Controller controller;
    protected final Space space;
    protected final String appName;
    private final DeviceTrioControllerButton dcb;
    private final LinkedList<Display> displayList = new LinkedList<Display>();
    private final LinkedList<Device> deviceList = new LinkedList<Device>();
    protected int repaintSleep = 0;
    private SimulationPanel simulationPanel;
    private int updateInterval = DEFAULT_UPDATE_INTERVAL;
    private HashMap<Box, IntegratorListenerAction> repaintActions = new HashMap<Box, IntegratorListenerAction>();
    private int graphicType = GRAPHIC_ONLY;

    public SimulationGraphic(Simulation simulation) {
        this(simulation, GRAPHIC_ONLY, "", DEFAULT_UPDATE_INTERVAL);
    }

    public SimulationGraphic(Simulation simulation,
                             int graphicType) {
        this(simulation, graphicType, "", DEFAULT_UPDATE_INTERVAL);
    }

    public SimulationGraphic(Simulation simulation,
                             String appName) {
        this(simulation, GRAPHIC_ONLY, appName, DEFAULT_UPDATE_INTERVAL);
    }

    public SimulationGraphic(Simulation simulation,
                             int graphicType,
                             String appName) {
        this(simulation, graphicType, appName, DEFAULT_UPDATE_INTERVAL);
    }

    public SimulationGraphic(Simulation simulation,
                             String appName,
                             int updateInterval) {
        this(simulation, GRAPHIC_ONLY, appName, updateInterval);
    }

    public SimulationGraphic(Simulation simulation,
                             int graphicType,
                             String appName,
                             int updateInterval) {
        this.simulation = simulation;
        this.controller = simulation.getController();
        this.space = simulation.getSpace();
        this.updateInterval = updateInterval;
        this.appName = appName;
        simulationPanel = new SimulationPanel(appName);
        if (graphicType == GRAPHIC_ONLY || graphicType == TABBED_PANE) {
            this.graphicType = graphicType;
        }
        switch (this.graphicType) {
            case GRAPHIC_ONLY:
                break;
            case TABBED_PANE:
                getPanel().graphicsPanel.add(getPanel().tabbedPane);
                break;
            default:
                throw new IllegalArgumentException("I don't understand graphicType " + graphicType);
        }
        dcb = new DeviceTrioControllerButton(simulation, this.space, this.controller);
        add(dcb);
        setupDisplayBox(simulation.getIntegrator(), new LinkedList<Box>());
    }

    public static JFrame makeAndDisplayFrame(JPanel panel, String title) {
        JFrame f = makeFrame(panel);
        f.setTitle(title);
        f.setVisible(true);
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
        return f;
    }

    public static JFrame makeAndDisplayFrame(JPanel panel) {
        JFrame f = makeFrame(panel);
        f.setVisible(true);
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
        return f;
    }

    private static JFrame makeFrame(JPanel panel) {
        JFrame f = new JFrame();
        f.setSize(700, 500);
        f.getContentPane().add(panel);
        f.pack();
        return f;
    }

    public Simulation getSimulation() {
        return simulation;
    }

    public final LinkedList<Display> displayList() {
        return displayList;
    }

    public final LinkedList<Device> deviceList() {
        return deviceList;
    }

    /**
     * A visual display of the simulation via a JPanel.
     */
    public SimulationPanel getPanel() {
        if (simulationPanel == null) simulationPanel = new SimulationPanel();
        return simulationPanel;
    }

    /**
     * Creates a DisplayBox for each Box handled by the given Integrator and/or its sub-integrators.  Boxs found are
     * added to boxList.  If a box handled by an Integrator is in BoxList, a new DisplayBox is not created.
     */
    private void setupDisplayBox(Integrator integrator, LinkedList<Box> boxList) {
        if (integrator instanceof IntegratorBox) {
            Box box = ((IntegratorBox) integrator).getBox();
            if (boxList.contains(box)) return;
            boxList.add(box);
            final DisplayBox display = new DisplayBox(simulation, box);
            add(display);

            /* For G3DSys: panel is invisible until set visible here.
             * This kind of looks like it could possibly have solved
             * the 'random gray panel on startup' bug. Have been
             * unable to reproduce after adding this, anyway.
             */
            /* setting to false and then true just in case that's enough
             * to fix it, since switching tabs on a gray startup will
             * always make the panel draw properly again
             */
            display.canvas.setVisible(false);
            display.canvas.setVisible(true);

            IntegratorListenerAction repaintAction = new IntegratorListenerAction(createDisplayBoxPaintAction(box));
            repaintAction.setInterval(updateInterval);
            integrator.getEventManager().addListener(repaintAction);
//	        integrator.addIntervalAction(repaintAction);
//	        integrator.setActionInterval(repaintAction, updateInterval);
            repaintActions.put(box, repaintAction);
        } else if (integrator instanceof IntegratorManagerMC) {
            Integrator[] subIntegrators = ((IntegratorManagerMC) integrator).getIntegrators();
            for (int i = 0; i < subIntegrators.length; i++) {
                setupDisplayBox(subIntegrators[i], boxList);
            }
        }

    }

    /**
     * setPaintInterval()
     * <p>
     * Sets the integrator interval between repaint actions to the value specified for the given Box.
     */
    public void setPaintInterval(Box box, int interval) {
        IntegratorListenerAction repaintAction = repaintActions.get(box);
        repaintAction.setInterval(interval);
    }

    public void setRepaintSleep(int newRepaintSleep) {
        repaintSleep = newRepaintSleep;
    }

    /**
     * getPaintAction()
     *
     * @return Returns the paint action associated with the given Box.
     */
    public IAction getPaintAction(Box box) {
        return repaintActions.get(box).getAction();
    }

    public void add(final Display display) {
        final Component component = display.graphic(null);
        if (component == null) return; //display is not graphic

        if (display instanceof DisplayTextBox || display instanceof DisplayTextBoxesCAE) {
            if (this.graphicType == GRAPHIC_ONLY) {
                getPanel().plotPanel.add(component, SimulationPanel.getVertGBC());
            } else {
                getPanel().metricPanel.add(component, SimulationPanel.getVertGBC());
                if (getPanel().metricPanel.getParent() == null) {
                    getPanel().tabbedPane.add(getPanel().metricPanel, "Metrics");

                    final JPopupMenu popupMenu = new JPopupMenu();
                    final JMenuItem detachMenuItem = new JMenuItem("detach");
                    final JMenuItem reattachMenuItem = new JMenuItem("reattach");
                    detachMenuItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            JFrame f = new JFrame();
                            f.setTitle("Metrics");
                            getPanel().tabbedPane.remove(getPanel().metricPanel);
                            f.add(getPanel().metricPanel);
                            f.setSize(getPanel().metricPanel.getWidth() + 5, getPanel().metricPanel.getHeight() + 40);
                            f.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
                            f.setVisible(true);
                            popupMenu.removeAll();
                            popupMenu.add(reattachMenuItem);
                        }
                    });
                    popupMenu.add(detachMenuItem);

                    reattachMenuItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            Container f = getPanel().metricPanel.getParent();
                            while (!(f instanceof JFrame)) {
                                f = f.getParent();
                                if (f == null) {
                                    return;
                                }
                            }
                            f.remove(getPanel().metricPanel);
                            f.setVisible(false);
                            getPanel().tabbedPane.add("Metrics", getPanel().metricPanel);
                            popupMenu.removeAll();
                            popupMenu.add(detachMenuItem);
                        }
                    });

                    getPanel().metricPanel.addMouseListener(new PopupListener(popupMenu));
                }
            }
        } else {
            if (this.graphicType == GRAPHIC_ONLY) {
                getPanel().graphicsPanel.add(component);
            } else {
                addAsTab(component, display.getLabel(), !(display instanceof DisplayBox));
            }
            //add a listener to update the tab label if the name of the display changes
            display.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
                public void propertyChange(java.beans.PropertyChangeEvent evt) {
                    if (evt.getPropertyName().equals("label")) {
                        int idx = getPanel().tabbedPane.indexOfComponent(component);
                        getPanel().tabbedPane.setTitleAt(idx, evt.getNewValue().toString());
                    }
                }
            });
        }
        displayList.add(display);
    }

    public void addAsTab(final Component component, final String label, boolean detachable) {
        getPanel().tabbedPane.add(label, component);

        if (detachable) {
            final JPopupMenu popupMenu = new JPopupMenu();
            final JMenuItem detachMenuItem = new JMenuItem("detach");
            final JMenuItem reattachMenuItem = new JMenuItem("reattach");
            detachMenuItem.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFrame f = new JFrame();
                    f.setTitle(label);
                    getPanel().tabbedPane.remove(component);
                    f.add(component);
                    f.setSize(component.getWidth() + 5, component.getHeight() + 40);
                    f.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
                    f.setVisible(true);
                    popupMenu.removeAll();
                    popupMenu.add(reattachMenuItem);
                }
            });
            popupMenu.add(detachMenuItem);

            reattachMenuItem.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    Container f = component.getParent();
                    while (!(f instanceof JFrame)) {
                        f = f.getParent();
                        if (f == null) {
                            return;
                        }
                    }
                    f.remove(component);
                    f.setVisible(false);
                    getPanel().tabbedPane.add(label, component);
                    popupMenu.removeAll();
                    popupMenu.add(detachMenuItem);
                }
            });

            component.addMouseListener(new PopupListener(popupMenu));
        }
    }

    public void remove(Display display) {
        final Component component = display.graphic(null);
        if (component == null) return; //display is not graphic
        if (display instanceof DisplayTextBox || display instanceof DisplayTextBoxesCAE) {
            if (this.graphicType == GRAPHIC_ONLY) {
                getPanel().plotPanel.remove(component);
            } else {
                getPanel().metricPanel.remove(component);
                // Remove metricPanel if nothing on it
                if (getPanel().metricPanel.getComponentCount() == 0) {
                    getPanel().tabbedPane.remove(getPanel().metricPanel);
                }
            }

        } else {
            if (this.graphicType == GRAPHIC_ONLY) {
                getPanel().graphicsPanel.remove(component);
            } else {
                getPanel().tabbedPane.remove(component);
            }
        }
        displayList.remove(display);
    }

    /**
     * Adds displays graphic to the simulation display pane
     */
    public void add(Device device) {
        Component component = device.graphic(null);
        if (device instanceof DeviceTable) {
            if (this.graphicType == GRAPHIC_ONLY) {
                getPanel().graphicsPanel.add(component);
            } else {
                getPanel().tabbedPane.add(component);
            }
        } else {
            if (device instanceof DeviceTrioControllerButton) {
                getPanel().graphicsPanel.add(component, BorderLayout.SOUTH);
            } else {
                getPanel().controlPanel.add(component, SimulationPanel.getVertGBC());
            }
        }
        deviceList.add(device);
    }

    public void remove(Device device) {
        final Component component = device.graphic(null);
        if (component == null) return; //display is not graphic
        if (device instanceof DeviceTable) {
            if (this.graphicType == GRAPHIC_ONLY) {
                getPanel().graphicsPanel.remove(component);
            } else {
                getPanel().tabbedPane.remove(component);
            }

        } else {
            if (device == dcb) {
                getPanel().graphicsPanel.remove(component);
            } else {
                getPanel().controlPanel.remove(component);
            }
        }
        deviceList.remove(device);
    }

    public DeviceTrioControllerButton getController() {
        return dcb;
    }

    protected IAction createDisplayBoxPaintAction(Box box) {
        IAction repaintAction = null;

        final DisplayBox display = getDisplayBox(box);
        if (display != null) {

            repaintAction = new IAction() {
                public void actionPerformed() {
                    display.repaint();
                    if (repaintSleep != 0) {
                        try {
                            Thread.sleep(repaintSleep);
                        } catch (InterruptedException ex) {
                        }
                    }
                }
            };

        }

        return repaintAction;
    }

    public final JFrame makeAndDisplayFrame() {
        return makeAndDisplayFrame(appName);
    }

    public final JFrame makeAndDisplayFrame(String title) {
        return makeAndDisplayFrame(getPanel(), title);
    }

    public DisplayBox getDisplayBox(Box box) {
        Iterator<Display> iterator = displayList.iterator();
        while (iterator.hasNext()) {
            Object display = iterator.next();
            if (display instanceof DisplayBox) {
                if (((DisplayBox) display).getBox().getIndex() == box.getIndex()) return (DisplayBox) display;
            }
        }
        return null;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
//        etomica.simulation.prototypes.SWMD2D sim = new etomica.simulation.prototypes.SWMD2D();
//        etomica.simulation.prototypes.LJMD2D sim = new etomica.simulation.prototypes.LJMD2D();
//        etomica.simulation.prototypes.HSMC2D sim = new etomica.simulation.prototypes.HSMC2D();
//        etomica.simulation.prototypes.SWMD3D sim = new etomica.simulation.prototypes.SWMD3D();
//        etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D();
//        final etomica.simulation.prototypes.HSMD3DNoNbr sim = new etomica.simulation.prototypes.HSMD3DNoNbr();
//        etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
        HSMD2D sim = new HSMD2D();
//        etomica.simulation.prototypes.HSMD2D_noNbr sim = new etomica.simulation.prototypes.HSMD2D_noNbr();
//        etomica.simulation.prototypes.GEMCWithRotation sim = new etomica.simulation.prototypes.GEMCWithRotation();
        Space space = Space.getInstance(2);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, GRAPHIC_ONLY);
        IAction repaintAction = simGraphic.getPaintAction(sim.getBox(0));

        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpecies(sim.getSpecies(0));
        nSelector.setBox(sim.getBox(0));
        nSelector.setPostAction(repaintAction);
        simGraphic.add(nSelector);
        simGraphic.getController().getReinitButton().setPostAction(repaintAction);

        simGraphic.makeAndDisplayFrame();
    }
}


