package etomica.graphics;
import java.awt.Component;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.api.IAction;
import etomica.api.IBox;
import etomica.api.IController;
import etomica.api.IIntegrator;
import etomica.api.ISimulation;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorManagerMC;
import etomica.simulation.SimulationContainer;
import etomica.space.Space;
import etomica.units.Pixel;

/**
 * General class for graphical presentation of the elements of a molecular simulation.
 *
 * @author David Kofke
 */
 
public class SimulationGraphic implements SimulationContainer {
    
    static {
        try {
            javax.swing.UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
        } catch(Exception e) {}
    }

    public static final int GRAPHIC_ONLY = 0;
    public static final int TABBED_PANE = 1;
    private static int DEFAULT_UPDATE_INTERVAL = 100;

    private SimulationPanel simulationPanel;
    private final DeviceTrioControllerButton dcb;
    protected final ISimulation simulation;
    private int updateInterval = DEFAULT_UPDATE_INTERVAL;
    private final LinkedList displayList = new LinkedList();
    private final LinkedList deviceList = new LinkedList();
    private HashMap repaintActions = new HashMap();
    private int graphicType = GRAPHIC_ONLY;
    protected final Space space;


    public SimulationGraphic(ISimulation simulation, Space space) {
    	this(simulation, GRAPHIC_ONLY, "", DEFAULT_UPDATE_INTERVAL, space);
    }

    public SimulationGraphic(ISimulation simulation, int graphicType, Space space) {
    	this(simulation, graphicType, "", DEFAULT_UPDATE_INTERVAL, space);
    }

    public SimulationGraphic(ISimulation simulation, String appName, Space space) {
    	this(simulation, GRAPHIC_ONLY, appName, DEFAULT_UPDATE_INTERVAL, space);
    }

    public SimulationGraphic(ISimulation simulation, int graphicType, String appName, Space space) {
    	this(simulation, graphicType, appName, DEFAULT_UPDATE_INTERVAL, space);
    }

    public SimulationGraphic(ISimulation simulation, String appName, int updateInterval, Space space) {
    	this(simulation, GRAPHIC_ONLY, appName, updateInterval, space);
    }

    public SimulationGraphic(ISimulation simulation, int graphicType, String appName, int updateInterval, Space space) {
        this.simulation = simulation;
        this.space = space;
        this.updateInterval = updateInterval;
        simulationPanel = new SimulationPanel(appName);
        if(graphicType == GRAPHIC_ONLY || graphicType == TABBED_PANE) {
            this.graphicType = graphicType;
        }
        switch(this.graphicType) {
        case GRAPHIC_ONLY:
        	break;
        case TABBED_PANE:
        	getPanel().graphicsPanel.add(getPanel().tabbedPane);
        	break;
        default:
            throw new IllegalArgumentException("I don't understand graphicType "+graphicType);
        }
        dcb = new DeviceTrioControllerButton(simulation, space);
        add(dcb);
        setupDisplayBox();
    }

    public ISimulation getSimulation() {return simulation;}
    
    public final LinkedList displayList() { return displayList;}
    public final LinkedList deviceList() { return deviceList; }
        
    /**
     * A visual display of the simulation via a JPanel.
     */
	public SimulationPanel getPanel() {
	   if(simulationPanel == null) simulationPanel = new SimulationPanel();
	   return simulationPanel;
	}
     
	private void setupDisplayBox() {
	    IController controller = simulation.getController();
	    IAction[] activities = controller.getPendingActions();
	    LinkedList boxList = new LinkedList();
	    for (int i=0; i<activities.length; i++) {
	        if (activities[i] instanceof ActivityIntegrate) {
	            IIntegrator integrator = ((ActivityIntegrate)activities[i]).getIntegrator();
	            setupDisplayBox(integrator, boxList);
	        }
	    }
	}

	/**
	  * Creates a DisplayBox for each Box handled by the given Integrator 
	  * and/or its sub-integrators.  Boxs found are added to boxList.  If 
	  * a box handled by an Integrator is in BoxList, a new DisplayBox is
	  * not created.
	  */
	private void setupDisplayBox(IIntegrator integrator, LinkedList boxList) {
	    if (integrator instanceof IntegratorBox) {
	        IBox box = ((IntegratorBox)integrator).getBox();
	        if (boxList.contains(box)) return;
	        boxList.add(box);
	        final DisplayBox display = new DisplayBox(box, space, new Pixel(30));
	        add(display);
	         
	        /* For G3DSys: panel is invisible until set visible here.
	         * This kind of looks like it could possibly have solved
	         * the 'random gray panel on startup' bug. Have been
	         * unable to reproduce after adding this, anyway.
	         */
	        if(display.canvas instanceof JComponent) {
	          /* setting to false and then true just in case that's enough
	           * to fix it, since switching tabs on a gray startup will
	           * always make the panel draw properly again
	           */
	          ((JComponent)display.canvas).setVisible(false);
	          ((JComponent)display.canvas).setVisible(true);
	        }
	         
            IAction repaintAction = createDisplayBoxPaintAction(box);
	        integrator.addIntervalAction(repaintAction);
	        integrator.setActionInterval(repaintAction, updateInterval);
	        repaintActions.put(box, repaintAction);
	    }
	    else if (integrator instanceof IntegratorManagerMC) {
	        IIntegrator[] subIntegrators = ((IntegratorManagerMC)integrator).getIntegrators();
	        for (int i=0; i<subIntegrators.length; i++) {
	            setupDisplayBox(subIntegrators[i], boxList);
	        }
	    }
	
	}

	/**
	  * setPaintInterval()
	  *
	  * Sets the integrator interval between repaint actions to the value
	  * specified for the given Box.
	  */
    public void setPaintInterval(IBox box, int interval) {

    	IAction repaintAction = (IAction)repaintActions.get(box);
	    IAction[] controllerActions = simulation.getController().getAllActions();

	    for (int i = 0;  i < controllerActions.length;  i++) {
	        if (controllerActions[i] instanceof ActivityIntegrate) {

	        	IIntegrator integrator = ((ActivityIntegrate)controllerActions[i]).getIntegrator();

	            if(box == ((IntegratorBox)integrator).getBox()) {

                    IAction[] boxActions = integrator.getIntervalActions();

            	    for (int j = 0;  j < boxActions.length;  j++) {
                        if(boxActions[j] == repaintAction) {
                            integrator.setActionInterval(boxActions[j], interval);
                        }
                    }
		        }
	        }
	    }
    }

	/**
	  * getPaintAction()
	  *
	  * @return Returns the paint action associated with the given Box.
	  */
    public IAction getPaintAction(IBox box) {
	    return (IAction)repaintActions.get(box);
    }
   
    public void add(Display display) {
        final Component component = display.graphic(null);
        if(component == null) return; //display is not graphic

        if(display instanceof DisplayTextBox || display instanceof DisplayTextBoxesCAE) {
            getPanel().plotPanel.add(component, SimulationPanel.getVertGBC());
        }
        else {
            if(this.graphicType == GRAPHIC_ONLY) {
        	    getPanel().graphicsPanel.add(component);
            }
            else {
            	getPanel().tabbedPane.add(display.getLabel(), component);
            }
            //add a listener to update the tab label if the name of the display changes
            display.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
                public void propertyChange(java.beans.PropertyChangeEvent evt) {
                	if(evt.getPropertyName().equals("label")) {
                        int idx = getPanel().tabbedPane.indexOfComponent(component);
                        getPanel().tabbedPane.setTitleAt(idx,evt.getNewValue().toString());
                    }
                }
            });
        }
        displayList.add(display);
    }

    public void remove(Display display) {
        final Component component = display.graphic(null);
    	if(component == null) return; //display is not graphic
    	if(display instanceof DisplayTextBox || display instanceof DisplayTextBoxesCAE) {
    	    getPanel().plotPanel.remove(component);
    	}
    	else {
    	    if(this.graphicType == GRAPHIC_ONLY) {
    	        getPanel().graphicsPanel.remove(component);
    	    }
    	    else {
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
        if(device instanceof DeviceTable) {
            if(this.graphicType == GRAPHIC_ONLY) {
        	    getPanel().graphicsPanel.add(component);
            }
            else {
            	getPanel().tabbedPane.add(component);
            }
        }
        else {
            getPanel().controlPanel.add(component,SimulationPanel.getVertGBC());
        }
        deviceList.add(device);
    }

    public void remove(Device device) {
        final Component component = device.graphic(null);
        if(component == null) return; //display is not graphic
        if(device instanceof DeviceTable) {
            if(this.graphicType == GRAPHIC_ONLY) {
                getPanel().graphicsPanel.remove(component);           	 
            }
            else {
            	getPanel().tabbedPane.remove(component);
            }

        }
        else {
        	getPanel().controlPanel.remove(component);
        }
        deviceList.remove(device);
    }

    public DeviceTrioControllerButton getController() { return dcb; }

    private IAction createDisplayBoxPaintAction(IBox box) {
    	IAction repaintAction = null;

    	final DisplayBox display = getDisplayBox(box);
    	if(display != null) {

    		repaintAction = new Action() {
    			public void actionPerformed() {
    			    display.repaint();
    			}
    		};

    	}

    	return repaintAction;
    }

    public final JFrame makeAndDisplayFrame() {
        return makeAndDisplayFrame(getPanel());
    }
    
    public final JFrame makeAndDisplayFrame(String title) {
    	return makeAndDisplayFrame(getPanel(), title);
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
        f.setSize(700,500);
        f.getContentPane().add(panel);
        f.pack();
        return f;
    }

    public static final java.awt.event.WindowAdapter WINDOW_CLOSER 
        = new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        };
        
    public DisplayBox getDisplayBox(IBox box) {
        Iterator iterator = displayList.iterator();
        while(iterator.hasNext()) {
            Object display = iterator.next();
            if(display instanceof DisplayBox) {
                if(((DisplayBox)display).getBox() == box) return (DisplayBox)display;
            }
        }
        return null;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
//        etomica.simulation.prototypes.SwMd2D sim = new etomica.simulation.prototypes.SwMd2D();
//        etomica.simulation.prototypes.LjMd2D sim = new etomica.simulation.prototypes.LjMd2D();
//        etomica.simulation.prototypes.HsMc2d sim = new etomica.simulation.prototypes.HsMc2d();
//        etomica.simulation.prototypes.SWMD3D sim = new etomica.simulation.prototypes.SWMD3D();
//        etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D();
//        final etomica.simulation.prototypes.HSMD3DNoNbr sim = new etomica.simulation.prototypes.HSMD3DNoNbr();
//        etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
        etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
//        etomica.simulation.prototypes.HSMD2D_noNbr sim = new etomica.simulation.prototypes.HSMD2D_noNbr();
//        etomica.simulation.prototypes.GEMCWithRotation sim = new etomica.simulation.prototypes.GEMCWithRotation();
        Space space = Space.getInstance(3);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, GRAPHIC_ONLY, space);
		IAction repaintAction = simGraphic.getPaintAction(sim.getBoxs()[0]);

        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim, space));
        nSelector.setSpecies(sim.getSpeciesManager().getSpecies()[0]);
        nSelector.setBox(sim.getBoxs()[0]);
        nSelector.setPostAction(repaintAction);
        simGraphic.add(nSelector);
        simGraphic.getController().getReinitButton().setPostAction(repaintAction);
        
        simGraphic.makeAndDisplayFrame();
    }
}


