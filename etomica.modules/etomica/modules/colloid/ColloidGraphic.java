package etomica.modules.colloid;

import java.util.ArrayList;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.listener.IntegratorListenerAction;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Pixel;

/**
 * Colloid module app.  Design by Alberto Striolo.
 * 
 * @author Andrew Schultz
 */
public class ColloidGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Colloid";
    private final static int REPAINT_INTERVAL = 10;
    protected ColloidSim sim;
    
    public ColloidGraphic(final ColloidSim simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

    	this.sim = simulation;

    	((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), sim.p2mm.getCoreDiameter());
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesColloid.getLeafType(), 7.5);
    	
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
        graftSelector.addOption("20", new GraftAction(20));
        graftSelector.setSelected(6);
        graftSelector.setLabel("# of grafted chains");
        add(graftSelector);
        
        DeviceSlider chainLengthSlider = new DeviceSlider(sim.getController(), sim, "chainLength");
        chainLengthSlider.setLabel("Chain length");
        chainLengthSlider.setShowBorder(true);
        chainLengthSlider.setMaximum(200);
        chainLengthSlider.setNMajor(4);
        chainLengthSlider.setPostAction(getPaintAction(sim.box));
        chainLengthSlider.setShowValues(true);
        add(chainLengthSlider);
        
        MeterTemperature meterTemperature = new MeterTemperature(sim.box, 3);
        DisplayTextBox displayTemperature = new DisplayTextBox();
        DataPumpListener tempPump = new DataPumpListener(meterTemperature, displayTemperature);
        sim.getIntegrator().getEventManager().addListener(tempPump);
        add(displayTemperature);
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
