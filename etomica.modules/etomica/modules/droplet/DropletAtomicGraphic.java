package etomica.modules.droplet;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.api.IAction;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.LatticeCubicFcc;
import etomica.modifier.Modifier;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Null;
import etomica.units.Pixel;

/**
 * Graphic UI for Droplet module.  Design by Ludwig Nitsche.
 *
 * @author Andrew Schultz
 */
public class DropletAtomicGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Droplet";
    private final static int REPAINT_INTERVAL = 1;
    protected DropletAtomic sim;
    protected final DeviceNSelector nSlider;
    
    public DropletAtomicGraphic(final DropletAtomic simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, _space.D() == 2 ? 10*REPAINT_INTERVAL : REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

        DisplayTimer displayTimer = new DisplayTimer(sim.integrator);
        add(displayTimer);

        //add meter and display for current kinetic temperature

        DeviceSlider radiusSlider = new DeviceSlider(sim.getController(), new Modifier() {

            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            public String getLabel() {
                return "Droplet Radius";
            }

            public double getValue() {
                return sim.dropRadius;
            }

            public void setValue(double newValue) {
                sim.dropRadius = newValue;
                sim.makeDropShape();

                if (sim.integrator.isInitialized()) {
                    sim.potentialMaster.getNeighborManager(sim.box).reset();
                    sim.integrator.reset();
                }

                getDisplayBox(sim.box).repaint();
            }
        });
        radiusSlider.setPrecision(1);
        radiusSlider.setMinimum(0.2);
        radiusSlider.setMaximum(0.8);
        radiusSlider.setNMajor(4);
        
        nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(2000);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        final ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                config.initializeCoordinates(sim.box);
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).getNeighborManager(sim.box).reset();
                sim.integrator.reset();

                getController().getSimRestart().actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
            
        });

        JPanel systemPanel = new JPanel(new GridBagLayout());
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        
        //************* Lay out components ****************//

        JTabbedPane tabbedPane = new JTabbedPane();
//        tabbedPane.add("System", systemPanel);
//        getPanel().controlPanel.add(tabbedPane, vertGBC);
        JPanel dropletPanel = new JPanel(new GridBagLayout());
//        numMoleculesPanel.add(nSlider.graphic(), vertGBC);
        dropletPanel.add(radiusSlider.graphic());
        tabbedPane.add("Droplet", dropletPanel);
//        JPanel potentialPanel = new JPanel(new GridBagLayout());
//        tabbedPane.add("Surfactant potential", potentialPanel);
        getPanel().controlPanel.add(tabbedPane);

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        DataFork peFork = new DataFork();
        DataPump pePump = new DataPump(meterPE, peFork);
        dataStreamPumps.add(pePump);
        sim.integrator.addIntervalAction(pePump);
        sim.integrator.setActionInterval(pePump, 10);
        AccumulatorHistory peHistory = new AccumulatorHistory();
        peFork.addDataSink(peHistory);
        peHistory.setTimeDataSource(timeCounter);
        AccumulatorAverageCollapsing peAvg = new AccumulatorAverageCollapsing();
        peFork.addDataSink(peAvg);
        DisplayTextBoxesCAE peBox = new DisplayTextBoxesCAE();
        peBox.setAccumulator(peAvg);
        peAvg.setPushInterval(1);
        add(peBox);

        MeterKineticEnergy meterKE = new MeterKineticEnergy();
        meterKE.setBox(sim.box);
        DataFork keFork = new DataFork();
        DataPump kePump = new DataPump(meterKE, keFork);
        dataStreamPumps.add(kePump);
        sim.integrator.addIntervalAction(kePump);
        sim.integrator.setActionInterval(kePump, 10);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keFork.addDataSink(keHistory);
        keHistory.setTimeDataSource(timeCounter);
        
        MeterEnergy meterE = new MeterEnergy(sim.potentialMaster, sim.box);
        DataFork eFork = new DataFork();
        DataPump ePump = new DataPump(meterE, eFork);
        dataStreamPumps.add(ePump);
        sim.integrator.addIntervalAction(ePump);
        sim.integrator.setActionInterval(ePump, 10);
        AccumulatorHistory eHistory = new AccumulatorHistory();
        eFork.addDataSink(eHistory);
        eHistory.setTimeDataSource(timeCounter);
        
        DisplayPlot ePlot = new DisplayPlot();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        eHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLabel("Energy");
        add(ePlot);
    }

    public static void main(String[] args) {
        Space sp = null;
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    sp = Space3D.getInstance();
                }
                else {
                	sp = Space2D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
        else {
        	sp = Space3D.getInstance();
        }

        DropletAtomic sim = new DropletAtomic(sp);
        DropletAtomicGraphic swGraphic = new DropletAtomicGraphic(sim, sp);
        swGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(2));
		SimulationGraphic.makeAndDisplayFrame
		        (swGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space3D.getInstance();
	        DropletAtomic sim = new DropletAtomic(sp);
            DropletAtomicGraphic ljmdGraphic = new DropletAtomicGraphic(sim, sp);
            ljmdGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(2));

		    getContentPane().add(ljmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }
}


