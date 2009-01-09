package etomica.modules.droplet;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.api.IAction;
import etomica.api.IBox;
import etomica.api.IData;
import etomica.api.IVectorMutable;
import etomica.config.ConfigurationLattice;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTensor;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.LatticeCubicFcc;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Pixel;

/**
 * Graphic UI for Droplet module.  Design by Ludwig Nitsche.
 *
 * @author Andrew Schultz
 */
public class DropletGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Droplet";
    private final static int REPAINT_INTERVAL = 1;
    protected Droplet sim;
    protected final DeviceNSelector nSlider;
    
    public DropletGraphic(final Droplet simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, _space.D() == 2 ? 10*REPAINT_INTERVAL : REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

        DisplayTimer displayTimer = new DisplayTimer(sim.integrator);
        add(displayTimer);

        //add meter and display for current kinetic temperature

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
                ((PotentialMasterList)sim.integrator.getPotential()).getNeighborManager(sim.box).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }

                getController().getSimRestart().actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
            
        });

        DeviceSlider cohesionEpsilon = new DeviceSlider(sim.getController());
        cohesionEpsilon.setShowBorder(true);
        cohesionEpsilon.setModifier(new ModifierGeneral(sim.p2, "epsilon"));
        cohesionEpsilon.setMaximum(5);
        cohesionEpsilon.setShowValues(true);
        cohesionEpsilon.setLabel("cohesion epsilon");
        
        JPanel systemPanel = new JPanel(new GridBagLayout());
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        
        //************* Lay out components ****************//

        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.add("System", systemPanel);
        getPanel().controlPanel.add(tabbedPane, vertGBC);
        JPanel numMoleculesPanel = new JPanel(new GridBagLayout());
        numMoleculesPanel.add(nSlider.graphic(), vertGBC);
        tabbedPane.add("# of molecules", numMoleculesPanel);
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(cohesionEpsilon.graphic(), vertGBC);
        tabbedPane.add("Surfactant potential", potentialPanel);

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

        Droplet sim = new Droplet(sp);
        DropletGraphic swGraphic = new DropletGraphic(sim, sp);
        swGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(9));
		SimulationGraphic.makeAndDisplayFrame
		        (swGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space3D.getInstance();
	        Droplet sim = new Droplet(sp);
            DropletGraphic ljmdGraphic = new DropletGraphic(sim, sp);
            ljmdGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));

		    getContentPane().add(ljmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }
    
    /**
     * Inner class to find the total pressure of the system from the pressure
     * tensor.
     */
    public static class DataProcessorTensorSplitter extends DataProcessor {

        public DataProcessorTensorSplitter() {
            data = new DataDoubleArray(3);
        }
        
        protected IData processData(IData inputData) {
            double[] x = data.getData();
            for (int i=0; i<x.length; i++) {
                x[i] = ((DataTensor)inputData).x.component(i,i);
            }
            return data;
        }

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataDoubleArray.DataInfoDoubleArray(inputDataInfo.getLabel(), inputDataInfo.getDimension(), new int[]{inputDataInfo.getLength()});
            return dataInfo;
        }

        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            return null;
        }

        private static final long serialVersionUID = 1L;
        protected final DataDoubleArray data;
    }
    
    public static class ModifierBoxSize implements Modifier {
        public ModifierBoxSize(ISpace space, IBox box, int dim, IAction reconfig) {
            this.box = box;
            this.dim = dim;
            this.reconfig = reconfig;
            size = space.makeVector();
        }
        
        public Dimension getDimension() {
            return Length.DIMENSION;
        }

        public String getLabel() {
            return "Box size";
        }

        public double getValue() {
            return box.getBoundary().getDimensions().x(dim);
        }

        public void setValue(double newValue) {
            if (newValue <= 0) {
                throw new IllegalArgumentException("Gotta be positive");
            }
            //newValue+=0.01;
            size.E(box.getBoundary().getDimensions());
            double oldValue = size.x(dim);
            size.setX(dim, newValue);
            if (dim == 1 && size.getD() == 3) {
                size.setX(2, newValue);
            }
            box.getBoundary().setDimensions(size);
            try {
                reconfig.actionPerformed();
            }
            catch (RuntimeException e) {
                // box is too small.  restore to original size
                size.setX(dim, oldValue);
                box.getBoundary().setDimensions(size);
                // and reconfig.  this shouldn't throw.
                reconfig.actionPerformed();
            }
        }
        
        protected final IBox box;
        protected final int dim;
        protected final IAction reconfig;
        protected final IVectorMutable size;
    }
}


