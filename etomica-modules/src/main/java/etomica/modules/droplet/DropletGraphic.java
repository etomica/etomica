/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTensor;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Pixel;
import etomica.units.SimpleUnit;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Time;
import g3dsys.images.Ellipse;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

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
    protected Ellipse ellipse = null;

    public DropletGraphic(final Droplet simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, _space.D() == 2 ? 10*REPAINT_INTERVAL : REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

    	getController().getSimRestart().setConfiguration(sim.config);
    	getController().getReinitButton().setPostAction(new IAction() {
    	    public void actionPerformed() {
    	        getDisplayBox(sim.box).repaint();
    	    }
    	});
    	
    	
        IAction recenterAction = new IAction() {
            public void actionPerformed() {
                IAtomList leafList = sim.box.getLeafList();
                center.E(0);
                for (int i=0; i<leafList.getAtomCount(); i++) {
                    center.PE(leafList.getAtom(i).getPosition());
                }
                center.TE(1.0/leafList.getAtomCount());

                for (int i=0; i<leafList.getAtomCount(); i++) {
                    leafList.getAtom(i).getPosition().ME(center);
                }
            }
            final Vector center = sim.getSpace().makeVector();
        };
        IntegratorListenerAction recenterActionListener = new IntegratorListenerAction(recenterAction);
        sim.integrator.getEventManager().addListener(recenterActionListener);
        recenterActionListener.setInterval(10);


        DisplayTimer displayTimer = new DisplayTimer(sim.integrator);
        // die unit symbol, die!
        displayTimer.setUnit(new SimpleUnit(Time.DIMENSION, 1, "time", "", true));
        displayTimer.setLabel("Simulation time");
        getPanel().controlPanel.add(displayTimer.graphic(), SimulationPanel.getVertGBC());
        displayTimer.setUpdateInterval(1);
        
        double timeStep = sim.integrator.getTimeStep();
        DeviceSlider timeStepSlider = new DeviceSlider(sim.getController(), new ModifierGeneral(sim.integrator, "timeStep"));
        timeStepSlider.setPrecision(1);
        timeStepSlider.setMinimum(0);
        timeStepSlider.setMaximum(1);
        timeStepSlider.setNMajor(4);
        timeStepSlider.setValue(timeStep);
        timeStepSlider.setShowValues(true);
        timeStepSlider.setLabel("time step");
        timeStepSlider.setShowBorder(true);
        

        //add meter and display for current kinetic temperature

        nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(2000);
        nSlider.setNMajor(4);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);

        final DiameterHashByType diameterHash = (DiameterHashByType)getDisplayBox(sim.box).getDiameterHash();
        diameterHash.setDiameter(sim.species.getLeafType(), Math.pow(2.0/nSlider.getValue(),1.0/3.0));
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.config.initializeCoordinates(sim.box);
                diameterHash.setDiameter(sim.species.getLeafType(), Math.pow(2.0/nSlider.getValue(),1.0/3.0));
                getDisplayBox(sim.box).repaint();
            }
        });
        
        DeviceSlider defSlider = new DeviceSlider(sim.getController(), new ModifierGeneral(sim.config, "deformation"));
        defSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.config.initializeCoordinates(sim.box);
                getDisplayBox(sim.box).repaint();
            }
        });
        defSlider.setPrecision(1);
        defSlider.setMinimum(-0.5);
        defSlider.setMaximum(0.5);
        defSlider.setNMajor(5);
        defSlider.setShowValues(true);
        defSlider.setLabel("initial deformation");
        defSlider.setShowBorder(true);

        final DeviceSlider cohesionEpsilon = new DeviceSlider(sim.getController());
        cohesionEpsilon.setShowBorder(true);
        cohesionEpsilon.setModifier(new ModifierGeneral(sim.p2, "epsilon"));
        cohesionEpsilon.setPrecision(2);
        cohesionEpsilon.setMinimum(0.4);
        cohesionEpsilon.setMaximum(1);
        cohesionEpsilon.setNMajor(3);
        cohesionEpsilon.setShowValues(true);
        cohesionEpsilon.setLabel("Cohesion Strength");

        DeviceSlider squeeze = new DeviceSlider(sim.getController());
        squeeze.setShowBorder(true);
        squeeze.setModifier(new ModifierGeneral(sim.p1Smash, "g"));
        squeeze.setPrecision(1);
        squeeze.setMaximum(5);
        squeeze.setNMajor(5);
        squeeze.setShowValues(true);
        squeeze.setLabel("Squeezing Force");
        final EllipseDisplayAction ellipseDisplayAction = new EllipseDisplayAction(this, 1);
        squeeze.setPostAction(new IAction() {
            public void actionPerformed() {
                ellipseDisplayAction.displayEllipse(sim.p1Smash.g);
            }
        });


        ModifierBoolean surfaceCohesionModifier = new ModifierBoolean(){
            public void setBoolean(boolean b) {
                sim.p2.setUseSurfaceOnly(b);
                if (!b) {
                    if (sim.p2.getEpsilon() > 1) {
                        sim.p2.setEpsilon(1);
                        cohesionEpsilon.doUpdate();
                    }
                    cohesionEpsilon.setPrecision(2);
                    cohesionEpsilon.setMinimum(0.4);
                    cohesionEpsilon.setMaximum(1);
                    cohesionEpsilon.setNMajor(3);
                }
                else {
                    cohesionEpsilon.setMaximum(20);
                    cohesionEpsilon.setMinimum(0);
                    cohesionEpsilon.setNMajor(4);
                    cohesionEpsilon.setPrecision(1);
                }
                cohesionEpsilon.doUpdate();
            }
            public boolean getBoolean() {return sim.p2.getUseSurfaceOnly();}
        };
        DeviceToggleButton surfaceCohesionButton = new DeviceToggleButton(sim.getController(), surfaceCohesionModifier, "Use Bulk Cohesion", "Use Surface Cohesion");
        
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        JPanel controlsPanel = new JPanel(new GridBagLayout());
        controlsPanel.add(timeStepSlider.graphic(), vertGBC);
        
        //************* Lay out components ****************//

        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.add("Controls", controlsPanel);
        getPanel().controlPanel.add(tabbedPane, vertGBC);
        JPanel numMoleculesPanel = new JPanel(new GridBagLayout());
        numMoleculesPanel.add(nSlider.graphic(), vertGBC);
        numMoleculesPanel.add(defSlider.graphic(), vertGBC);
        tabbedPane.add("Configuration", numMoleculesPanel);
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(surfaceCohesionButton.graphic(), vertGBC);
        potentialPanel.add(cohesionEpsilon.graphic(), vertGBC);
        potentialPanel.add(squeeze.graphic(), vertGBC);
        tabbedPane.add("Potential", potentialPanel);
        
        DataSplitter splitter = new DataSplitter();
        DataPump deformationPump = new DataPump(sim.meterDeformation, splitter);
        dataStreamPumps.add(deformationPump);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(deformationPump));

        DataSourceCountTime timer = new DataSourceCountTime(sim.integrator);
        
        DataFork radiusFork = new DataFork();
        splitter.setDataSink(0, radiusFork);
        AccumulatorHistory radiusHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
        radiusHistory.setTimeDataSource(timer);
        radiusFork.addDataSink(radiusHistory);
        DisplayPlot radiusPlot = new DisplayPlot();
        radiusPlot.setLabel("Radius");
        radiusHistory.setDataSink(radiusPlot.getDataSet().makeDataSink());
        radiusPlot.setDoLegend(false);
        add(radiusPlot);

        DataFork deformationFork = new DataFork();
        splitter.setDataSink(1, deformationFork);
        AccumulatorHistory deformationHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
        deformationHistory.setTimeDataSource(timer);
        deformationFork.addDataSink(deformationHistory);
        DisplayPlot deformPlot = new DisplayPlot();
        deformPlot.setLabel("Deformation");
        deformationHistory.setDataSink(deformPlot.getDataSet().makeDataSink());
        deformPlot.setDoLegend(false);
        add(deformPlot);
        
        ColorSchemeDropletSurface colorScheme2 = new ColorSchemeDropletSurface(sim.liquidFilter);
        getDisplayBox(sim.box).setColorScheme(colorScheme2);
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
        DropletGraphic simGraphic = new DropletGraphic(sim, sp);
		SimulationGraphic.makeAndDisplayFrame
		        (simGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space3D.getInstance();
	        Droplet sim = new Droplet(sp);
            DropletGraphic simGraphic = new DropletGraphic(sim, sp);
            simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));

		    getContentPane().add(simGraphic.getPanel());
	    }
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

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            dataInfo = new DataDoubleArray.DataInfoDoubleArray(inputDataInfo.getLabel(), inputDataInfo.getDimension(), new int[]{inputDataInfo.getLength()});
            return dataInfo;
        }

        protected final DataDoubleArray data;
    }
    
    public static class ModifierBoxSize implements Modifier {
        public ModifierBoxSize(Space space, Box box, int dim, IAction reconfig) {
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
            return box.getBoundary().getBoxSize().getX(dim);
        }

        public void setValue(double newValue) {
            if (newValue <= 0) {
                throw new IllegalArgumentException("Gotta be positive");
            }
            //newValue+=0.01;
            size.E(box.getBoundary().getBoxSize());
            double oldValue = size.getX(dim);
            size.setX(dim, newValue);
            if (dim == 1 && size.getD() == 3) {
                size.setX(2, newValue);
            }
            box.getBoundary().setBoxSize(size);
            try {
                reconfig.actionPerformed();
            }
            catch (RuntimeException e) {
                // box is too small.  restore to original size
                size.setX(dim, oldValue);
                box.getBoundary().setBoxSize(size);
                // and reconfig.  this shouldn't throw.
                reconfig.actionPerformed();
            }
        }
        
        protected final Box box;
        protected final int dim;
        protected final IAction reconfig;
        protected final Vector size;
    }
}


