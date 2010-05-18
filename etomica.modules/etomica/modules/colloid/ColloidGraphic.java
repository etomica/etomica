package etomica.modules.colloid;

import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.action.IAction;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataDump;
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
import etomica.listener.IntegratorListenerAction;
import etomica.math.geometry.Plane;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.CriterionPositionWall;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Length;
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
        
        DeviceSlider chainLengthSlider = new DeviceSlider(sim.getController(), sim, "chainLength");
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
        
        DeviceThermoSlider thermoSlider = new DeviceThermoSlider(sim.getController());
        thermoSlider.setIntegrator(sim.integrator);
        thermoSlider.setMaximum(5);
        thermoSlider.setShowValues(true);
        add(thermoSlider);
        
        DeviceSlider boxSizeSlider = new DeviceSlider(sim.getController());
        boxSizeSlider.setModifier(new Modifier() {
            
            public void setValue(double newValue) {
                if (newValue == getValue()) return;
                IVectorMutable v = space.makeVector();
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

        JTabbedPane potentialTabs = new JTabbedPane();
        JPanel monomerPanel = new JPanel();
        DeviceBox monomerRangeBox = new DeviceBox();
        monomerRangeBox.setController(sim.getController());
        monomerRangeBox.setLabel("Wall range");
        monomerRangeBox.setModifier(new WallRangeModifier(sim.p1WallMonomer, sim.criterionWallMonomer));
        monomerPanel.add(monomerRangeBox.graphic());
        
        DeviceBox monomerEpsilonBox = new DeviceBox();
        monomerEpsilonBox.setController(sim.getController());
        monomerEpsilonBox.setLabel("Wall epsilon");
        monomerEpsilonBox.setModifier(new ModifierGeneral(sim.p1WallMonomer, "epsilon"));
        monomerPanel.add(monomerEpsilonBox.graphic());
        
        potentialTabs.add(monomerPanel, "monomer");

        JPanel colloidPanel = new JPanel();
        DeviceBox colloidRangeBox = new DeviceBox();
        colloidRangeBox.setController(sim.getController());
        colloidRangeBox.setLabel("Wall range");
        colloidRangeBox.setModifier(new WallRangeModifier(sim.p1WallColloid, null));
        colloidPanel.add(colloidRangeBox.graphic());

        DeviceBox colloidEpsilonBox = new DeviceBox();
        colloidEpsilonBox.setController(sim.getController());
        colloidEpsilonBox.setLabel("Wall epsilon");
        colloidEpsilonBox.setModifier(new ModifierGeneral(sim.p1WallColloid, "epsilon"));
        colloidPanel.add(colloidEpsilonBox.graphic());
        
        potentialTabs.add(colloidPanel, "colloid");
        
        getPanel().controlPanel.add(potentialTabs,SimulationPanel.getVertGBC());


        MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(space);
        densityProfileMeter.setProfileDim(1);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(10);
        DataPumpListener profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, 10);
        sim.integrator.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        MeterProfileByVolume colloidDensityProfileMeter = new MeterProfileByVolume(space);
        colloidDensityProfileMeter.setProfileDim(1);
        colloidDensityProfileMeter.setBox(sim.box);
        meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.speciesColloid);
        colloidDensityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed colloidDensityProfileAvg = new AccumulatorAverageFixed(10);
        colloidDensityProfileAvg.setPushInterval(10);
        profilePump = new DataPumpListener(colloidDensityProfileMeter, colloidDensityProfileAvg, 10);
        sim.integrator.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        colloidDensityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag()}, "monomer");
        profilePlot.setLegend(new DataTag[]{colloidDensityProfileAvg.getTag()}, "colloid");
        profilePlot.setDoLegend(true);
        profilePlot.setLabel("Density");

        add(profilePlot);
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
