package etomica.modules.interfacial;

import java.awt.GridBagConstraints;
import java.util.ArrayList;

import etomica.api.IAction;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IVector;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.Data;
import etomica.data.DataFork;
import etomica.data.DataGroupSplitter;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTensor;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Energy;
import etomica.units.Pixel;
import etomica.units.Unit;
import etomica.units.systems.LJ;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.Constants.CompassDirection;

/**
 * Graphic UI for interfacial tension module.  Design by Heath Turner.
 *
 * @author Andrew Schultz
 */
public class InterfacialGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Interfacial Tension";
    private final static int REPAINT_INTERVAL = 20;
    private DeviceThermoSlider temperatureSelect;
    protected Interfacial sim;
    protected final DeviceNSelector nSlider;
    
    public InterfacialGraphic(final Interfacial simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, _space.D() == 2 ? 10*REPAINT_INTERVAL : REPAINT_INTERVAL, _space);

        ArrayList dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

        LJ unitSystem = new LJ();
        Unit tUnit = Energy.DIMENSION.getUnit(unitSystem);

        sim.activityIntegrate.setSleepPeriod(1);
        final double expansionFac = 4;
        
        final DeviceButton expandButton = new DeviceButton(sim.getController());
        IAction expandAction = new IAction() {
            public void actionPerformed() {
                IVector dim = sim.box.getBoundary().getDimensions();
                dim.setX(0, expansionFac * dim.x(0));
                sim.box.getBoundary().setDimensions(dim);
                ((PotentialMasterList)sim.integrator.getPotential()).getNeighborManager(sim.box).reset();
                nSlider.setEnabled(false);
                expandButton.getButton().setEnabled(false);
                getDisplayBox(sim.box).repaint();
            }
        };
        expandButton.setAction(expandAction);
        expandButton.setLabel("Expand");
        add(expandButton);
        
        getController().getReinitButton().setPreAction(new IAction() {
            public void actionPerformed() {
                IVector dim = sim.box.getBoundary().getDimensions();
                dim.setX(0, dim.x(0) / expansionFac);
                sim.box.getBoundary().setDimensions(dim);
                nSlider.setEnabled(true);
                expandButton.getButton().setEnabled(true);
            }
        });
        
        IAction recenterAction = new IAction() {
            public void actionPerformed() {
                double dx = 0.5;
                int leftN1 = 0, leftP1 = 0;
                IAtomSet leafAtoms = sim.box.getLeafList();
                int nTot = leafAtoms.getAtomCount();
                for (int i=0; i<nTot; i++) {
                    IVector pos = ((IAtomPositioned)leafAtoms.getAtom(i)).getPosition();
                    if (pos.x(0) < -dx) {
                        leftN1++;
                        leftP1++;
                    }
                    else if (pos.x(0) < dx) {
                        leftP1++;
                    }
                }
                double target = 0.5 * nTot; 
                // interpolate to find median
                double median = -dx + (target - leftN1) * (2*dx) / (leftP1 - leftN1);
                // make sure we didn't extrapolate
                if (median < -dx) {
                    median = -dx;
                }
                else if (median > dx) {
                    median = dx;
                }

                for (int i=0; i<nTot; i++) {
                    IVector pos = ((IAtomPositioned)leafAtoms.getAtom(i)).getPosition();
                    pos.setX(0, pos.x(0) - median);
                }
                ((PotentialMasterList)sim.integrator.getPotential()).getNeighborManager(sim.box).reset();
            }
        };
        sim.integrator.addIntervalAction(recenterAction);
        sim.integrator.setActionInterval(recenterAction, 1000);

	    //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));
//        sim.integrator.addListener(new IntervalActionAdapter(this.getDisplayBoxPaintAction(sim.box)));

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
        DisplayTimer displayTimer = new DisplayTimer(sim.integrator);
        add(displayTimer);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        sim.integrator.addIntervalAction(temperaturePump);
        sim.integrator.setActionInterval(temperaturePump, 10);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		final DisplayTextBox tBox = new DisplayTextBox();
		temperatureFork.setDataSinks(new DataSink[]{tBox,temperatureHistory});
        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);

		dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);

		// Number density box
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
	    final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        sim.integrator.addIntervalAction(densityPump);
        sim.integrator.setActionInterval(densityPump, 10);
        dataStreamPumps.add(densityPump);
	    densityBox.setLabel("Number Density");
	    
//	      MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential(), sim.box);
//        AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        energyHistory.setTimeDataSource(timeCounter);
//        DataPump energyPump = new DataPump(eMeter, energyHistory);
//        sim.integrator.addIntervalAction(energyPump);
//        sim.integrator.setActionInterval(energyPump, 60);
//        energyHistory.setPushInterval(5);
//        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotential());
        peMeter.setBox(sim.box);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        sim.integrator.addIntervalAction(pePump);
        sim.integrator.setActionInterval(pePump, 60);
        peHistory.setPushInterval(1);
        dataStreamPumps.add(pePump);
		
        DisplayPlot ePlot = new DisplayPlot();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());

        ePlot.getPlot().setTitle("Energy History");
		ePlot.setDoLegend(true);
		ePlot.setLabel("Energy");
		
        MeterVirialTensorFromForceSum pMeter = new MeterVirialTensorFromForceSum(space, sim.forceSum);
        DataProcessorTensorSplitter tensorSplitter = new DataProcessorTensorSplitter();
        final DataPump pPump = new DataPump(pMeter, tensorSplitter);
        DataFork virialFork = new DataFork();
        tensorSplitter.setDataSink(virialFork);
        final DataSplitter splitter = new DataSplitter();
        virialFork.addDataSink(splitter);
        final AccumulatorAverageCollapsing[] pAccumulator = new AccumulatorAverageCollapsing[space.D()];
        final DisplayTextBoxesCAE[] pDisplay = new DisplayTextBoxesCAE[space.D()];
        String[] comp = new String[]{"x", "y", "z"};
        for (int i=0; i<pAccumulator.length; i++) {
            pAccumulator[i] = new AccumulatorAverageCollapsing();
            splitter.setDataSink(i, pAccumulator[i]);
            pAccumulator[i].setPushInterval(10);
            pDisplay[i] = new DisplayTextBoxesCAE();
            pDisplay[i].setLabel(comp[i]+" Virial");
            pDisplay[i].setAccumulator(pAccumulator[i]);
        }
        sim.integrator.addIntervalAction(pPump);
        sim.integrator.setActionInterval(pPump, 20);
        dataStreamPumps.add(pPump);
        
        DataProcessorInterfacialTension interfacialTension = new DataProcessorInterfacialTension(space);
        interfacialTension.setBox(sim.box);
        virialFork.addDataSink(interfacialTension);
        final AccumulatorAverageCollapsing tensionAvg = new AccumulatorAverageCollapsing();
        interfacialTension.setDataSink(tensionAvg);
        tensionAvg.setPushInterval(10);
        DisplayTextBoxesCAE tensionDisplay = new DisplayTextBoxesCAE();
        tensionDisplay.setAccumulator(tensionAvg);

        MeterVirialProfileFromForceSum virialProfileMeter = new MeterVirialProfileFromForceSum(sim.forceSum);
        DataFork virialProfileFork = new DataFork();
        DataPump virialProfilePump = new DataPump(virialProfileMeter, virialProfileFork);
        DataGroupSplitter virialSplitter = new DataGroupSplitter();
        virialProfileFork.addDataSink(virialSplitter);
        sim.integrator.addIntervalAction(virialProfilePump);
        sim.integrator.setActionInterval(virialProfilePump, 20);
        AccumulatorAverageFixed[] virialProfileAvg = new AccumulatorAverageFixed[space.D()];
        DisplayPlot virialPlot = new DisplayPlot();
        for (int i=0; i<space.D(); i++) {
            virialProfileAvg[i] = new AccumulatorAverageFixed(100);
            virialProfileAvg[i].setPushInterval(10);
            virialSplitter.setDataSink(i, virialProfileAvg[i]);
            virialProfileAvg[i].addDataSink(virialPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
            virialPlot.setLegend(new DataTag[]{virialProfileAvg[i].getTag()}, comp[i]+" Virial");
        }
        virialPlot.setLabel("Virial Profile");
        add(virialPlot);
        dataStreamPumps.add(virialProfilePump);
        
        DataProcessorInterfacialTensionProfile interfacialTensionProfile = new DataProcessorInterfacialTensionProfile(sim.forceSum);
        interfacialTensionProfile.setBox(sim.box);
        virialProfileFork.addDataSink(interfacialTensionProfile);
        AccumulatorAverageFixed tensionProfileAvg = new AccumulatorAverageFixed(100);
        interfacialTensionProfile.setDataSink(tensionProfileAvg);
        tensionProfileAvg.setPushInterval(10);
        DisplayPlot tensionPlot = new DisplayPlot();
        tensionProfileAvg.addDataSink(tensionPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        tensionPlot.setLabel("Tension Profile");
        add(tensionPlot);
        
        MeterDensityProfileFromForceSum densityProfileMeter = new MeterDensityProfileFromForceSum(sim.forceSum);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed();
        densityProfileAvg.setPushInterval(10);
        DataPump profilePump = new DataPump(densityProfileMeter, densityProfileAvg);
        sim.integrator.addIntervalAction(profilePump);
        sim.integrator.setActionInterval(profilePump, 10);
        dataStreamPumps.add(profilePump);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        profilePlot.setLabel("Density");
        add(profilePlot);


        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);

        nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(space.D() == 3 ? 2048 : 500);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                int n = (int)nSlider.getValue();
                if(n == 0) {
                	sim.integrator.setThermostatInterval(400);
                }
                else {
                    sim.integrator.setThermostatInterval((400+(n-1))/n);
                }
                
                if (oldN < n) {
                    config.initializeCoordinates(sim.box);
                }
                oldN = n;
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
                getController().getSimRestart().actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
            
            ConfigurationLattice config = new ConfigurationLattice((space.D() == 2) ? new LatticeOrthorhombicHexagonal() : new LatticeCubicFcc(), space);
            int oldN = sim.box.getMoleculeList().getAtomCount();
        });

        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getDisplayBox(sim.box).setScale(0.7);

        //temperature selector
        temperatureSelect = new DeviceThermoSlider(sim.getController());
        temperatureSelect.setPrecision(2);
        temperatureSelect.setMinimum(0.0);
        if (space.D() == 3) {
            // Tc = 1.312
            temperatureSelect.setMaximum(1.5);
        }
        else {
            // Tc = 0.515
            temperatureSelect.setMaximum(0.6);
        }
        temperatureSelect.setSliderMajorValues(3);
	    temperatureSelect.setUnit(tUnit);
	    temperatureSelect.setIntegrator(sim.integrator);
	    temperatureSelect.setIsothermalButtonsVisibility(false);

        IAction resetAction = new IAction() {
        	public void actionPerformed() {

                // Reset density (Density is set and won't change, but
        		// do this anyway)
        		densityPump.actionPerformed();
        		densityBox.repaint();

        		// Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
//                tBox.putData(temperatureHistory.getData());
                tBox.repaint();

                // IS THIS WORKING?
                pPump.actionPerformed();
                for (int i=0; i<space.D(); i++) {
                    pDisplay[i].putData(pAccumulator[i].getData());
                    pDisplay[i].repaint();
                }
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        add(nSlider);

    	add(ePlot);
    	add(densityBox);
    	add(tBox);
    	for (int i=0; i<space.D(); i++) {
    	    add(pDisplay[i]);
    	}
    	add(peDisplay);
        add(tensionDisplay);

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

        Interfacial sim = new Interfacial(sp);
        InterfacialGraphic ljmdGraphic = new InterfacialGraphic(sim, sp);
        ljmdGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(9));
		SimulationGraphic.makeAndDisplayFrame
		        (ljmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space2D.getInstance();
	        Interfacial sim = new Interfacial(sp);
            InterfacialGraphic ljmdGraphic = new InterfacialGraphic(sim, sp);
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
        
        protected Data processData(Data inputData) {
            double[] x = data.getData();
            for (int i=0; i<x.length; i++) {
                x[i] = ((DataTensor)inputData).x.component(i,i);
            }
            return data;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            dataInfo = new DataDoubleArray.DataInfoDoubleArray(inputDataInfo.getLabel(), inputDataInfo.getDimension(), new int[]{inputDataInfo.getLength()});
            return dataInfo;
        }

        public DataPipe getDataCaster(IDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            return null;
        }

        private static final long serialVersionUID = 1L;
        protected final DataDoubleArray data;
    }
}


