package etomica.modules.catalysis;

 import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.IDataSink;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergyFromIntegrator;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Pixel;
import etomica.util.Constants.CompassDirection;

/**
 * Catalysis graphical app.
 * Design by Ken Benjamin
 * 
 * @author Andrew Schultz
 */
public class CatalysisGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Square-Well Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    public DeviceBox sigBox, epsBox, lamBox, massBox;
//    public double lambda, epsilon, mass, sigma;
//    protected Unit eUnit, dUnit, pUnit, mUnit;
    protected Catalysis sim;
    
    public CatalysisGraphic(final Catalysis simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

    	this.sim = simulation;

//        Unit tUnit = Kelvin.UNIT;

//        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
//        mUnit = new UnitRatio(Gram.UNIT, Mole.UNIT);
//        lambda = 2.0;
//        epsilon = eUnit.toSim(space.D() == 3 ? 1000 : 1500);
//        mass = space.D() == 3 ? 131 : 40;
//        sigma = 4.0;
//        if (sim.getSpace().D() == 2) {
//            dUnit = new UnitRatio(Mole.UNIT, 
//                                    new MKS().area());
//            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
//            double[] exponents = new double[] {1.0, 1.0};
//            pUnit = new CompoundUnit(units, exponents);
//        }
//        else {
//            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
//            pUnit = Bar.UNIT;
//        }

        getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));

        sim.activityIntegrate.setSleepPeriod(0);

        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        lamBox = new DeviceBox();
        massBox = new DeviceBox();

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");
        
        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController());
//        tempSlider.setUnit(Kelvin.UNIT);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(5.0);
        tempSlider.setSliderMajorValues(3);
//        tempSlider.setUnit(tUnit);
        tempSlider.setAdiabatic();
        tempSlider.setIntegrator(sim.integrator);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(lamBox.graphic());
        parameterPanel.add(massBox.graphic());
        potentialPanel.add(parameterPanel,vertGBC);
        
        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");

//        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter();
//        sigModifier.setValue(sigma);
//        ModifierGeneral epsModifier = new ModifierGeneral(sim.potentialSW, "epsilon");
//        ModifierGeneral lamModifier = new ModifierGeneral(sim.potentialSW, "lambda");
//        ModifierGeneral massModifier = new ModifierGeneral(sim.species.getLeafType().getElement(),"mass");
//        sigBox.setModifier(sigModifier);
//        sigBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
//        epsBox.setUnit(eUnit);
//        epsBox.setModifier(epsModifier);
//        lamBox.setModifier(lamModifier);
//        massBox.setModifier(massModifier);
//        massBox.setUnit(mUnit);
//        sigBox.setController(sim.getController());
//        epsBox.setController(sim.getController());
//        lamBox.setController(sim.getController());
//        massBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeRadical colorScheme = new ColorSchemeRadical(sim, sim.interactionTracker.getAgentManager());
        colorScheme.setColor(sim.speciesO.getLeafType(),java.awt.Color.BLUE);
        colorScheme.setColor(sim.speciesC.getLeafType(),java.awt.Color.RED);
        colorScheme.setColor(sim.speciesSurface.getLeafType(),java.awt.Color.GRAY);
        colorScheme.setRadicalColor(sim.speciesO.getLeafType(),java.awt.Color.CYAN);
        colorScheme.setRadicalColor(sim.speciesC.getLeafType(),java.awt.Color.PINK);
        colorScheme.setFullBondColor(sim.speciesC.getLeafType(),new Color(0.5f, 0.0f, 0.5f));
        getDisplayBox(sim.box).setColorScheme(colorScheme);
		
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature(sim, sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(1);
        final AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
        temperatureAverage.setPushInterval(20);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		temperatureFork.setDataSinks(new IDataSink[]{temperatureAverage,temperatureHistory});
        final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
        tBox.setAccumulator(temperatureAverage);
		dataStreamPumps.add(temperaturePump);
//        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature (K)");
		tBox.setLabelPosition(CompassDirection.NORTH);

		// Number density box
	    MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
	    final DisplayTextBox densityBox = new DisplayTextBox();
//	    densityBox.setUnit(dUnit);
        final DataPumpListener densityPump = new DataPumpListener(densityMeter, densityBox, 10);
        sim.integrator.getEventManager().addListener(densityPump);
        dataStreamPumps.add(densityPump);
	    densityBox.setLabel("Density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        final AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        final DataPumpListener energyPump = new DataPumpListener(eMeter, energyHistory, 100);
        sim.integrator.getEventManager().addListener(energyPump);
        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        final DataPumpListener pePump = new DataPumpListener(peMeter, peFork, 100);
        sim.integrator.getEventManager().addListener(pePump);
        dataStreamPumps.add(pePump);

		MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        final DataPumpListener kePump = new DataPumpListener(keMeter, keHistory, 100);
        sim.integrator.getEventManager().addListener(kePump);
        dataStreamPumps.add(kePump);
        
        final DisplayPlot ePlot = new DisplayPlot();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History (J/mol)");
		ePlot.setDoLegend(true);
		ePlot.setLabel("Energy");
//		ePlot.setXUnit(Picosecond.UNIT);
		
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (J/mol)");

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        SimulationRestart simRestart = new SimulationRestart(sim, space, sim.getController());
        simRestart.setConfiguration(sim.config);
        nSlider.setResetAction(simRestart);
        nSlider.setSpecies(sim.speciesO);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(sim.getSpace().D() == 3 ? 500 : 168);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                final int n = (int)nSlider.getValue() > 0 ? (int)nSlider.getValue() : 1;
                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40/n);
//                eExcludeOverlap.numAtoms = n;
//                peExcludeOverlap.numAtoms = n;
//                keExcludeOverlap.numAtoms = n;

                getDisplayBox(sim.box).repaint();
                sim.config.initializeCoordinates(sim.box);
            }
        });
        JPanel nSliderPanel = new JPanel(new GridLayout(0,1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSlider.setShowBorder(false);
        nSlider.setNMajor(4);
        nSliderPanel.add(nSlider.graphic());
        gbc2.gridx = 0;  gbc2.gridy = 1;
        statePanel.add(nSliderPanel, gbc2);

        //************* Lay out components ****************//

        getDisplayBox(sim.box).setScale(0.7);


	    ActionListener isothermalListener = new ActionListener() {
	        public void actionPerformed(ActionEvent event) {
                // we can't tell if we're isothermal here...  :(
                // if we're adiabatic, we'll re-set the temperature elsewhere
                resetDataAction.actionPerformed();
            }
        };
		tempSlider.setSliderPostAction(resetDataAction);
        tempSlider.addRadioGroupActionListener(isothermalListener);

        IAction resetAction = new IAction() {
        	public void actionPerformed() {
        	    sim.integrator.reset();

        	    // Reset density (Density is set and won't change, but
        		// do this anyway)
        		densityPump.actionPerformed();
        		densityBox.repaint();

        		// Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
                tBox.putData(temperatureAverage.getData());
                tBox.repaint();

                // IS THIS WORKING?
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        		
        		displayCycles.putData(meterCycles.getData());
        		displayCycles.repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        
        getPanel().controlPanel.add(setupPanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

    	add(ePlot);
    	add(displayCycles);
    	add(densityBox);
    	add(tBox);
    	add(peDisplay);
    	
        java.awt.Dimension d = ePlot.getPlot().getPreferredSize();
        d.width -= 50;
        ePlot.getPlot().setSize(d);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }

        CatalysisGraphic swmdGraphic = new CatalysisGraphic(new Catalysis(space), space);
		SimulationGraphic.makeAndDisplayFrame
		        (swmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            String dimStr = getParameter("dim");
            int dim = 3;
            if (dimStr != null) {
                dim = Integer.valueOf(dimStr).intValue();
            }
            Space sp = Space.getInstance(dim);
            CatalysisGraphic swmdGraphic = new CatalysisGraphic(new Catalysis(sp), sp);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
