package etomica.modules.catalysis;

 import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.IDataSink;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Bar;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.util.HistoryCollapsing;
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

        Unit tUnit = Kelvin.UNIT;

        getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));

        sim.activityIntegrate.setSleepPeriod(0);

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController());
        tempSlider.setIsothermalButtonsVisibility(false);
        tempSlider.setUnit(Kelvin.UNIT);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(1000);
        tempSlider.setSliderMajorValues(4);
        tempSlider.setUnit(tUnit);
        tempSlider.setAdiabatic();
        tempSlider.setIntegrator(sim.integrator);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        
        //display of box, timer
        ColorSchemeRadical colorScheme = new ColorSchemeRadical(sim, sim.interactionTracker.getAgentManager());
        colorScheme.setColor(sim.speciesO.getLeafType(),java.awt.Color.RED);
        colorScheme.setColor(sim.speciesC.getLeafType(),java.awt.Color.BLUE);
        colorScheme.setColor(sim.speciesSurface.getLeafType(),java.awt.Color.GRAY);
        colorScheme.setRadicalColor(sim.speciesO.getLeafType(),java.awt.Color.PINK);
        colorScheme.setRadicalColor(sim.speciesC.getLeafType(),java.awt.Color.CYAN);
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
        tBox.setUnit(tUnit);
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
	    
//		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
//        final AccumulatorHistory energyHistory = new AccumulatorHistory();
//        energyHistory.setTimeDataSource(timeCounter);
//        final DataPumpListener energyPump = new DataPumpListener(eMeter, energyHistory, 100);
//        sim.integrator.getEventManager().addListener(energyPump);
//        dataStreamPumps.add(energyPump);
//		
//		MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//        final AccumulatorHistory peHistory = new AccumulatorHistory();
//        peHistory.setTimeDataSource(timeCounter);
//        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
//        peAccumulator.setPushInterval(2);
//        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
//        final DataPumpListener pePump = new DataPumpListener(peMeter, peFork, 100);
//        sim.integrator.getEventManager().addListener(pePump);
//        dataStreamPumps.add(pePump);
//
//		MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
//        final AccumulatorHistory keHistory = new AccumulatorHistory();
//        keHistory.setTimeDataSource(timeCounter);
//        final DataPumpListener kePump = new DataPumpListener(keMeter, keHistory, 100);
//        sim.integrator.getEventManager().addListener(kePump);
//        dataStreamPumps.add(kePump);
//        
//        final DisplayPlot ePlot = new DisplayPlot();
//        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
//        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
//        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

//        ePlot.getPlot().setTitle("Energy History (J/mol)");
//		ePlot.setDoLegend(true);
//		ePlot.setLabel("Energy");
//		ePlot.setXUnit(Picosecond.UNIT);
		
//        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
//        peDisplay.setAccumulator(peAccumulator);
//        peDisplay.setLabel("Potential Energy (J/mol)");
        
        MeterDensityCO meterDensityCO = new MeterDensityCO(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        DataFork densityCOFork = new DataFork();
        DataPumpListener densityCOPump = new DataPumpListener(meterDensityCO, densityCOFork, 100);
        sim.integrator.getEventManager().addListener(densityCOPump);
        dataStreamPumps.add(densityCOPump);
        AccumulatorHistory densityCOHistory = new AccumulatorHistory(new HistoryCollapsing());
        densityCOFork.addDataSink(densityCOHistory);
        
        MeterDensityO2 meterDensityO2 = new MeterDensityO2(sim.box, sim.speciesO, sim.interactionTracker.getAgentManager());
        DataFork densityO2Fork = new DataFork();
        DataPumpListener densityO2Pump = new DataPumpListener(meterDensityO2, densityO2Fork, 100);
        sim.integrator.getEventManager().addListener(densityO2Pump);
        dataStreamPumps.add(densityO2Pump);
        AccumulatorHistory densityO2History = new AccumulatorHistory(new HistoryCollapsing());
        densityO2Fork.addDataSink(densityO2History);

        MeterDensityCO2 meterDensityCO2 = new MeterDensityCO2(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        DataFork densityCO2Fork = new DataFork();
        DataPumpListener densityCO2Pump = new DataPumpListener(meterDensityCO2, densityCO2Fork, 100);
        sim.integrator.getEventManager().addListener(densityCO2Pump);
        dataStreamPumps.add(densityCO2Pump);
        AccumulatorHistory densityCO2History = new AccumulatorHistory(new HistoryCollapsing());
        densityCO2Fork.addDataSink(densityCO2History);

        DisplayPlot densityHistoryPlot = new DisplayPlot();
        densityCOHistory.setDataSink(densityHistoryPlot.getDataSet().makeDataSink());
        densityO2History.setDataSink(densityHistoryPlot.getDataSet().makeDataSink());
        densityCO2History.setDataSink(densityHistoryPlot.getDataSet().makeDataSink());
        densityHistoryPlot.setLabel("Density");
        densityHistoryPlot.setLegend(new DataTag[]{meterDensityCO.getTag()}, "CO");
        densityHistoryPlot.setLegend(new DataTag[]{meterDensityO2.getTag()}, "O2");
        densityHistoryPlot.setLegend(new DataTag[]{meterDensityCO2.getTag()}, "CO2");
        densityHistoryPlot.setUnit(new UnitRatio(Mole.UNIT, Liter.UNIT));
        densityHistoryPlot.getPlot().setYLabel("Density (mol/L)");
        
        DataSourceWallPressureCatalysis dataSourcePressure = new DataSourceWallPressureCatalysis(space, sim.speciesC, sim.speciesO, sim.interactionTracker.getAgentManager());
        dataSourcePressure.setIntegrator(sim.integrator);
        DataSplitter pressureSplitter = new DataSplitter();
        DataPumpListener pressurePump = new DataPumpListener(dataSourcePressure, pressureSplitter, 100);
        sim.integrator.getEventManager().addListener(pressurePump);
        dataStreamPumps.add(pressurePump);
        AccumulatorHistory pressureCOHistory = new AccumulatorHistory(new HistoryCollapsing());
        pressureSplitter.setDataSink(0, pressureCOHistory);
        AccumulatorHistory pressureO2History = new AccumulatorHistory(new HistoryCollapsing());
        pressureSplitter.setDataSink(1, pressureO2History);
        AccumulatorHistory pressureCO2History = new AccumulatorHistory(new HistoryCollapsing());
        pressureSplitter.setDataSink(2, pressureCO2History);
        
        DisplayPlot pressurePlot = new DisplayPlot();
        pressureCOHistory.setDataSink(pressurePlot.getDataSet().makeDataSink());
        pressureO2History.setDataSink(pressurePlot.getDataSet().makeDataSink());
        pressureCO2History.setDataSink(pressurePlot.getDataSet().makeDataSink());
        pressurePlot.setLegend(new DataTag[]{pressureCOHistory.getTag()}, "CO");
        pressurePlot.setLegend(new DataTag[]{pressureO2History.getTag()}, "O2");
        pressurePlot.setLegend(new DataTag[]{pressureCO2History.getTag()}, "CO2");
        pressurePlot.setUnit(Bar.UNIT);
        pressurePlot.setLabel("Pressure");
        pressurePlot.getPlot().setYLabel("Pressure (bar)");

        final DeviceSlider nSliderCO = new DeviceSlider(sim.getController());
        nSliderCO.setMinimum(0);
        nSliderCO.setMaximum(500);
        nSliderCO.setShowBorder(true);
        nSliderCO.setShowValues(true);
        nSliderCO.setModifier(new ModifierGeneral(sim.config, "numCO"));
        nSliderCO.setLabel("Number of CO");
        nSliderCO.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.config.initializeCoordinates(sim.box);
                getDisplayBox(sim.box).repaint();
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).reset();
                sim.integrator.reset();
            }
        });

        final DeviceSlider nSliderO2 = new DeviceSlider(sim.getController());
        nSliderO2.setMinimum(0);
        nSliderO2.setMaximum(500);
        nSliderO2.setShowBorder(true);
        nSliderO2.setShowValues(true);
        nSliderO2.setModifier(new ModifierGeneral(sim.config, "numO2"));
        nSliderO2.setLabel("Number of O2");
        nSliderO2.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.config.initializeCoordinates(sim.box);
                getDisplayBox(sim.box).repaint();
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).reset();
                sim.integrator.reset();
            }
        });

        JPanel nSliderPanel = new JPanel(new GridLayout(0,1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSliderPanel.add(nSliderCO.graphic());
        nSliderPanel.add(nSliderO2.graphic());
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
//                peDisplay.putData(peAccumulator.getData());
//                peDisplay.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        		
        		displayCycles.putData(meterCycles.getData());
        		displayCycles.repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        
        getPanel().controlPanel.add(statePanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

//    	add(ePlot);
    	add(displayCycles);
    	add(densityBox);
    	add(tBox);
//    	add(peDisplay);
    	add(densityHistoryPlot);
    	add(pressurePlot);
    	
//        java.awt.Dimension d = ePlot.getPlot().getPreferredSize();
//        d.width -= 50;
//        ePlot.getPlot().setSize(d);
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
