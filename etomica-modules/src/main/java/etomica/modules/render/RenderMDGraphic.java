/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.render;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomPair;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.IDataSink;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ActionConfigWindow;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modules.render.RenderMD.RenderMDParam;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Null;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Unit;
import etomica.units.systems.LJ;
import etomica.util.Constants.CompassDirection;

public class RenderMDGraphic extends SimulationGraphic {

    private final static String APP_NAME = "MD of Self-Assembly";
    private final static int REPAINT_INTERVAL = 20;
    private DeviceThermoSlider temperatureSelect;
    protected RenderMD sim;
    
    private boolean showConfig = false;

    public RenderMDGraphic(final RenderMD simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
    	
    	RenderMDParam params = new RenderMDParam();

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        
    	this.sim = simulation;

        LJ unitSystem = new LJ();
        Unit tUnit = new PrefixedUnit(Prefix.CENTI,Energy.DIMENSION.getUnit(unitSystem));

        sim.activityIntegrate.setSleepPeriod(0);
        
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getAtomType(0),0.1);

       
	    //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType(sim));
//        sim.integrator.addListener(new IntervalActionAdapter(this.getDisplayBoxPaintAction(sim.box)));

	    //meters and displays
        
        final IAction resetDataAction = new IAction(){
            public void actionPerformed() {
                getController().getSimRestart().getDataResetAction();
            }
        };

				
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(10);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
		final DisplayTextBox tBox = new DisplayTextBox();
		temperatureFork.setDataSinks(new IDataSink[]{tBox,temperatureHistory});
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
        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
        sim.integrator.getEventManager().addListener(densityPumpListener);
        densityPumpListener.setInterval(10);
        dataStreamPumps.add(densityPump);
	    densityBox.setLabel("Number Density");
	    
		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataPump energyPump = new DataPump(eMeter, energyHistory);
        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        sim.integrator.getEventManager().addListener(energyPumpListener);
        energyPumpListener.setInterval(60);
        energyHistory.setPushInterval(5);
        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());
        peMeter.setBox(sim.box);
        AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        sim.integrator.getEventManager().addListener(pePumpListener);
        pePumpListener.setInterval(60);
        peHistory.setPushInterval(5);
        dataStreamPumps.add(pePump);
		
		MeterKineticEnergy keMeter = new MeterKineticEnergy();
        keMeter.setBox(sim.box);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        DataFork keFork = new DataFork();
        DataPump kePump = new DataPump(keMeter, keFork);
        keFork.addDataSink(keHistory);
        final AccumulatorAverage keAvg = new AccumulatorAverageCollapsing();
        keFork.addDataSink(keAvg);
        IntegratorListenerAction kePumpListener = new IntegratorListenerAction(kePump);
        sim.integrator.getEventManager().addListener(kePumpListener);
        kePumpListener.setInterval(60);
        keHistory.setPushInterval(5);
        dataStreamPumps.add(kePump);
        
        DisplayPlot ePlot = new DisplayPlot();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History");
		ePlot.setDoLegend(true);
		ePlot.setLabel("Energy");
		
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        
        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getDisplayBox(sim.box).setScale(0.7);

        //temperature selector
        temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
        temperatureSelect.setPrecision(3);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(1.0);
        temperatureSelect.setSliderMajorValues(3);
//	    temperatureSelect.setUnit(tUnit);
//	    temperatureSelect.setIsothermal();

        final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
                resetDataAction.actionPerformed();
		    }
		};

		temperatureSelect.setSliderPostAction(temperatureAction);
        temperatureSelect.setRadioGroupPostAction(temperatureAction);

		DeviceSlider lambdaSlider = new DeviceSlider(null);
        lambdaSlider = new DeviceSlider(null);
        lambdaSlider.setShowValues(true);
        lambdaSlider.setEditValues(true);
        lambdaSlider.setPrecision(2);
        lambdaSlider.setMinimum(1);
        lambdaSlider.setMaximum(6);
        lambdaSlider.setNMajor(6);
        lambdaSlider.setValue(params.lambda);
        
        lambdaSlider.setModifier(new Modifier(){

            public Dimension getDimension() {
                return Null.DIMENSION;
            }

            public String getLabel() {
                return "Tightness";
            }

            public double getValue() {
                return sim.potentialBonded.getLambda();
            }

            public void setValue(double value) {
                sim.potentialBonded.setLambda(value);
                
            }
            
        });

        // panel for lambda control / display
        JPanel lambdaSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        lambdaSlider.setShowBorder(false);
        lambdaSliderPanel.add(lambdaSlider.graphic());


        DeviceSlider epsilonSlider = new DeviceSlider(null);
        epsilonSlider.setShowValues(true);
        epsilonSlider.setEditValues(true);
        epsilonSlider.setMinimum(0);
        epsilonSlider.setMaximum(10000);
        epsilonSlider.setNMajor(3);
        epsilonSlider.setNMinor(1);
        epsilonSlider.setValue(new RenderMDParam().lambda);
        
        epsilonSlider.setModifier(new Modifier(){

            public Dimension getDimension() {
                return Energy.DIMENSION;
            }

            public String getLabel() {
                return "Well depth";
            }

            public double getValue() {
                return sim.potentialBonded.getEpsilon();
            }

            public void setValue(double value) {
                sim.potentialBonded.setEpsilon(value);
                
            }
            
        });

        // panel for epsilon control / display
        JPanel epsilonSliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        epsilonSlider.setShowBorder(false);
        epsilonSliderPanel.add(epsilonSlider.graphic());
 
        // Add all state page sub panels onto a single panel
        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(temperatureSelect.graphic(), gbc2);
        gbc2.gridx = 0;  gbc2.gridy = 1;
        statePanel.add(lambdaSliderPanel, gbc2);
        gbc2.gridx = 0;  gbc2.gridy = 2;
        statePanel.add(epsilonSliderPanel, gbc2);

        getPanel().controlPanel.add(statePanel, vertGBC);

        // show config button
        DeviceButton configButton = new DeviceButton(sim.getController());
        configButton.setLabel("Show Config");
        configButton.setAction(new ActionConfigWindow(sim.box));

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

                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

//        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);

        if(showConfig == true) {
            add(configButton);
        }

    	add(densityBox);
    	add(tBox);
    	add(peDisplay);
        //sim.integrator.setTemperature(new RenderMDParam().temperature);
    	temperatureSelect.setTemperature(new RenderMDParam().temperature);

//    	Integer paintInterval= ((IntegratorListenerAction)getPaintAction(sim.box)).getInterval();
    	if(params.drawBonds) {
        	IntegratorListenerAction bondAction = new IntegratorListenerAction(new IAction() {
                Map<IAtom,Set<IAtom>> bondedSet = sim.criterion.bondedSet;
        	    DisplayBoxCanvasG3DSys dis = ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas);
                Map<IAtomList,Double> bondedMap = sim.potentialBonded.bondMap;
                AtomPair pair = new AtomPair();
                ArrayList allBonds = new ArrayList();
    
        	    public void actionPerformed() {
        	        IAtomList atoms = simulation.box.getLeafList();
        	        double lambda2 = sim.potentialBonded.getLambda()*sim.potentialBonded.getLambda();
        	        int nAtoms = atoms.getAtomCount();
        	        int nBonds = allBonds.size();
        	        for(int i=0; i<nBonds; i++) {
        	            dis.releaseBond(allBonds.get(i));
        	        }
        	        allBonds.clear();
        	        for(int i=0; i<nAtoms; i++) {
        	            pair.atom0 = atoms.getAtom(i);
        	            Set<IAtom> aSet = bondedSet.get(pair.atom0);
        	            Iterator<IAtom> iterator = aSet.iterator();
        	            while(iterator.hasNext()) {
        	                pair.atom1 = iterator.next();
        	                double bondLength2 = bondedMap.get(pair).doubleValue();
        	                double separation2 = pair.atom1.getPosition().Mv1Squared(pair.atom0.getPosition());
        	                if(separation2 > 0.999*bondLength2 && separation2 < 1.001*bondLength2*lambda2) {
        	                    allBonds.add(dis.makeBond(pair, null));
        	                }
        	            }
        	            
        	        }
        	    }
        	}, 20);
        	sim.integrator.getEventManager().addListener(bondAction);
    	}
    	
    	
//    	Set<IAtomList> bondPairs = sim.potentialBonded.bondMap.keySet();
//    	Iterator<IAtomList> iterator = bondPairs.iterator();
//    	while(iterator.hasNext()) {
//    	    IAtomList pair = iterator.next();
////            g3dsys.images.Bond b = ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).makeBond((AtomPair)pair,null);
//    	    ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).makeBond((AtomPair)pair,null);
//    	}
//        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).releaseAllBonds();
    	
    }

    public static void main(String[] args) {
        Space sp = Space3D.getInstance();
        RenderMDParam params = new RenderMDParam();
        final RenderMD sim = new RenderMD(sp, params);

        RenderMDGraphic rendermdGraphic = new RenderMDGraphic(sim, sp);
		SimulationGraphic.makeAndDisplayFrame
		        (rendermdGraphic.getPanel(), APP_NAME);
    }
        

}


