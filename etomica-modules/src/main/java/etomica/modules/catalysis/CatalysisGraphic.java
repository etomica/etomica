/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

 import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
 import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Fraction;
import etomica.units.Kelvin;
import etomica.units.Length;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.Volume;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.history.HistoryCollapsingAverage;

/**
 * Catalysis graphical app.
 * Design by Ken Benjamin
 * 
 * @author Andrew Schultz
 */
public class CatalysisGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Catalysis";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    protected Catalysis sim;

    public CatalysisGraphic(final Catalysis simulation, Space _space, boolean showParams) {

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
        double sigmaS = 2*sim.potentialCS.getCoreDiameter() - sim.potentialCC.getCoreDiameter();
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesSurface.getLeafType(), sigmaS);
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesO.getLeafType(), sim.potentialOO.getCoreDiameter());
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesC.getLeafType(), sim.potentialCC.getCoreDiameter());

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump= new DataPump(meterCycles,displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");

        DisplayTextBox displayVolume = new DisplayTextBox();
        displayVolume.putDataInfo(new DataInfoDouble("Volume", Volume.DIMENSION));
        DataDouble dataV = new DataDouble();
        dataV.x = sim.box.getBoundary().volume();
        displayVolume.putData(dataV);

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setUnit(Kelvin.UNIT);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(1000);
        tempSlider.setSliderMajorValues(4);
        tempSlider.setUnit(tUnit);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        //display of box, timer
        ColorSchemeRadical colorScheme = new ColorSchemeRadical(sim, sim.interactionTracker.getAgentManager());
        colorScheme.setColor(sim.speciesO.getLeafType(),java.awt.Color.RED);
        colorScheme.setColor(sim.speciesC.getLeafType(),java.awt.Color.BLUE);
        colorScheme.setColor(sim.speciesSurface.getLeafType(),Color.GRAY);
        colorScheme.setRadicalColor(sim.speciesO.getLeafType(),Color.PINK);
        colorScheme.setRadicalColor(sim.speciesC.getLeafType(),Color.CYAN);
        colorScheme.setFullBondColor(sim.speciesC.getLeafType(),Color.YELLOW);
        getDisplayBox(sim.box).setColorScheme(colorScheme);
		
		MeterTemperature thermometer = new MeterTemperature(sim, sim.box, space.D());
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        final DataPumpListener temperaturePump = new DataPumpListener(thermometer,temperatureHistory);
        sim.integrator.getEventManager().addListener(temperaturePump);
		dataStreamPumps.add(temperaturePump);
        DisplayPlot temperatureHistoryPlot = new DisplayPlot();
        temperatureHistory.setDataSink(temperatureHistoryPlot.getDataSet().makeDataSink());
        temperatureHistoryPlot.setLabel("Temperature");
        temperatureHistoryPlot.setDoLegend(false);
        temperatureHistoryPlot.setUnit(Kelvin.UNIT);
        temperatureHistoryPlot.getPlot().setYLabel("Temperature (K)");
        
        Unit dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);

	    final MeterDensityCO meterDensityCO = new MeterDensityCO(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
	    DataFork densityCOFork = new DataFork();
        final DataPumpListener densityCOPump = new DataPumpListener(meterDensityCO, densityCOFork, 100);
        sim.integrator.getEventManager().addListener(densityCOPump);
        dataStreamPumps.add(densityCOPump);
        AccumulatorHistory densityCOHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
        densityCOFork.addDataSink(densityCOHistory);
        densityCOHistory.setTimeDataSource(meterCycles);
        DisplayTextBox densityCOBox = new DisplayTextBox();
        densityCOBox.setUnit(dUnit);
        densityCOFork.addDataSink(densityCOBox);
        densityCOBox.setLabel("CO Density");
        
        final MeterDensityO2 meterDensityO2 = new MeterDensityO2(sim.box, sim.speciesO, sim.interactionTracker.getAgentManager());
        DataFork densityO2Fork = new DataFork();
        final DataPumpListener densityO2Pump = new DataPumpListener(meterDensityO2, densityO2Fork, 100);
        sim.integrator.getEventManager().addListener(densityO2Pump);
        dataStreamPumps.add(densityO2Pump);
        AccumulatorHistory densityO2History = new AccumulatorHistory(new HistoryCollapsingDiscard());
        densityO2Fork.addDataSink(densityO2History);
        densityO2History.setTimeDataSource(meterCycles);
        DisplayTextBox densityO2Box = new DisplayTextBox();
        densityO2Box.setUnit(dUnit);
        densityO2Fork.addDataSink(densityO2Box);
        densityO2Box.setLabel("O2 Density");

        final MeterDensityCO2 meterDensityCO2 = new MeterDensityCO2(sim.box, sim.speciesC, sim.interactionTracker.getAgentManager());
        DataFork densityCO2Fork = new DataFork();
        final DataPumpListener densityCO2Pump = new DataPumpListener(meterDensityCO2, densityCO2Fork, 100);
        sim.integrator.getEventManager().addListener(densityCO2Pump);
        dataStreamPumps.add(densityCO2Pump);
        AccumulatorHistory densityCO2History = new AccumulatorHistory(new HistoryCollapsingDiscard());
        densityCO2Fork.addDataSink(densityCO2History);
        densityCO2History.setTimeDataSource(meterCycles);
        DisplayTextBox densityCO2Box = new DisplayTextBox();
        densityCO2Box.setUnit(dUnit);
        densityCO2Fork.addDataSink(densityCO2Box);
        densityCO2Box.setLabel("CO2 Density");

        DisplayPlot densityHistoryPlot = new DisplayPlot();
        densityCOHistory.setDataSink(densityHistoryPlot.getDataSet().makeDataSink());
        densityO2History.setDataSink(densityHistoryPlot.getDataSet().makeDataSink());
        densityCO2History.setDataSink(densityHistoryPlot.getDataSet().makeDataSink());
        densityHistoryPlot.setLabel("Density");
        densityHistoryPlot.setLegend(new DataTag[]{meterDensityCO.getTag()}, "CO");
        densityHistoryPlot.setLegend(new DataTag[]{meterDensityO2.getTag()}, "O2");
        densityHistoryPlot.setLegend(new DataTag[]{meterDensityCO2.getTag()}, "CO2");
        densityHistoryPlot.setUnit(dUnit);
        densityHistoryPlot.getPlot().setYLabel("Density (mol/L)");

        final DeviceSlider nSliderCO = new DeviceSlider(sim.getController());
        nSliderCO.setMinimum(0);
        nSliderCO.setMaximum(500);
        nSliderCO.setShowBorder(true);
        nSliderCO.setShowValues(true);
        nSliderCO.setEditValues(true);
        nSliderCO.setModifier(new ModifierGeneral(sim.config, "numCO"));
        nSliderCO.setLabel("Number of CO");
        nSliderCO.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.interactionTracker.reset();
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
        nSliderO2.setEditValues(true);
        nSliderO2.setModifier(new ModifierGeneral(sim.config, "numO2"));
        nSliderO2.setLabel("Number of O2");
        nSliderO2.setPostAction(new IAction() {
            public void actionPerformed() {
                sim.interactionTracker.reset();
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


		tempSlider.setSliderPostAction(resetDataAction);
        tempSlider.setRadioGroupPostAction(resetDataAction);

        final IAction resetDisplayDataAction = new IAction() {
            public void actionPerformed() {
                temperaturePump.actionPerformed();

                displayCycles.putData(meterCycles.getData());
                displayCycles.repaint();
            }
        };
        final IAction resetAction = new IAction() {
        	public void actionPerformed() {
        	    sim.interactionTracker.reset();
        	    
        	    sim.integrator.reset();

        	    resetDisplayDataAction.actionPerformed();

        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };
        
        ((SimulationRestart)getController().getReinitButton().getAction()).setConfiguration(sim.config);

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetDisplayDataAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        
        if (showParams) {
            JTabbedPane controlsTabs = new JTabbedPane();
            JPanel mainControls = new JPanel(new GridBagLayout());
            controlsTabs.add(mainControls, "State");
            
            mainControls.add(statePanel, vertGBC);
            mainControls.add(delaySlider.graphic(), vertGBC);
            getPanel().controlPanel.add(controlsTabs, vertGBC);

            // OO
            {
            JPanel controlsOO = new JPanel(new GridBagLayout());
            controlsTabs.add(controlsOO, "OO");
            DeviceBox sigmaOBox = new DeviceBox();
            sigmaOBox.setController(sim.getController());
            sigmaOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Length.DIMENSION;
                }

                public String getLabel() {
                    return "sigma";
                }

                public double getValue() {
                    return sim.potentialOO.getCoreDiameter();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    double oldValue = sim.potentialOO.getCoreDiameter();
                    sim.potentialOO.setCoreDiameter(newValue);
                    ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesO.getLeafType(), newValue);
                    sim.potentialCO.setCoreDiameter(0.5*(sim.potentialCC.getCoreDiameter()+newValue));
                    double sigmaOS = sim.potentialOS.getCoreDiameter();
                    sigmaOS += 0.5*(newValue-oldValue);
                    sim.potentialOS.setCoreDiameter(sigmaOS);
                    sim.config.initializeCoordinates(sim.box);
                    resetAction.actionPerformed();
                }
            });
            controlsOO.add(sigmaOBox.graphic(), vertGBC);
            
            DeviceBox epsOBox = new DeviceBox();
            epsOBox.setController(sim.getController());
            epsOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "epsilon";
                }

                public double getValue() {
                    return sim.potentialOO.getEpsilon();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialOO.setEpsilon(newValue);
                    sim.potentialCO.setEpsilon(Math.sqrt(sim.potentialCC.getEpsilon()*newValue));
                }
            });
            epsOBox.setUnit(Kelvin.UNIT);
            epsOBox.doUpdate();
            controlsOO.add(epsOBox.graphic(), vertGBC);
            
            DeviceBox lambdaOBox = new DeviceBox();
            lambdaOBox.setController(sim.getController());
            lambdaOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Fraction.DIMENSION;
                }

                public String getLabel() {
                    return "lambda";
                }

                public double getValue() {
                    return sim.potentialOO.getLambda();
                }

                public void setValue(double newValue) {
                    if (newValue <= 1) throw new RuntimeException("value must be greater than 1");
                    sim.potentialOO.setLambda(newValue);
                }
            });
            lambdaOBox.doUpdate();
            controlsOO.add(lambdaOBox.graphic(), vertGBC);

            DeviceBox nSitesOBox = new DeviceBox();
            nSitesOBox.setController(sim.getController());
            nSitesOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Quantity.DIMENSION;
                }

                public String getLabel() {
                    return "# of sites";
                }

                public double getValue() {
                    return sim.potentialOO.getNumSurfaceSites();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialOO.setNumSurfaceSites((int)newValue);
                    sim.potentialOS.setMinRadicalSites((int)newValue);
                    sim.config.initializeCoordinates(sim.box);
                    resetAction.actionPerformed();
                }
            });
            nSitesOBox.setInteger(true);
            nSitesOBox.setPrecision(0);
            nSitesOBox.doUpdate();
            controlsOO.add(nSitesOBox.graphic(), vertGBC);

            DeviceBox epsBarrierOBox = new DeviceBox();
            epsBarrierOBox.setController(sim.getController());
            epsBarrierOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "barrier";
                }

                public double getValue() {
                    return sim.potentialOO.getBarrier();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialOO.setBarrier(newValue);
                }
            });
            epsBarrierOBox.setUnit(Kelvin.UNIT);
            epsBarrierOBox.doUpdate();
            controlsOO.add(epsBarrierOBox.graphic(), vertGBC);

            DeviceBox epsBondingOBox = new DeviceBox();
            epsBondingOBox.setController(sim.getController());
            epsBondingOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "eps bonding";
                }

                public double getValue() {
                    return sim.potentialOO.getEpsilonBonding();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialOO.setEpsilonBonding(newValue);
                }
            });
            epsBondingOBox.setUnit(Kelvin.UNIT);
            epsBondingOBox.doUpdate();
            controlsOO.add(epsBondingOBox.graphic(), vertGBC);

            DeviceBox epsOSBox = new DeviceBox();
            epsOSBox.setController(sim.getController());
            epsOSBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "surface epsilon";
                }

                public double getValue() {
                    return sim.potentialOS.getEpsilon();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialOS.setEpsilon(newValue);
                }
            });
            epsOSBox.setUnit(Kelvin.UNIT);
            epsOSBox.doUpdate();
            controlsOO.add(epsOSBox.graphic(), vertGBC);
            }
            
            // CC
            {
            JPanel controlsCC = new JPanel(new GridBagLayout());
            controlsTabs.add(controlsCC, "CC");
            DeviceBox sigmaCBox = new DeviceBox();
            sigmaCBox.setController(sim.getController());
            sigmaCBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Length.DIMENSION;
                }

                public String getLabel() {
                    return "sigma";
                }

                public double getValue() {
                    return sim.potentialCC.getCoreDiameter();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    double oldValue = sim.potentialCC.getCoreDiameter();
                    sim.potentialCC.setCoreDiameter(newValue);
                    ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesC.getLeafType(), newValue);
                    sim.potentialCO.setCoreDiameter(0.5*(sim.potentialOO.getCoreDiameter()+newValue));
                    double sigmaCS = sim.potentialCS.getCoreDiameter();
                    sigmaCS += 0.5*(newValue-oldValue);
                    sim.potentialCS.setCoreDiameter(sigmaCS);
                    sim.config.initializeCoordinates(sim.box);
                    resetAction.actionPerformed();
                }
            });
            controlsCC.add(sigmaCBox.graphic(), vertGBC);
            
            DeviceBox epsCBox = new DeviceBox();
            epsCBox.setController(sim.getController());
            epsCBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "epsilon";
                }

                public double getValue() {
                    return sim.potentialCC.getEpsilon();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialCC.setEpsilon(newValue);
                    sim.potentialCO.setEpsilon(Math.sqrt(sim.potentialOO.getEpsilon()*newValue));
                }
            });
            epsCBox.setUnit(Kelvin.UNIT);
            epsCBox.doUpdate();
            controlsCC.add(epsCBox.graphic(), vertGBC);

            DeviceBox lambdaCBox = new DeviceBox();
            lambdaCBox.setController(sim.getController());
            lambdaCBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Fraction.DIMENSION;
                }

                public String getLabel() {
                    return "lambda";
                }

                public double getValue() {
                    return sim.potentialCC.getLambda();
                }

                public void setValue(double newValue) {
                    if (newValue <= 1) throw new RuntimeException("value must be greater than 1");
                    sim.potentialCC.setLambda(newValue);
                }
            });
            lambdaCBox.doUpdate();
            controlsCC.add(lambdaCBox.graphic(), vertGBC);

            DeviceBox epsCSBox = new DeviceBox();
            epsCSBox.setController(sim.getController());
            epsCSBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "surface epsilon";
                }

                public double getValue() {
                    return sim.potentialCS.getEpsilon();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialCS.setEpsilon(newValue);
                }
            });
            epsCSBox.setUnit(Kelvin.UNIT);
            epsCSBox.doUpdate();
            controlsCC.add(epsCSBox.graphic(), vertGBC);
            }


            // CO
            {
            JPanel controlsCO = new JPanel(new GridBagLayout());
            controlsTabs.add(controlsCO, "CO");
            
            DeviceBox lambdaCOBox = new DeviceBox();
            lambdaCOBox.setController(sim.getController());
            lambdaCOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Fraction.DIMENSION;
                }

                public String getLabel() {
                    return "lambda";
                }

                public double getValue() {
                    return sim.potentialCO.getLambda();
                }

                public void setValue(double newValue) {
                    if (newValue <= 1) throw new RuntimeException("value must be greater than 1");
                    sim.potentialCO.setLambda(newValue);
                }
            });
            lambdaCOBox.doUpdate();
            controlsCO.add(lambdaCOBox.graphic(), vertGBC);

            DeviceBox nSitesCOBox = new DeviceBox();
            nSitesCOBox.setController(sim.getController());
            nSitesCOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Quantity.DIMENSION;
                }

                public String getLabel() {
                    return "# of sites";
                }

                public double getValue() {
                    return sim.reactionManagerCO.getnReactCO();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.reactionManagerCO.setnReactCO((int)newValue);
                    sim.potentialCS.setMinRadicalSites((int)newValue);
                    sim.config.initializeCoordinates(sim.box);
                    resetAction.actionPerformed();
                }
            });
            nSitesCOBox.setInteger(true);
            nSitesCOBox.setPrecision(0);
            nSitesCOBox.doUpdate();
            controlsCO.add(nSitesCOBox.graphic(), vertGBC);

            DeviceBox uRadCOBox = new DeviceBox();
            uRadCOBox.setController(sim.getController());
            uRadCOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "U rad";
                }

                public double getValue() {
                    return sim.reactionManagerCO.getuReactCO();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.reactionManagerCO.setuReactCO(newValue);
                }
            });
            uRadCOBox.setUnit(Kelvin.UNIT);
            controlsCO.add(uRadCOBox.graphic(), vertGBC);

            DeviceBox uRadCORevBox = new DeviceBox();
            uRadCORevBox.setController(sim.getController());
            uRadCORevBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "U derad";
                }

                public double getValue() {
                    return sim.reactionManagerCO.getuReactCORev();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.reactionManagerCO.setuReactCORev(newValue);
                }
            });
            uRadCORevBox.setUnit(Kelvin.UNIT);
            controlsCO.add(uRadCORevBox.graphic(), vertGBC);

            DeviceBox epsBarrierCOBox = new DeviceBox();
            epsBarrierCOBox.setController(sim.getController());
            epsBarrierCOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "barrier";
                }

                public double getValue() {
                    return sim.potentialCO.getBarrier();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialCO.setBarrier(newValue);
                }
            });
            epsBarrierCOBox.setUnit(Kelvin.UNIT);
            epsBarrierCOBox.doUpdate();
            controlsCO.add(epsBarrierCOBox.graphic(), vertGBC);

            DeviceBox epsBondingCOBox = new DeviceBox();
            epsBondingCOBox.setController(sim.getController());
            epsBondingCOBox.setModifier(new Modifier() {

                public Dimension getDimension() {
                    return Energy.DIMENSION;
                }

                public String getLabel() {
                    return "eps bonding";
                }

                public double getValue() {
                    return sim.potentialCO.getEpsilonBonding();
                }

                public void setValue(double newValue) {
                    if (newValue <= 0) throw new RuntimeException("value must be positive");
                    sim.potentialCO.setEpsilonBonding(newValue);
                }
            });
            epsBondingCOBox.setUnit(Kelvin.UNIT);
            epsBondingCOBox.doUpdate();
            controlsCO.add(epsBondingCOBox.graphic(), vertGBC);
            }

        }
        else {
            getPanel().controlPanel.add(statePanel, vertGBC);
            getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
            
            DeviceSelector surfaceSelector = new DeviceSelector(sim.getController());
            surfaceSelector.setLabel("Catalyst type");
            surfaceSelector.addOption("Catalyst A", new IAction() {
                public void actionPerformed() {
                    sim.potentialCS.setEpsilon(Kelvin.UNIT.toSim(2000));
                    sim.potentialOS.setEpsilon(Kelvin.UNIT.toSim(500));
                    sim.reactionManagerCO.setuReactCO(Kelvin.UNIT.toSim(80));
                    sim.potentialCO.setBarrier(Kelvin.UNIT.toSim(400));
                    sim.potentialOO.setBarrier(Kelvin.UNIT.toSim(200));
                }
            });
            surfaceSelector.addOption("Catalyst B", new IAction() {
                public void actionPerformed() {
                    sim.potentialCS.setEpsilon(Kelvin.UNIT.toSim(3000));
                    sim.potentialOS.setEpsilon(Kelvin.UNIT.toSim(700));
                    sim.reactionManagerCO.setuReactCO(Kelvin.UNIT.toSim(300));
                    sim.potentialCO.setBarrier(Kelvin.UNIT.toSim(4000));
                    sim.potentialOO.setBarrier(Kelvin.UNIT.toSim(400));
                }
            });

            getPanel().controlPanel.add(surfaceSelector.graphic(), vertGBC);
        }

        add(displayCycles);
        add(displayVolume);
        add(densityCOBox);
        add(densityO2Box);
        add(densityCO2Box);
        add(temperatureHistoryPlot);
        add(densityHistoryPlot);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        boolean showParams = false;
        int nCellsZ = 20;
        if(args.length != 0) {
            try {
                showParams = Integer.parseInt(args[0]) != 0;
            } catch(NumberFormatException e) {}
            
            if(args.length > 1) {
                try {
                    nCellsZ = Integer.parseInt(args[1]);
                } catch(NumberFormatException e) {}
            }
        }

        CatalysisGraphic swmdGraphic = new CatalysisGraphic(new Catalysis(space, nCellsZ), space, showParams);
		SimulationGraphic.makeAndDisplayFrame
		        (swmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            String paramStr = getParameter("showParams");
            boolean showParams = false;
            if (paramStr != null) {
                showParams = Integer.valueOf(paramStr) != 0;
            }
            int nCellsZ = 20;
            String nCellsZStr = getParameter("nCellsZ");
            if (nCellsZStr != null) {
                nCellsZ = Integer.valueOf(nCellsZStr);
            }
            Space sp = Space3D.getInstance();
            CatalysisGraphic swmdGraphic = new CatalysisGraphic(new Catalysis(sp, nCellsZ), sp, showParams);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
