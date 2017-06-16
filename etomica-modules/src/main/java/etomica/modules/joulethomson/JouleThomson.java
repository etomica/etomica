/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.joulethomson;

import etomica.action.AtomAction;
import etomica.action.BoxInflate;
import etomica.action.SimulationRestart;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryScrolling;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.*;
import etomica.integrator.IntegratorGear4NPH;
import etomica.integrator.IntegratorGear4NPH.MeterEnthalpy;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.listener.IntegratorListenerAction;
import etomica.math.function.Function;
import etomica.modifier.ModifierFunctionWrapper;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.*;
import etomica.units.systems.MKS;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;

/**
 * Joule Thomson Experiment
 * Ed Maginn [University of Notre Dame], David Kofke [University at Buffalo]
 * Permits study of a system undergoing isothermal or isenthalpic expansion and
 * compression. Isenthalpic expansion is known as  Joule-Thomson process, and 
 * is the basis of many common refrigeration processes.
 */
public class JouleThomson extends SimulationGraphic {

    final static String APP_NAME = "Joule Thomson Experiment";
    final static int REPAINT_INTERVAL = 25;

    DeviceButton restart;
    DeviceBox timeBox;
    JouleThomsonSim sim;


    public JouleThomson(JouleThomsonSim simulation, Space _space) {
        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());
        sim = simulation;
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        final Unit pUnit, dUnit, dadUnit;
        if (sim.getSpace().D() == 2) {
            Unit areaUnit = new MKS().area();
            dUnit = new UnitRatio(Mole.UNIT, 
                                    areaUnit);
            dadUnit = new UnitRatio(new PrefixedUnit(Prefix.DECI, Mole.UNIT),
                                    areaUnit);
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        }
        else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            dadUnit = new UnitRatio(new PrefixedUnit(Prefix.DECI, Mole.UNIT),
                    Liter.UNIT);
            pUnit = Bar.UNIT;
        }

        final Unit tUnit = Kelvin.UNIT;
        final Unit hUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);

	    //colorscheme to color atoms blue to red according to their velocity
	    DeviceSlider scaleSlider = null;
        if(sim.getSpace().D() == 2) 
            getDisplayBox(sim.box).setColorScheme(new ColorSchemeTemperature(100, 500));
        else {
            ColorSchemeByType colorScheme = new ColorSchemeByType(sim);
            colorScheme.setColor(sim.species.getLeafType(), Color.blue);
            getDisplayBox(sim.box).setColorScheme(colorScheme);
        }
        ModifierFunctionWrapper scaleModifier = new ModifierFunctionWrapper(getDisplayBox(sim.box), "scale");
        scaleModifier.setFunction(new Function.Linear(0.01, 0.0));
        scaleSlider = new DeviceSlider(sim.getController());
        scaleSlider.setModifier(scaleModifier);

        // For 2D displays, if the graphic size slider is adjusted, repaint
        // the display
        scaleSlider.getSlider().addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
            	getDisplayBox(sim.box).graphic().repaint();
            }
        });
        scaleSlider.setMinimum(10);
        scaleSlider.setMaximum(100);
        scaleSlider.getSlider().setValue(100);
        scaleSlider.getSlider().setMajorTickSpacing(10);
        scaleSlider.getSlider().setMinorTickSpacing(5);
        scaleSlider.getSlider().setOrientation(1);
        scaleSlider.getSlider().setLabelTable(scaleSlider.getSlider().createStandardLabels(10));
        scaleSlider.setLabel("Scale (%)");
        scaleSlider.setShowBorder(true);
	    //slider for scale of display
	      
        DeviceToggleRadioButtons isothermalToggleButton = new DeviceToggleRadioButtons(new IntegratorGear4NPH.EnsembleToggler(sim.integrator, space.D()), "Ensemble", "NPT", "NPH");

		final DeviceSlider pSlider = new DeviceSlider(sim.getController());
		pSlider.setUnit(pUnit);
        pSlider.setModifier(new ModifierGeneral(sim.integrator,"targetP"));
		pSlider.setLabel("Setpoint Pressure ("+pUnit.symbol()+")");
		pSlider.setMinimum(0);
		pSlider.setMaximum(200);
		pSlider.setNMajor(5);
		pSlider.setShowBorder(true);
		
		final DeviceSlider tSlider = new DeviceSlider(sim.getController());
        tSlider.setModifier(new ModifierGeneral(sim.integrator,"targetT"));
		tSlider.setUnit(tUnit);
		tSlider.setLabel("Setpoint Temperature ("+tUnit.symbol()+")");
		tSlider.setNMajor(4);
		tSlider.setShowBorder(true);

        //meter and display for density
		DataSourceCountTime time = new DataSourceCountTime(sim.integrator);
        final MeterDensity densityMeter = new MeterDensity(sim.getSpace());
        densityMeter.setBox(sim.box);
        AccumulatorHistory densityHistoryAcc = new AccumulatorHistory();
        densityHistoryAcc.setTimeDataSource(time);
        DataPump densityPump = new DataPump(densityMeter, densityHistoryAcc);
        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
        sim.integratorJT.getEventManager().addListener(densityPumpListener);
        densityPumpListener.setInterval(20);

        sim.activityIntegrate.setSleepPeriod(1);

        //set-pressure history
        IEtomicaDataSource targetPressureDataSource = new DataSourceScalar("Set-Pressure",Pressure.DIMENSION) {
            public double getDataAsScalar() {
                return pUnit.toSim(pSlider.getValue());
            }
        };
        AccumulatorHistory targetPressureHistory = new AccumulatorHistory();
        targetPressureHistory.setTimeDataSource(time);
        DataPump pump = new DataPump(targetPressureDataSource, targetPressureHistory);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        sim.integratorJT.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);

        //set-pressure history
        IEtomicaDataSource targetTemperatureDataSource = new DataSourceScalar("Set-Temperature",Temperature.DIMENSION) {
            public double getDataAsScalar() {
                return tUnit.toSim(tSlider.getValue());
            }
        };
        tSlider.setValue(tSlider.getUnit().fromSim(sim.integrator.getTemperature()));
        AccumulatorHistory targetTemperatureHistory = new AccumulatorHistory();
        targetTemperatureHistory.setTimeDataSource(time);
        pump = new DataPump(targetTemperatureDataSource, targetTemperatureHistory);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integratorJT.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);

        //plot of pressure-density-temperature setpoints and averages
        DisplayPlot plotT = new DisplayPlot();
        DisplayPlot plotH = new DisplayPlot();
        DisplayPlot plotP = new DisplayPlot();
        DisplayPlot plotD = new DisplayPlot();
        densityHistoryAcc.addDataSink(plotD.getDataSet().makeDataSink());
        //add actual temperature history
        targetTemperatureHistory.addDataSink(plotT.getDataSet().makeDataSink());
        // add actual pressure history
        targetPressureHistory.addDataSink(plotP.getDataSet().makeDataSink());
        plotD.setLabel("Density");
        plotD.setUnit(new DataTag[]{densityHistoryAcc.getTag()}, dadUnit);
        plotT.setLabel("Temperature");
        plotT.setUnit(new DataTag[]{targetTemperatureHistory.getTag()}, tUnit);
        plotP.setLabel("Pressure");
        plotP.setUnit(new DataTag[]{targetPressureHistory.getTag()}, pUnit);
        plotH.setLabel("Enthalpy");

        MeterPressure meterPressure = new MeterPressure(sim.getSpace());
        meterPressure.setIntegrator(sim.integrator);
        DataFork pressureFork = new DataFork();
        AccumulatorHistory pressureHistory = new AccumulatorHistory(new HistoryScrolling());
        pressureHistory.setTimeDataSource(time);
        pump = new DataPump(meterPressure, pressureFork);
        pressureFork.addDataSink(pressureHistory);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integratorJT.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);
        pressureHistory.addDataSink(plotP.getDataSet().makeDataSink());
        plotP.setUnit(new DataTag[]{pressureHistory.getTag()}, pUnit);
        AccumulatorAverageCollapsing pressureAverage = new AccumulatorAverageCollapsing();
        pressureFork.addDataSink(pressureAverage);
        pressureAverage.setPushInterval(10);
        
        MeterTemperature meterTemperature = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(time);
        pump = new DataPump(meterTemperature, temperatureFork);
        temperatureFork.addDataSink(temperatureHistory);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integratorJT.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);
        temperatureHistory.addDataSink(plotT.getDataSet().makeDataSink());
        plotT.setUnit(new DataTag[]{temperatureHistory.getTag()}, tUnit);
        AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
        temperatureFork.addDataSink(temperatureAverage);
        temperatureAverage.setPushInterval(10);

        MeterEnthalpy meterEnthalpy = new MeterEnthalpy(sim.integrator, space.D());
        DataFork enthalpyFork = new DataFork();
        AccumulatorHistory enthalpyHistory = new AccumulatorHistory();
        enthalpyHistory.setTimeDataSource(time);
        pump = new DataPump(meterEnthalpy, enthalpyFork);
        enthalpyFork.addDataSink(enthalpyHistory);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);
        enthalpyHistory.addDataSink(plotH.getDataSet().makeDataSink());
        plotH.setUnit(new DataTag[]{enthalpyHistory.getTag()}, hUnit);
        AccumulatorAverageCollapsing enthalpyAverage = new AccumulatorAverageCollapsing();
        enthalpyFork.addDataSink(enthalpyAverage);
        enthalpyAverage.setPushInterval(10);
        
        SpeciesChooser speciesChooser = new SpeciesChooser(this, space);
	    speciesChooser.setSpecies("Methane");
	    
        AccumulatorAverageCollapsing densityAverage = new AccumulatorAverageCollapsing();
        pump = new DataPump(densityMeter, densityAverage);
        densityAverage.setPushInterval(10);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integratorJT.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);

        DisplayTextBoxesCAE densityDisplay = new DisplayTextBoxesCAE();
        densityDisplay.setAccumulator(densityAverage);
        densityDisplay.setUnit(dadUnit);

        DisplayTextBoxesCAE temperatureDisplay = new DisplayTextBoxesCAE();
        temperatureDisplay.setAccumulator(temperatureAverage);
        temperatureDisplay.setUnit(tUnit);

        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
        pressureDisplay.setAccumulator(pressureAverage);
        pressureDisplay.setUnit(pUnit);

        DisplayTextBoxesCAE enthalpyDisplay = new DisplayTextBoxesCAE();
        enthalpyDisplay.setAccumulator(enthalpyAverage);
        enthalpyDisplay.setUnit(hUnit);

        timeBox = new DeviceBox();
        timeBox.setLabel("Time Step");
        timeBox.setModifier(new ModifierGeneral(sim.integrator, "timeStep"));

        //panel for sliders
        JPanel sliderPanel = new JPanel(new java.awt.GridLayout(0,1));
        sliderPanel.add(pSlider.graphic(null));
        sliderPanel.add(tSlider.graphic(null));

        //species selector
        JPanel potentialPanel = new JPanel();
        potentialPanel.add(speciesChooser);
	    potentialPanel.setBorder(new TitledBorder("Species Selection"));

	    // Add objects to control panel
	    add(isothermalToggleButton);
	    add(timeBox);
        getPanel().controlPanel.add(potentialPanel, vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);

        // Panel that graphic sits on configuration tab
	    JPanel displayBoxPanel = new JPanel(new BorderLayout());
	    displayBoxPanel.add(getDisplayBox(sim.box).graphic(),BorderLayout.CENTER);

	    // Add widgets specific to 2D application
	    if(sim.getSpace().D() < 3) {
	    	displayBoxPanel.add(scaleSlider.getSlider(),BorderLayout.EAST);
	    }

        // Add objects to tabbed pane (display area)
	    getPanel().tabbedPane.add(getDisplayBox(sim.box).getLabel(), displayBoxPanel);
        add(densityDisplay);
        add(temperatureDisplay);
        add(pressureDisplay);
        add(enthalpyDisplay);
	    getPanel().tabbedPane.add(plotD.getLabel(), plotD.graphic());
	    getPanel().tabbedPane.add(plotT.getLabel(), plotT.graphic());
        getPanel().tabbedPane.add(plotP.getLabel(), plotP.graphic());
        getPanel().tabbedPane.add(plotH.getLabel(), plotH.graphic());

	    getPanel().toolbar.addContributor("Colin Tedlock");

	    getController().getReinitButton().setPostAction(getPaintAction(sim.box));
    }

    public static void main(String[] args) {
        int dim = 2;
        if (args.length > 0) {
            try {
                dim = Integer.parseInt(args[0]);
            } catch (NumberFormatException ex) {
                throw new NumberFormatException(args[0] + " is not a number!");
            }
        }

        Space sp = Space.getInstance(dim);
        SimulationGraphic simPanel = new JouleThomson(new JouleThomsonSim(sp), sp);
        SimulationGraphic.makeAndDisplayFrame(simPanel.getPanel(), APP_NAME);
    } //end of main

    //inner class the defines a drop-down menu to select LJ parameters to mimic
    //several real substances
    //the field "names" must be static to pass to super, and that requires this inner class to be static
    public static class SpeciesChooser extends JComboBox implements java.awt.event.ItemListener {

        public static final String[] names = new String[]
            {"Methane", "Argon", "Krypton", "Nitrogen", "Ideal gas"};
        private final Space space;
        double[] sigma   = new double[] {
                4.010,
                3.499,
                3.846,
                3.694,
                                         2.};
        double[] epsilon = new double[] {
                Kelvin.UNIT.toSim(142.87),
                                         Kelvin.UNIT.toSim(118.13),
                Kelvin.UNIT.toSim(162.74),
                                         Kelvin.UNIT.toSim(96.26),
                                         Kelvin.UNIT.toSim(100.)};
        double[] mass    = new double[] {
                16.0,
                39.95,
                83.80,
                28.0,
                                         1.0};
        JouleThomsonSim sim;
        JouleThomson simGraphic;
        double currentMass = mass[0];
        double currentEps = epsilon[0];
        double currentSig = sigma[0];
        AtomAction updateMass = new AtomAction() {
            public void actionPerformed(IAtom a) {((ElementSimple)a.getType().getElement()).setMass(currentMass);}
        };
        SimulationRestart simRestart;

        SpeciesChooser(JouleThomson simGraphic, Space _space) {
            super(names);
            this.space = _space;
            this.simGraphic = simGraphic;
            setPreferredSize(new java.awt.Dimension(150,30));
            setOpaque(false);
            addItemListener(this);
            sim = simGraphic.sim;
            currentEps = sim.potential.getEpsilon();
            currentSig = sim.potential.getSigma();
            sigma[0] = currentSig;
            epsilon[0] = currentEps;
            mass[0] = currentMass;
            simRestart = new SimulationRestart(sim, space, sim.getController());
        }

        public void itemStateChanged(java.awt.event.ItemEvent evt) {
            if(evt.getStateChange() == java.awt.event.ItemEvent.DESELECTED) return;
            String speciesName = (String)evt.getItemSelectable().getSelectedObjects()[0];
            setSpecies(speciesName);

            // Display the new box
            simGraphic.getDisplayBox(sim.box).repaint();
        }

        public void setSpecies(String speciesName) {
            for(int i=0; i<names.length; i++) {
                if(speciesName.equals(names[i])) {
                    double tStep = sim.integrator.getTimeStep();
//                    tStep *= sigma[i]/currentSig*Math.sqrt(mass[i]*currentEps/(currentMass*epsilon[i]));
                    tStep *= sigma[i]/currentSig*Math.sqrt(mass[i]/(currentMass));
                    sim.integrator.setTimeStep(tStep);
                    sim.integratorNVE.setTimeStep(tStep);
                    currentMass = mass[i];
                    ((ElementSimple)sim.species.getLeafType().getElement()).setMass(currentMass);
                    currentEps = epsilon[i];
                    currentSig = sigma[i];
                    ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), sigma[i]);
                    sim.potential.setEpsilon(epsilon[i]);
                    sim.potential.setSigma(sigma[i]);
                }
            }
    	    double targetDensity = 0.5;
    	    if(currentSig > 0.5) {
    	        BoxInflate inflater = new BoxInflate(sim.box, space);
    	        inflater.setTargetDensity(targetDensity/Math.pow(currentSig,sim.getSpace().D()));
    	        inflater.actionPerformed();
    	        double size = currentSig*Math.pow(sim.box.getLeafList().getAtomCount()/targetDensity,1.0/sim.getSpace().D());
                Vector v = sim.getSpace().makeVector();
                v.E(size);
                sim.box.getBoundary().setBoxSize(v);
                SpaceLattice lattice;
                if (space.D() == 2) {
                    lattice = new LatticeOrthorhombicHexagonal(space);
                }
                else {
                    lattice = new LatticeCubicFcc(space);
                }
                Configuration config = new ConfigurationLattice(lattice, space);
                config.initializeCoordinates(sim.box);
            }
            if(speciesName.equals("Ideal gas")) {
                if (sim.integrator.getPotentialMaster().getPotentials().length > 0) {
                    sim.integrator.getPotentialMaster().removePotential(sim.potential);
                }
            }
            else {
                if (sim.integrator.getPotentialMaster().getPotentials().length == 0) {
                    sim.integrator.getPotentialMaster().addPotential(sim.potential, new AtomType[]{sim.species.getLeafType(), sim.species.getLeafType()});
                }
            }
            simRestart.actionPerformed();
        }
    }

    public static class Applet extends javax.swing.JApplet {

	    public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space.getInstance(3);
	        SimulationGraphic simPanel = new JouleThomson(new JouleThomsonSim(sp), sp);
		    getContentPane().add(simPanel.getPanel());
	    }
    }//end of Applet

}
