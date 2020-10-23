/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.swmd;

import etomica.action.BoxImposePbc;
import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.*;
import etomica.data.types.DataDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.math.DoubleRange;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2HardSphere;
import etomica.potential.P2Ideal;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.statmech.MaxwellBoltzmann;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.units.systems.MKS;
import etomica.util.Constants.CompassDirection;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

public class SwmdGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Square-Well Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    public ItemListener potentialChooserListener;
    public JComboBox potentialChooser;
    protected P2HardSphere potentialHS;
    protected P2SquareWell potentialSW;
    protected P2Ideal potentialIdeal;
    public DeviceBox sigBox, epsBox, lamBox, massBox;
    public double lambda, epsilon, mass, sigma;
    protected Unit tUnit, eUnit, dUnit, pUnit, mUnit;
    protected Swmd sim;
    
    private boolean showConfig = true;

    public SwmdGraphic(final Swmd simulation, Space _space) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

        this.sim = simulation;

        tUnit = Kelvin.UNIT;

        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
        mUnit = new UnitRatio(Gram.UNIT, Mole.UNIT);
        lambda = 2.0;
        epsilon = eUnit.toSim(space.D() == 3 ? 1000 : 1500);
        mass = space.D() == 3 ? 131 : 40;
        sigma = 4.0;
        if (sim.getSpace().D() == 2) {
            dUnit = new UnitRatio(Mole.UNIT,
                    new MKS().area());
            Unit[] units = new Unit[]{Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[]{1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        } else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;
        }

        if (sim.getSpace().D() == 2) {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(400 / sim.box.getBoundary().getBoxSize().getX(1)));
        } else {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(40 / sim.box.getBoundary().getBoxSize().getX(1)));
        }

        sim.getController().setSleepPeriod(0);
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

        //combo box to select potentials
        final String idealGas = "Ideal gas";
        final String repulsionOnly = "Repulsion only";
        final String repulsionAttraction = "Repulsion and attraction";
        potentialChooser = new JComboBox(new String[]{
                idealGas, repulsionOnly, repulsionAttraction});

        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        lamBox = new DeviceBox();
        massBox = new DeviceBox();
        // Unselectable because "Ideal gas" is selected initially
        potentialChooser.setSelectedIndex(0);
        epsBox.setEditable(false);
        lamBox.setEditable(false);

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump = new DataPump(meterCycles, displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setUnit(tUnit);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(1500.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setAdiabatic();

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(potentialChooser, vertGBC);
        JPanel parameterPanel = new JPanel(new GridLayout(0, 1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(lamBox.graphic());
        parameterPanel.add(massBox.graphic());
        potentialPanel.add(parameterPanel, vertGBC);

        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");

        potentialSW = new P2SquareWell(sim.getSpace(), sigma, lambda, epsilon, true);
        potentialHS = new P2HardSphere(sim.getSpace(), sigma, true);
        potentialIdeal = new P2Ideal(sim.getSpace());

        if (potentialChooserListener != null) potentialChooser.removeItemListener(potentialChooserListener);

        potentialChooserListener = new ItemListener() {
            public void itemStateChanged(ItemEvent evt) {
                if (evt.getStateChange() == ItemEvent.DESELECTED) return;
                setPotential((String) evt.getItem());
                if (evt.getItem() == idealGas ||
                        evt.getItem() == repulsionOnly) {
                    epsBox.setEditable(false);
                    lamBox.setEditable(false);
                } else {
                    epsBox.setEditable(true);
                    lamBox.setEditable(true);
                }
            }
        };
        potentialChooser.addItemListener(potentialChooserListener);
        setPotential((String) potentialChooser.getSelectedItem());

        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter();
        sigModifier.setValue(sigma);
        ModifierGeneral epsModifier = new ModifierGeneral(potentialSW, "epsilon");
        ModifierGeneral lamModifier = new ModifierGeneral(potentialSW, "lambda");
        ModifierGeneral massModifier = new ModifierGeneral(sim.species.getLeafType().getElement(), "mass");
        sigBox.setModifier(sigModifier);
        sigBox.setLabel("Core Diameter (" + Angstrom.UNIT.symbol() + ")");
        epsBox.setUnit(eUnit);
        epsBox.setModifier(epsModifier);
        lamBox.setModifier(lamModifier);
        massBox.setModifier(massModifier);
        massBox.setUnit(mUnit);
        sigBox.setController(sim.getController());
        epsBox.setController(sim.getController());
        lamBox.setController(sim.getController());
        massBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getLeafType(), Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());

        //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListener = new IntegratorListenerAction(rdfMeter);
        sim.integrator.getEventManager().addListener(rdfMeterListener);
        rdfMeterListener.setInterval(10);
        rdfMeter.getXDataSource().setXMax(12.0);
        rdfMeter.setBox(sim.box);
        DisplayPlotXChart rdfPlot = new DisplayPlotXChart();
        DataPump rdfPump = new DataPump(rdfMeter, rdfPlot.getDataSet().makeDataSink());
        IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
        sim.integrator.getEventManager().addListener(rdfPumpListener);
        rdfPumpListener.setInterval(10);
        getController().getResetAveragesButton().setPostAction(new IAction() {
            public void actionPerformed() {
                rdfMeter.reset();
            }
        });

        rdfPlot.setDoLegend(false);
        rdfPlot.getPlot().setTitle("Radial Distribution Function");
        rdfPlot.setLabel("RDF");

        //velocity distribution
        double vMin = 0;
        double vMax = 10;
        DataSourceRmsVelocity meterVelocity = new DataSourceRmsVelocity(new HistogramSimple(100, new DoubleRange(0, vMax)));
        meterVelocity.setIterator(new AtomIteratorLeafAtoms(sim.box));
        AccumulatorAverage rmsAverage = new AccumulatorAverageFixed(10);
        DataPump velocityPump = new DataPump(meterVelocity, rmsAverage);
        IntegratorListenerAction velocityPumpListener = new IntegratorListenerAction(velocityPump);
        sim.integrator.getEventManager().addListener(velocityPumpListener);
        velocityPumpListener.setInterval(10);
        rmsAverage.setPushInterval(1);
        dataStreamPumps.add(velocityPump);

        final DisplayPlotXChart vPlot = new DisplayPlotXChart();
        rmsAverage.addDataSink(vPlot.getDataSet().makeDataSink(), new StatType[]{rmsAverage.AVERAGE});
        vPlot.setLegend(new DataTag[]{meterVelocity.getTag()}, "measured");
        vPlot.setDoLegend(false);
        vPlot.getPlot().setTitle("Velocity Distribution");
        vPlot.setDoLegend(true);
        vPlot.setLabel("Velocity");

        final MaxwellBoltzmann.Distribution mbDistribution = new MaxwellBoltzmann.Distribution(sim.getSpace(), sim.integrator.getTemperature(), sim.species.getLeafType().getMass());
        final DataSourceFunction mbSource = new DataSourceFunction("Maxwell Boltzmann",
                Null.DIMENSION, mbDistribution, 100, "Speed (A/ps)", new DimensionRatio(Length.DIMENSION, Time.DIMENSION));
        DataSourceUniform mbX = mbSource.getXSource();
        mbX.setTypeMax(LimitType.HALF_STEP);
        mbX.setTypeMin(LimitType.HALF_STEP);
        mbX.setNValues(meterVelocity.getDataInfo().getLength());
        mbX.setXMin(vMin);
        mbX.setXMax(vMax);
        mbSource.update();
        DataPump mbPump = new DataPump(mbSource, vPlot.getDataSet().makeDataSink());
        IntegratorListenerAction mbPumpListener = new IntegratorListenerAction(mbPump);
        sim.integrator.getEventManager().addListener(mbPumpListener);
        mbPumpListener.setInterval(10);

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

        MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer, temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(1);
        final AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
        temperatureAverage.setPushInterval(20);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
        temperatureFork.setDataSinks(new IDataSink[]{temperatureAverage, temperatureHistory});
        final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
        tBox.setAccumulator(temperatureAverage);
        dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
        tBox.setLabel("Measured Temperature (K)");
        tBox.setLabelPosition(CompassDirection.NORTH);

        // Number density box
        MeterDensity densityMeter = new MeterDensity(sim.box);
        final DisplayTextBox densityBox = new DisplayTextBox();
        densityBox.setUnit(dUnit);
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
        sim.integrator.getEventManager().addListener(densityPumpListener);
        densityPumpListener.setInterval(1);
        dataStreamPumps.add(densityPump);
        densityBox.setLabel("Density");

        MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        final AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        final DataSinkExcludeOverlap eExcludeOverlap = new DataSinkExcludeOverlap();
        eExcludeOverlap.setDataSink(energyHistory);
        final DataPump energyPump = new DataPump(eMeter, eExcludeOverlap);
        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        sim.integrator.getEventManager().addListener(energyPumpListener);
        dataStreamPumps.add(energyPump);

        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        final DataSinkExcludeOverlap peExcludeOverlap = new DataSinkExcludeOverlap();
        peExcludeOverlap.setDataSink(peFork);
        final DataPump pePump = new DataPump(peMeter, peExcludeOverlap);
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        sim.integrator.getEventManager().addListener(pePumpListener);
        dataStreamPumps.add(pePump);

        MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        final DataSinkExcludeOverlap keExcludeOverlap = new DataSinkExcludeOverlap();
        keExcludeOverlap.setDataSink(keHistory);
        final DataPump kePump = new DataPump(keMeter, keExcludeOverlap);
        IntegratorListenerAction kePumpListener = new IntegratorListenerAction(kePump);
        sim.integrator.getEventManager().addListener(kePumpListener);
        dataStreamPumps.add(kePump);
        int numAtoms = sim.box.getLeafList().size();
        energyPumpListener.setInterval(numAtoms > 120 ? 1 : 120 / numAtoms);
        kePumpListener.setInterval(numAtoms > 120 ? 1 : 120 / numAtoms);
        pePumpListener.setInterval(numAtoms > 120 ? 1 : 120 / numAtoms);

        final DisplayPlotXChart ePlot = new DisplayPlotXChart();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History (J/mol)");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy");
        ePlot.setUnit(eUnit);
        ePlot.setXUnit(Picosecond.UNIT);

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        final DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator);
        sim.integrator.getEventManager().addListener(pPump);
        pAccumulator.setPushInterval(50);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setLabel(sim.getSpace().D() == 3 ? "Pressure (bar)" : "Pressure (bar-nm)");
        pDisplay.setAccumulator(pAccumulator);
        pDisplay.setUnit(pUnit);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (J/mol)");
        peDisplay.setUnit(eUnit);

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setResetAction(new SimulationRestart(sim));
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(sim.getSpace().D() == 3 ? 500 : 168);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        ChangeListener nListener = new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                final int n = (int) nSlider.getValue() > 0 ? (int) nSlider.getValue() : 1;
                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40 / n);
                eExcludeOverlap.numAtoms = n;
                peExcludeOverlap.numAtoms = n;
                keExcludeOverlap.numAtoms = n;

                getDisplayBox(sim.box).repaint();
            }
        };
        nSlider.getSlider().addChangeListener(nListener);
        nListener.stateChanged(null);
        JPanel nSliderPanel = new JPanel(new GridLayout(0, 1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSlider.setShowBorder(false);
        nSlider.setNMajor(4);
        nSliderPanel.add(nSlider.graphic());
        gbc2.gridx = 0;
        gbc2.gridy = 1;
        statePanel.add(nSliderPanel, gbc2);

        //************* Lay out components ****************//

        getDisplayBox(sim.box).setScale(0.7);


        final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
                resetDataAction.actionPerformed();
                mbDistribution.setTemperature(tUnit.toSim(tempSlider.getTemperature()));
                mbSource.update();
                vPlot.doUpdate();
            }
        };
        tempSlider.setSliderPostAction(temperatureAction);
        tempSlider.setRadioGroupPostAction(temperatureAction);

        // show config button
        DeviceButton configButton = new DeviceButton(sim.getController());
        configButton.setLabel("Show Config");
        configButton.setAction(new ActionConfigWindow(sim.box));
        DeviceButton velocityButton = new DeviceButton(sim.getController());
        velocityButton.setLabel("Show Velocities");
        velocityButton.setAction(new ActionVelocityWindow(sim.box));

        final IAction resetAction = new IAction() {
            public void actionPerformed() {
                sim.integrator.reset();

                rdfMeter.reset();

                // Reset density (Density is set and won't change, but
                // do this anyway)
                densityPump.actionPerformed();

                // Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
                tBox.putData(temperatureAverage.getData());

                // IS THIS WORKING?
                pPump.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                peDisplay.putData(peAccumulator.getData());

                getDisplayBox(sim.box).graphic().repaint();

                displayCycles.putData(meterCycles.getData());
            }
        };

        this.getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
                sim.integrator.doThermostat();
                sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SINGLE);
                resetAction.actionPerformed();
            }
        });
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController());

        getPanel().controlPanel.add(setupPanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
        if (showConfig == true) {
            add(configButton);
            add(velocityButton);
        }

        add(rdfPlot);
        add(vPlot);
        add(ePlot);
        add(displayCycles);
        add(densityBox);
        add(tBox);
        add(pDisplay);
        add(peDisplay);
    }

    public void setPotential(String potentialDesc) {
        final boolean HS = potentialDesc.equals("Repulsion only"); 
        final boolean SW = potentialDesc.equals("Repulsion and attraction"); 
        sim.getController().submitActionInterrupt(new IAction() {
            public void actionPerformed() {
                if (HS) {
                    potentialHS.setBox(sim.box);
                    sim.potentialWrapper.setWrappedPotential(potentialHS);
                }
                else if (SW) {
                    potentialSW.setBox(sim.box);
                    sim.potentialWrapper.setWrappedPotential(potentialSW);
                }
                else {
                    potentialIdeal.setBox(sim.box);
                    sim.potentialWrapper.setWrappedPotential(potentialIdeal);
                }
                try {
                    sim.integrator.reset();
                } catch(ConfigurationOverlapException e) {}
            }
        });
    }

    protected class ModifierAtomDiameter implements Modifier {

        public void setValue(double d) {
            if (d > 4.0) {
                throw new IllegalArgumentException("diameter can't exceed 4.0A");
            }
            //assume one type of atom
            ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), d);
            potentialHS.setCollisionDiameter(d);
            potentialSW.setCoreDiameter(d);
            new BoxImposePbc(sim.box, space).actionPerformed();
            try {
                sim.integrator.reset();
            }
            catch (ConfigurationOverlapException e){
                // can happen when increasing diameter
            }
            sigma = d;
            getDisplayBox(sim.box).repaint();
        }

        public double getValue() {
            return sigma;
        }

        public Dimension getDimension() {
            return Length.DIMENSION;
        }
        
        public String getLabel() {
            return "Atom Diameter";
        }
        
        public String toString() {
            return getLabel();
        }
    }
    
    public static class DataSinkExcludeOverlap extends DataProcessor {

        public DataSinkExcludeOverlap() {
            myData = new DataDouble();
        }
        
        public IData processData(IData data) {
            if (Double.isInfinite(data.getValue(0))) {
                return null;
            }
            myData.E(data);
            myData.TE(1.0/numAtoms);
            return myData;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            return inputDataInfo;
        }
        
        public int numAtoms;
        protected final DataDouble myData;
    }

    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }

        SwmdGraphic swmdGraphic = new SwmdGraphic(new Swmd(space), space);
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
            SwmdGraphic swmdGraphic = new SwmdGraphic(new Swmd(sp), sp);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
