/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ljmd;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.*;
import etomica.graphics.*;
import etomica.integrator.ActionZeroMomentum;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.math.DoubleRange;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.statmech.MaxwellBoltzmann;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Time;
import etomica.util.Constants.CompassDirection;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.awt.*;
import java.util.ArrayList;

public class LjmdGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Lennard-Jones Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 20;
    private final DeviceThermoSlider temperatureSelect;
    protected Ljmd sim;

    private final boolean showConfig = false;

    public LjmdGraphic(final Ljmd simulation) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        this.sim = simulation;

        sim.getController().setSleepPeriod(1);
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getLeafType(), Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeByType());
//        sim.integrator.addListener(new IntervalActionAdapter(this.getDisplayBoxPaintAction(sim.box)));

        //meters and displays
        final MeterRDF rdfMeter = new MeterRDF(sim.getSpace());
        IntegratorListenerAction rdfMeterListener = new IntegratorListenerAction(rdfMeter);
        sim.integrator.getEventManager().addListener(rdfMeterListener);
        rdfMeterListener.setInterval(10);
        rdfMeter.getXDataSource().setXMax(4.0);
        rdfMeter.setBox(sim.box);
        DisplayPlotXChart rdfPlot = new DisplayPlotXChart();
        DataPump rdfPump = new DataPump(rdfMeter, rdfPlot.getDataSet().makeDataSink());
        IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
        sim.integrator.getEventManager().addListener(rdfPumpListener);
        rdfPumpListener.setInterval(10);
        dataStreamPumps.add(rdfPump);

        rdfPlot.setDoLegend(false);
        rdfPlot.getPlot().setTitle("Radial Distribution Function");
        rdfPlot.setLabel("RDF");

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                rdfMeter.reset();
            }
        };

        //velocity distribution
        double vMin = 0;
        double vMax = 4;
        DataSourceRmsVelocity meterVelocity = new DataSourceRmsVelocity(sim.box, new HistogramSimple(100, new DoubleRange(0, 4)));
        AccumulatorAverage rmsAverage = new AccumulatorAverageFixed(10);
        DataPump velocityPump = new DataPump(meterVelocity, rmsAverage);
        IntegratorListenerAction velocityPumpListener = new IntegratorListenerAction(velocityPump);
        sim.integrator.getEventManager().addListener(velocityPumpListener);
        velocityPumpListener.setInterval(10);
        rmsAverage.setPushInterval(10);
        dataStreamPumps.add(velocityPump);

        final DisplayPlotXChart vPlot = new DisplayPlotXChart();
        rmsAverage.addDataSink(vPlot.getDataSet().makeDataSink(), new StatType[]{rmsAverage.AVERAGE});
        vPlot.setDoLegend(false);
        vPlot.getPlot().setTitle("Velocity Distribution");
        vPlot.setDoLegend(true);
        vPlot.setLabel("Velocity");

        final MaxwellBoltzmann.Distribution mbDistribution = new MaxwellBoltzmann.Distribution(sim.getSpace(), sim.integrator.getTemperature(), sim.species.getLeafType().getMass());
        final DataSourceFunction mbSource = new DataSourceFunction("Maxwell Boltzmann Distribution",
                Null.DIMENSION, mbDistribution, 100, "Speed", new DimensionRatio(Length.DIMENSION, Time.DIMENSION));
        DataSourceUniform mbX = mbSource.getXSource();
        mbX.setTypeMax(LimitType.HALF_STEP);
        mbX.setTypeMin(LimitType.HALF_STEP);
        mbX.setNValues(meterVelocity.getDataInfo().getLength());
        mbX.setXMin(vMin);
        mbX.setXMax(vMax);
        mbSource.update();
        final DataPump mbPump = new DataPump(mbSource, vPlot.makeSink("boltzmann"));
        vPlot.getSeries("boltzmann")
                .setMarker(SeriesMarkers.NONE);

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

        MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer, temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(10);
        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
        temperatureHistory.setTimeDataSource(timeCounter);
        final DisplayTextBox tBox = new DisplayTextBox();
        temperatureFork.setDataSinks(new IDataSink[]{tBox, temperatureHistory});
        tBox.setLabel("Measured Temperature");
        tBox.setLabelPosition(CompassDirection.NORTH);

        dataStreamPumps.add(temperaturePump);
        tBox.setLabel("Measured Temperature");
        tBox.setLabelPosition(CompassDirection.NORTH);

        // Number density box
        MeterDensity densityMeter = new MeterDensity(sim.box);
        final DisplayTextBox densityBox = new DisplayTextBox();
        final DataPump densityPump = new DataPump(densityMeter, densityBox);
        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
        sim.integrator.getEventManager().addListener(densityPumpListener);
        densityPumpListener.setInterval(10);
        dataStreamPumps.add(densityPump);
        densityBox.setLabel("Number Density");

        MeterEnergyFromIntegrator eMeter = new MeterEnergyFromIntegrator(sim.integrator);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataPump energyPump = new DataPump(eMeter, energyHistory);
        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        sim.integrator.getEventManager().addListener(energyPumpListener);
        energyPumpListener.setInterval(60);
        energyHistory.setPushInterval(5);
        dataStreamPumps.add(energyPump);

        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
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

        MeterKineticEnergy keMeter = new MeterKineticEnergy(sim.box);
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

        DisplayPlotXChart ePlot = new DisplayPlotXChart();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy");

        MeterPressureFromIntegrator pMeter = new MeterPressureFromIntegrator(sim.integrator);
        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
        final DataPump pPump = new DataPump(pMeter, pAccumulator);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pPump));
        pAccumulator.setPushInterval(10);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
        pDisplay.setAccumulator(pAccumulator);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);

        IAction mbAction = new IAction() {
            public void actionPerformed() {
                if (!sim.integrator.isIsothermal()) {
                    // in adiabatic mode, we want the average kinetic temperature, which we'll
                    // steal from our kinetic energy accumulator
                    double temp = keAvg.getData().getValue(keAvg.AVERAGE.index);
                    temp *= 2.0 / (sim.box.getLeafList().size() * space.D());
                    mbDistribution.setTemperature(temp);
                    mbSource.update();
                }
                mbPump.actionPerformed();
            }
        };

        IntegratorListenerAction mbActionListener = new IntegratorListenerAction(mbAction);
        sim.integrator.getEventManager().addListener(mbActionListener);
        mbActionListener.setInterval(100);

        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(space.D() == 2 ? 225 : 4000);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        nSlider.setPostAction(new IAction() {
            final ConfigurationLattice config = new ConfigurationLattice((space.D() == 2) ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space);
            int oldN = sim.box.getMoleculeList().size();

            public void actionPerformed() {
                int n = (int) nSlider.getValue();
                if (n == 0) {
                    sim.integrator.setThermostatInterval(400);
                } else {
                    sim.integrator.setThermostatInterval(400 / n + 1);
                }

                if (oldN < n) {
                    config.initializeCoordinates(sim.box);
                }
                oldN = n;
                sim.integrator.reset();
                if (!sim.integrator.isIsothermal() && n > 1) {
                    ActionZeroMomentum.zeroMomenta(sim.box);
                }
                resetDataAction.actionPerformed();
                getDisplayBox(sim.box).repaint();
            }

        });

        //************* Lay out components ****************//

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getDisplayBox(sim.box).setScale(0.7);

        //temperature selector
        temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
        temperatureSelect.setPrecision(1);
        temperatureSelect.setMinimum(0.0);
        temperatureSelect.setMaximum(10.0);
        temperatureSelect.setSliderMajorValues(3);
        temperatureSelect.setAdiabatic();

        final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
                resetDataAction.actionPerformed();
                mbDistribution.setTemperature(temperatureSelect.getTemperature());
                mbSource.update();
                mbPump.actionPerformed();
                vPlot.doUpdate();
            }
        };

        temperatureSelect.setSliderPostAction(temperatureAction);
        temperatureSelect.setRadioGroupPostAction(temperatureAction);

        // show config button
        DeviceButton configButton = new DeviceButton(sim.getController());
        configButton.setLabel("Show Config");
        configButton.setAction(new ActionConfigWindow(sim.box));

        IAction resetAction = new IAction() {
            public void actionPerformed() {
                rdfMeter.reset();

                // Reset density (Density is set and won't change, but
                // do this anyway)
                densityPump.actionPerformed();

                // Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
//                tBox.putData(temperatureHistory.getData());

                // IS THIS WORKING?
                pPump.actionPerformed();
                pDisplay.putData(pAccumulator.getData());
                peDisplay.putData(peAccumulator.getData());

                getDisplayBox(sim.box).graphic().repaint();
            }
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        getPanel().controlPanel.add(temperatureSelect.graphic(), vertGBC);
        add(nSlider);
        if (showConfig == true) {
            add(configButton);
        }

        add(rdfPlot);
        add(vPlot);
        add(ePlot);
        add(densityBox);
        add(tBox);
        add(pDisplay);
        add(peDisplay);
    }

    public static void main(String[] args) {
        Space sp = null;
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    sp = Space3D.getInstance();
                } else {
                    sp = Space2D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        } else {
            sp = Space3D.getInstance();
        }

        LjmdGraphic ljmdGraphic = new LjmdGraphic(new Ljmd(sp));
        SimulationGraphic.makeAndDisplayFrame
                (ljmdGraphic.getPanel(), APP_NAME);
    }
}


