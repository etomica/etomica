/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.*;
import etomica.data.types.DataDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2Electrostatic;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.units.systems.MKS;
import etomica.util.Constants.CompassDirection;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.util.ArrayList;

/**
 * Graphical and data collection bits for reverse osmosis simulation.
 * 
 * @author Andrew Schultz
 */
public class ReverseOsmosisWaterGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Osmosis";
    private final static int REPAINT_INTERVAL = 2;
    protected DeviceThermoSlider tempSlider;
    protected DeviceSlider solventChamberDensitySlider, soluteChamberDensitySlider;
    protected DeviceSlider soluteMoleFractionSlider;
    protected DeviceBox sigBox, epsBox, massBox, tetherBox;
    protected DeviceSlider membraneThicknessSlider, membraneWidthSlider;
    protected DeviceSlider soluteChargeSlider;
    protected Unit eUnit, dUnit, pUnit;
    protected ReverseOsmosisWater sim;
    
    public ReverseOsmosisWaterGraphic(final ReverseOsmosisWater simulation) {

        super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = getController().getSimRestart().getDataResetAction();

        this.sim = simulation;
        sim.integrator.printInterval = 1000;

        Unit tUnit = Kelvin.UNIT;

        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);

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

        ((SimulationRestart) getController().getReinitButton().getAction()).setConfiguration(sim.configMembrane);

        sim.getController().setSleepPeriod(0);
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

        //combo box to select potentials
        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        massBox = new DeviceBox();
        tetherBox = new DeviceBox();
        membraneThicknessSlider = new DeviceSlider(sim.getController());
        membraneWidthSlider = new DeviceSlider(sim.getController());
        soluteChargeSlider = new DeviceSlider(sim.getController());

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump = new DataPump(meterCycles, displayCycles);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
        displayCycles.setLabel("Simulation time");

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setUnit(Kelvin.UNIT);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(500.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setUnit(tUnit);
        tempSlider.setAdiabatic();
        tempSlider.setSliderPostAction(resetDataAction);
        tempSlider.setRadioGroupPostAction(resetDataAction);

        ModifierGeneral modifier = new ModifierGeneral(sim.configMembrane, "solventChamberDensity");
        solventChamberDensitySlider = new DeviceSlider(sim.getController(), modifier);
        solventChamberDensitySlider.setPrecision(1);
        solventChamberDensitySlider.setMaximum(80);
        solventChamberDensitySlider.setNMajor(3);
        solventChamberDensitySlider.setShowValues(true);
        solventChamberDensitySlider.setUnit(dUnit);
        solventChamberDensitySlider.doUpdate();
        solventChamberDensitySlider.setShowBorder(true);
        solventChamberDensitySlider.setLabel("Solvent Density (mol/L)");

        modifier = new ModifierGeneral(sim.configMembrane, "solutionChamberDensity");
        soluteChamberDensitySlider = new DeviceSlider(sim.getController(), modifier);
        soluteChamberDensitySlider.setPrecision(1);
        soluteChamberDensitySlider.setMaximum(80);
        soluteChamberDensitySlider.setNMajor(3);
        soluteChamberDensitySlider.setShowValues(true);
        soluteChamberDensitySlider.setUnit(dUnit);
        soluteChamberDensitySlider.doUpdate();
        soluteChamberDensitySlider.setShowBorder(true);
        soluteChamberDensitySlider.setLabel("Solution Density (mol/L)");

        modifier = new ModifierGeneral(sim.configMembrane, "soluteMoleFraction");
        soluteMoleFractionSlider = new DeviceSlider(sim.getController(), modifier);
        soluteMoleFractionSlider.setPrecision(2);
        soluteMoleFractionSlider.setMaximum(1);
        soluteMoleFractionSlider.setNMajor(3);
        soluteMoleFractionSlider.setValue(sim.configMembrane.getSoluteMoleFraction());
        soluteMoleFractionSlider.setShowValues(true);
        soluteMoleFractionSlider.setShowBorder(true);
        soluteMoleFractionSlider.setLabel("Solute Mole Fraction");

        JPanel configPanel = new JPanel(new GridBagLayout());
        configPanel.add(solventChamberDensitySlider.graphic(), vertGBC);
        configPanel.add(soluteChamberDensitySlider.graphic(), vertGBC);
        configPanel.add(soluteMoleFractionSlider.graphic(), vertGBC);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.setBorder(new TitledBorder(null, "Potential Selection", TitledBorder.CENTER, TitledBorder.TOP));
        JPanel parameterPanel = new JPanel(new GridLayout(0, 1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(massBox.graphic());
        parameterPanel.add(tetherBox.graphic());
        potentialPanel.add(parameterPanel, vertGBC);
        potentialPanel.add(membraneThicknessSlider.graphic(), vertGBC);
        potentialPanel.add(membraneWidthSlider.graphic(), vertGBC);

        // inline class to make water's OO potential look like a P2LennardJones
        P2LennardJones p2LJOO = new P2LennardJones(sim.getSpace()) {
            public double getSigma() {
                return sim.potentialWater.getSigma();
            }

            public double getEpsilon() {
                return sim.potentialWater.getEpsilon();
            }
        };
        DiameterHashByType diameterManager = (DiameterHashByType) getDisplayBox(sim.box).getDiameterHash();
        diameterManager.setDiameter(sim.speciesMembrane.getLeafType(), sim.potentialMM.getSigma());
        ModifierAtomDiameter sigMembraneModifier = new ModifierAtomDiameter(sim.speciesMembrane, sim.potentialMM,
                new P2LennardJones[]{sim.potentialMO, sim.potentialMNa, sim.potentialMCl},
                new P2LennardJones[]{p2LJOO, sim.potentialLJNaNa, sim.potentialLJClCl}, diameterManager);
        ModifierEpsilon epsModifier = new ModifierEpsilon(sim.potentialMM,
                new P2LennardJones[]{sim.potentialMO, sim.potentialMNa, sim.potentialMCl},
                new P2LennardJones[]{p2LJOO, sim.potentialLJNaNa, sim.potentialLJClCl});
        ModifierGeneral massModifier = new ModifierGeneral(sim.speciesMembrane.getLeafType().getElement(), "mass");
        ModifierGeneral tetherModifier = new ModifierGeneral(sim.potentialTether, "epsilon");
        Modifier membraneThicknessModifier = new Modifier() {
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "membrane thickness";
            }

            public double getValue() {
                return sim.configMembrane.getNumMembraneLayers();
            }

            public void setValue(double newValue) {
                sim.configMembrane.setNumMembraneLayers((int) newValue);
            }
        };
        Modifier membraneWidthModifier = new Modifier() {
            public Dimension getDimension() {
                return Quantity.DIMENSION;
            }

            public String getLabel() {
                return "membrane width";
            }

            public double getValue() {
                return sim.configMembrane.getMembraneWidth();
            }

            public void setValue(double newValue) {
                sim.configMembrane.setMembraneWidth((int) newValue);
            }
        };
        sigBox.setModifier(sigMembraneModifier);
        sigBox.setLabel("Core Diameter (" + Angstrom.UNIT.symbol() + ")");
        sigBox.doUpdate();
        epsBox.setUnit(eUnit);
        epsBox.setModifier(epsModifier);
        epsBox.setPostAction(resetDataAction);
        massBox.setModifier(massModifier);
        massBox.setUnit(Dalton.UNIT);
        massBox.setPostAction(resetDataAction);
        tetherBox.setUnit(eUnit);
        tetherBox.setModifier(tetherModifier);
        tetherBox.setLabel("Tether Constant (" + eUnit.symbol() + ")");
        tetherBox.setPostAction(resetDataAction);
        sigBox.setController(sim.getController());
        epsBox.setController(sim.getController());
        massBox.setController(sim.getController());
        tetherBox.setController(sim.getController());
        membraneThicknessSlider.setLabel("Membrane Thickness");
        membraneThicknessSlider.setShowBorder(true);
        membraneThicknessSlider.setMinimum(1);
        membraneThicknessSlider.setMaximum(4);
        membraneThicknessSlider.setModifier(membraneThicknessModifier);
        membraneWidthSlider.setLabel("Membrane Width");
        membraneWidthSlider.setShowBorder(true);
        membraneWidthSlider.setMinimum(2);
        membraneWidthSlider.setMaximum(4);
        membraneWidthSlider.setModifier(membraneWidthModifier);

        ModifierSoluteCharge soluteChargeModifier = new ModifierSoluteCharge(sim.potentialQNaNa, sim.potentialQClCl,
                sim.potentialQNaCl, new P2Electrostatic[]{sim.potentialQHNa, sim.potentialQONa},
                new P2Electrostatic[]{sim.potentialQHCl, sim.potentialQOCl});
        soluteChargeSlider.setLabel("Solute Charge");
        soluteChargeSlider.setShowBorder(true);
        soluteChargeSlider.setMinimum(0);
        soluteChargeSlider.setMaximum(2);
        soluteChargeSlider.setUnit(Electron.UNIT);
        soluteChargeSlider.setModifier(soluteChargeModifier);
        soluteChargeSlider.setPostAction(resetDataAction);

        JPanel solutePanel = new JPanel(new GridBagLayout());
        solutePanel.add(soluteChargeSlider.graphic(), vertGBC);

        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Membrane");
        setupPanel.add(configPanel, "Configuration");
        setupPanel.add(solutePanel, "Solute");

        IAction neighborRangeReset = new IAction() {
            public void actionPerformed() {
//                ((PotentialMasterList)sim.integrator.getPotential()).reset();
//                double nbrRange = ((PotentialMasterList)sim.integrator.getPotential()).getMaxPotentialRange();
//                nbrRange *= 1.2;
//                ((PotentialMasterList)sim.integrator.getPotential()).setRange(nbrRange);
//                ((PotentialMasterList)sim.integrator.getPotential()).reset();

                resetDataAction.actionPerformed();
                try {
                    sim.integrator.reset();
                } catch (ConfigurationOverlapException e) {
                }
                getDisplayBox(sim.box).repaint();
            }
        };
        sigBox.setPostAction(neighborRangeReset);

        //display of box, timer
        ColorSchemeByType colorScheme = (ColorSchemeByType) getDisplayBox(sim.box).getColorScheme();
        colorScheme.setColor(sim.speciesSodium.getLeafType(), Color.BLUE);
        colorScheme.setColor(sim.speciesChlorine.getLeafType(), Color.GREEN);
        colorScheme.setColor(sim.speciesSolvent.getOxygenType(), Color.RED);
        colorScheme.setColor(sim.speciesSolvent.getHydrogenType(), Color.WHITE);
        colorScheme.setColor(sim.speciesMembrane.getLeafType(), Color.CYAN);

        //meters and displays
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

        MeterTemperature thermometer = new MeterTemperature(sim, sim.box, space.D());
        thermometer.setKineticEnergyMeter(new MeterKineticEnergyFromIntegrator(sim.integrator));
        final DisplayTextBox tBox = new DisplayTextBox();
        final DataPump temperaturePump = new DataPump(thermometer, tBox);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(10);
        dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
        tBox.setLabel("Measured Temperature");
        tBox.setLabelPosition(CompassDirection.NORTH);

        MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
        eMeter.setKinetic(new MeterKineticEnergyFromIntegrator(sim.integrator));
        final AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataSinkExcludeOverlap foo = new DataSinkExcludeOverlap(sim.box);
        final DataPump energyPump = new DataPump(eMeter, foo);
        foo.setDataSink(energyHistory);
        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        sim.integrator.getEventManager().addListener(energyPumpListener);
        energyPumpListener.setInterval(10);
        dataStreamPumps.add(energyPump);

        MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster(), sim.box);
        final AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        foo = new DataSinkExcludeOverlap(sim.box);
        final DataPump pePump = new DataPump(peMeter, foo);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        foo.setDataSink(peFork);
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        sim.integrator.getEventManager().addListener(pePumpListener);
        pePumpListener.setInterval(10);
        dataStreamPumps.add(pePump);

        MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        foo = new DataSinkExcludeOverlap(sim.box);
        final DataPump kePump = new DataPump(keMeter, foo);
        foo.setDataSink(keHistory);
        IntegratorListenerAction kePumpListener = new IntegratorListenerAction(kePump);
        sim.integrator.getEventManager().addListener(kePumpListener);
        kePumpListener.setInterval(10);
        dataStreamPumps.add(kePump);

        MeterFlux meterFlux = new MeterFlux(sim, space);
        double xLength = sim.box.getBoundary().getBoxSize().getX(0);
        meterFlux.setBoundaries(0, new double[]{-0.25 * xLength, 0.25 * xLength}, new int[]{1, -1});
        meterFlux.setIntegrator(sim.integrator);
        meterFlux.setBox(sim.box);
        meterFlux.setSpecies(new ISpecies[]{sim.speciesSodium, sim.speciesChlorine, sim.speciesSolvent});
        AccumulatorHistory fluxHistory = new AccumulatorHistory(new HistoryCollapsingAverage(20));
        fluxHistory.setTimeDataSource(timeCounter);
        DataPump fluxPump = new DataPump(meterFlux, fluxHistory);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(fluxPump));
        dataStreamPumps.add(fluxPump);

        MeterOsmoticPressure meterOsmoticPressure = new MeterOsmoticPressure(sim.torqueSum, sim.box);
        AccumulatorAverageCollapsing pressureAvg = new AccumulatorAverageCollapsing();
        pressureAvg.setPushInterval(10);
        DataFork pressureFork = new DataFork();
        DataPump pressurePump = new DataPump(meterOsmoticPressure, pressureFork);
        pressureFork.addDataSink(pressureAvg);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pressurePump));
        dataStreamPumps.add(pressurePump);
        AccumulatorHistory pressureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        pressureFork.addDataSink(pressureHistory);
        pressureHistory.setTimeDataSource(timeCounter);

        MeterNMolecules densitySodium = new MeterNMolecules();
        MeterNMolecules densityChlorine = new MeterNMolecules();
        MeterNMolecules densitySolvent = new MeterNMolecules();
        densitySodium.setBox(sim.box);
        densityChlorine.setBox(sim.box);
        densitySolvent.setBox(sim.box);
        densitySodium.setSpecies(sim.speciesSodium);
        densityChlorine.setSpecies(sim.speciesChlorine);
        densitySolvent.setSpecies(sim.speciesSolvent);
        MeterProfileByVolume meterProfileSodium = new MeterProfileByVolume(sim.getSpace());
        MeterProfileByVolume meterProfileChlorine = new MeterProfileByVolume(sim.getSpace());
        MeterProfileByVolume meterProfileSolvent = new MeterProfileByVolume(sim.getSpace());
        meterProfileSodium.setBox(sim.box);
        meterProfileSodium.setDataSource(densitySodium);
        meterProfileChlorine.setBox(sim.box);
        meterProfileChlorine.setDataSource(densityChlorine);
        meterProfileSolvent.setBox(sim.box);
        meterProfileSolvent.setDataSource(densitySolvent);
        meterProfileSodium.setProfileDim(0);
        meterProfileChlorine.setProfileDim(0);
        meterProfileSolvent.setProfileDim(0);
        AccumulatorAverageFixed profSodiumAvg = new AccumulatorAverageFixed(10);
        profSodiumAvg.setPushInterval(10);
        DataPump profPumpSodium = new DataPump(meterProfileSodium, profSodiumAvg);
        IntegratorListenerAction profPumpSodiumListener = new IntegratorListenerAction(profPumpSodium);
        sim.integrator.getEventManager().addListener(profPumpSodiumListener);
        profPumpSodiumListener.setInterval(2);
        AccumulatorAverageFixed profChlorineAvg = new AccumulatorAverageFixed(10);
        profChlorineAvg.setPushInterval(10);
        DataPump profPumpChlorine = new DataPump(meterProfileChlorine, profChlorineAvg);
        IntegratorListenerAction profPumpChlorineListener = new IntegratorListenerAction(profPumpChlorine);
        sim.integrator.getEventManager().addListener(profPumpChlorineListener);
        profPumpChlorineListener.setInterval(2);
        AccumulatorAverageFixed profSolventAvg = new AccumulatorAverageFixed(10);
        profSolventAvg.setPushInterval(10);
        DataPump profPumpSolvent = new DataPump(meterProfileSolvent, profSolventAvg);
        IntegratorListenerAction profPumpSolventListener = new IntegratorListenerAction(profPumpSolvent);
        sim.integrator.getEventManager().addListener(profPumpSolventListener);
        profPumpSolventListener.setInterval(2);
        dataStreamPumps.add(profPumpSodium);
        dataStreamPumps.add(profPumpChlorine);
        dataStreamPumps.add(profPumpSolvent);

        final DisplayPlotXChart ePlot = new DisplayPlotXChart();
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
//        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");

        ePlot.getPlot().setTitle("Energy History (J/mol)");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy");
        ePlot.setXUnit(Picosecond.UNIT);

        final DisplayPlotXChart profPlot = new DisplayPlotXChart();
        profSodiumAvg.addDataSink(profPlot.getDataSet().makeDataSink(),
                new AccumulatorAverage.StatType[]{profSolventAvg.AVERAGE});
        profPlot.setLegend(new DataTag[]{meterProfileSodium.getTag()}, "Sodium");
        profChlorineAvg.addDataSink(profPlot.getDataSet().makeDataSink(),
                new AccumulatorAverage.StatType[]{profSolventAvg.AVERAGE});
        profPlot.setLegend(new DataTag[]{meterProfileChlorine.getTag()}, "Chlorine");
        profSolventAvg.addDataSink(profPlot.getDataSet().makeDataSink(),
                new AccumulatorAverage.StatType[]{profSolventAvg.AVERAGE});
        profPlot.setLegend(new DataTag[]{meterProfileSolvent.getTag()}, "Solvent");
        profPlot.getPlot().setTitle("DensityProfiles (mol/L)");
        profPlot.setDoLegend(true);
        profPlot.setLabel("Profiles");
        profPlot.setXUnit(Angstrom.UNIT);
        profPlot.setUnit(dUnit);

        final DisplayPlotXChart pPlot = new DisplayPlotXChart();
        pressureHistory.setDataSink(pPlot.getDataSet().makeDataSink());
        pPlot.setLabel("Osmotic Pressure");
        pPlot.setUnit(Bar.UNIT);

        final DisplayPlotXChart fluxPlot = new DisplayPlotXChart();
        fluxHistory.setDataSink(fluxPlot.getDataSet().makeDataSink());
        fluxPlot.setLabel("Flux");

//        MeterPressureHard pMeter = new MeterPressureHard(sim.getSpace());
//        pMeter.setIntegrator(sim.integrator);
//        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
//        final DataPump pPump = new DataPump(pMeter, pAccumulator);
//        sim.integrator.addIntervalAction(pPump);
//        pAccumulator.setPushInterval(50);
//        dataStreamPumps.add(pPump);
//
//        final DisplayTextBoxesCAE pDisplay = new DisplayTextBoxesCAE();
//        pDisplay.setLabel(sim.getSpace().D() == 3 ? "Pressure (bar)" : "Pressure (bar-nm)");
//        pDisplay.setAccumulator(pAccumulator);
//        pDisplay.setUnit(pUnit);
        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);
        peDisplay.setLabel("Potential Energy (J/mol)");

        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
        pressureDisplay.setAccumulator(pressureAvg);
        pressureDisplay.setUnit(Bar.UNIT);
        pressureDisplay.setLabel("Osmotic Pressure (bar)");

//        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
//        nSlider.setResetAction(new SimulationRestart(sim));
//        nSlider.setSpecies(sim.species);
//        nSlider.setBox(sim.box);
//        nSlider.setMinimum(0);
//        nSlider.setMaximum(sim.getSpace().D() == 3 ? 500 : 168);
//        nSlider.setLabel("Number of Atoms");
//        nSlider.setShowBorder(true);
//        nSlider.setShowValues(true);
//        // add a listener to adjust the thermostat interval for different
//        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
//        // don't need as much thermostating.
//        ChangeListener nListener = new ChangeListener() {
//            public void stateChanged(ChangeEvent evt) {
//                final int n = (int)nSlider.getValue() > 0 ? (int)nSlider.getValue() : 1;
////                sim.activityIntegrate.setSleepPeriod(n > 400 ? 0 : 1);
//                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40/n);
////                sim.integrator.setActionInterval(energyPump, 1);
////                sim.integrator.setActionInterval(kePump, 1);
////                sim.integrator.setActionInterval(pePump, 1);
//                eExcludeOverlap.numAtoms = n;
//                peExcludeOverlap.numAtoms = n;
//                keExcludeOverlap.numAtoms = n;
//
//                getDisplayBox(sim.box).repaint();
//            }
//        };
//        nSlider.getSlider().addChangeListener(nListener);
//        nListener.stateChanged(null);
//        JPanel nSliderPanel = new JPanel(new GridLayout(0,1));
//        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
//        nSlider.setShowBorder(false);
//        nSlider.setNMajor(4);
//        nSliderPanel.add(nSlider.graphic());
//        gbc2.gridx = 0;  gbc2.gridy = 1;
//        statePanel.add(nSliderPanel, gbc2);
//
        //************* Lay out components ****************//

        getDisplayBox(sim.box).setScale(0.7);


        final IAction resetAction = new IAction() {
            public void actionPerformed() {
                try {
                    sim.integrator.reset();
                } catch (ConfigurationOverlapException e) {
                }

                // Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();

                // IS THIS WORKING?
//                pPump.actionPerformed();
//                pDisplay.putData(pAccumulator.getData());
//                pDisplay.repaint();
                peDisplay.putData(peAccumulator.getData());

                getDisplayBox(sim.box).graphic().repaint();

                displayCycles.putData(meterCycles.getData());
            }
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);
        IAction reconfigAction = new IAction() {
            public void actionPerformed() {
                sim.configMembrane.initializeCoordinates(sim.box);
                resetDataAction.actionPerformed();
                resetAction.actionPerformed();
            }
        };
        solventChamberDensitySlider.setPostAction(reconfigAction);
        soluteChamberDensitySlider.setPostAction(reconfigAction);
        soluteMoleFractionSlider.setPostAction(reconfigAction);
        membraneThicknessSlider.setPostAction(reconfigAction);
        membraneWidthSlider.setPostAction(reconfigAction);

        soluteChamberDensitySlider.setValue(16);
        solventChamberDensitySlider.setValue(16);

        getPanel().controlPanel.add(setupPanel, vertGBC);

        add(ePlot);
        add(profPlot);
        add(pPlot);
        add(fluxPlot);
        add(displayCycles);
        add(tBox);
        add(peDisplay);
        add(pressureDisplay);
    }

    protected static class ModifierAtomDiameter implements Modifier {
        
        public ModifierAtomDiameter(SpeciesSpheresMono species, P2LennardJones potential,
                P2LennardJones[] crossPotentials, P2LennardJones[] otherPurePotentials, DiameterHashByType diameterManager) {
            this.species = species;
            this.potential = potential;
            this.crossPotentials = crossPotentials;
            this.otherPurePotentials = otherPurePotentials;
            this.diameterHash = diameterManager;
        }

        public void setValue(double d) {
            if (d > 4.5) {
                throw new IllegalArgumentException("diameter can't exceed 4.0A");
            }

            diameterHash.setDiameter(species.getLeafType(), d);
            potential.setSigma(d);
            for (int i=0; i<crossPotentials.length; i++) {
                double otherSigma = otherPurePotentials[i].getSigma();
                crossPotentials[i].setSigma(0.5*(d+otherSigma));
            }
        }

        public double getValue() {
            return potential.getSigma();
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
        
        protected final SpeciesSpheresMono species;
        protected final P2LennardJones potential;
        protected final P2LennardJones[] crossPotentials;
        protected final P2LennardJones[] otherPurePotentials;
        protected final DiameterHashByType diameterHash;
    }
    
    protected static class ModifierEpsilon implements Modifier {
        
        public ModifierEpsilon(P2LennardJones potential,
                P2LennardJones[] otherPurePotentials, P2LennardJones[] mixPotentials) {
            this.potential = potential;
            this.otherPurePotentials = otherPurePotentials;
            this.mixPotentials = mixPotentials;
        }

        public void setValue(double epsilon) {
            potential.setEpsilon(epsilon);
            for (int i=0; i<mixPotentials.length; i++) {
                double otherEpsilon = otherPurePotentials[i].getEpsilon();
                mixPotentials[i].setEpsilon(Math.sqrt(epsilon*otherEpsilon));
            }
        }

        public double getValue() {
            return potential.getEpsilon();

        }

        public Dimension getDimension() {
            return Energy.DIMENSION;
        }
        
        public String getLabel() {
            return "Epsilon";
        }
        
        protected final P2LennardJones potential;
        protected final P2LennardJones[] otherPurePotentials;
        protected final P2LennardJones[] mixPotentials;
    }
    
    protected static class ModifierSoluteCharge implements Modifier {
        
        public ModifierSoluteCharge(P2Electrostatic pSolutePP, P2Electrostatic pSoluteMM,
                P2Electrostatic pSolutePM,
                P2Electrostatic[] mixPotentialsP, P2Electrostatic[] mixPotentialsM) {
            this.pSolutePP = pSolutePP;
            this.pSoluteMM = pSolutePM;
            this.pSolutePM = pSolutePM;
            this.mixPotentialsP = mixPotentialsP;
            this.mixPotentialsM = mixPotentialsM;
        }

        public void setValue(double newCharge) {
            pSolutePP.setCharge1(newCharge);
            pSolutePP.setCharge2(newCharge);
            pSolutePM.setCharge1(newCharge);
            pSolutePM.setCharge2(-newCharge);
            pSoluteMM.setCharge1(-newCharge);
            pSoluteMM.setCharge2(-newCharge);
            for (int i=0; i<mixPotentialsP.length; i++) {
                mixPotentialsP[i].setCharge2(newCharge);
            }
            for (int i=0; i<mixPotentialsM.length; i++) {
                mixPotentialsM[i].setCharge2(-newCharge);
            }
        }

        public double getValue() {
            return pSolutePP.getCharge1();

        }

        public Dimension getDimension() {
            return Charge.DIMENSION;
        }
        
        public String getLabel() {
            return "Charge";
        }
        
        protected final P2Electrostatic pSolutePP, pSoluteMM, pSolutePM;
        protected final P2Electrostatic[] mixPotentialsP, mixPotentialsM;
    }
    
    public static class DataSinkExcludeOverlap extends DataProcessor {

        public DataSinkExcludeOverlap(Box box) {
            myData = new DataDouble();
            this.box = box;
        }
        
        public IData processData(IData data) {
            myData.E(data);
            myData.TE(1.0/box.getLeafList().size());
            return myData;
        }

        protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
            return inputDataInfo;
        }
        
        protected final DataDouble myData;
        protected final Box box;
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        ReverseOsmosisWaterGraphic reverseOsmosisGraphic = new ReverseOsmosisWaterGraphic(new ReverseOsmosisWater());
		SimulationGraphic.makeAndDisplayFrame
		        (reverseOsmosisGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            ReverseOsmosisWaterGraphic reverseOsmosisGraphic = new ReverseOsmosisWaterGraphic(new ReverseOsmosisWater());

		    getContentPane().add(reverseOsmosisGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}


