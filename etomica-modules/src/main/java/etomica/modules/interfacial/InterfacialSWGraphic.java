/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.interfacial;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.action.IAction;
import etomica.api.*;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.MoleculePositionGeometricCenterPBC;
import etomica.atom.DiameterHashByType;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataDump;
import etomica.data.DataFork;
import etomica.data.DataGroupSplitter;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataProcessorChemicalPotential;
import etomica.data.DataProcessorInterfacialTension;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourcePositionedBoltzmannFactor;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSink;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterProfile;
import etomica.data.meter.MeterProfileByAtoms;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTensor;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.Modifier;
import etomica.modules.interfacial.DataSourceTensorVirialHardProfile.DataSourceVirialProfile;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Pixel;
import etomica.units.Unit;
import etomica.units.systems.LJ;
import etomica.util.Constants.CompassDirection;
import etomica.data.history.HistoryCollapsingAverage;

/**
 * Graphic UI for interfacial tension module.  Design by Heath Turner.
 *
 * @author Andrew Schultz
 */
public class InterfacialSWGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Interfacial Tension";
    private final static int REPAINT_INTERVAL = 20;
    private DeviceThermoSlider temperatureSelect;
    protected InterfacialSW sim;
    protected final DeviceNSelector nSlider;
    protected final DeviceSlider nSurfactantSlider;
    protected final DeviceSlider xSlider;
    protected final MeterProfileByVolume densityProfileMeter, surfactantProfileMeter;
    protected final MeterProfileByAtoms orientationProfileMeter;
    protected final MeterProfile muProfileMeter;
    protected boolean isExpanded;
    protected final DataPumpListener surfactantProfilePump, orientationProfilePump;
    protected final DisplayPlot profilePlot, orientationPlot;
    
    public InterfacialSWGraphic(final InterfacialSW simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, _space.D() == 2 ? 10*REPAINT_INTERVAL : REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

        LJ unitSystem = new LJ();
        Unit tUnit = Energy.DIMENSION.getUnit(unitSystem);
        
        ColorSchemeByType colorScheme = (ColorSchemeByType)getDisplayBox(sim.box).getColorScheme();
        colorScheme.setColor(sim.leafType, Color.RED);
        colorScheme.setColor(sim.headType, new Color(190, 0, 190));
        colorScheme.setColor(sim.tailType, Color.BLUE);

        final double expansionFac = 3;
        isExpanded = false;

        final DeviceButton expandButton = new DeviceButton(sim.getController());
        IAction expandAction = new IAction() {
            public void actionPerformed() {
                dim.E(sim.box.getBoundary().getBoxSize());
                dim.setX(0, expansionFac * dim.getX(0));
                sim.box.getBoundary().setBoxSize(dim);
                
                int numSurfactants = (int)nSurfactantSlider.getValue();
                Box pretendBox = new Box(space);
                dim.setX(0, dim.getX(0)/expansionFac);
                pretendBox.getBoundary().setBoxSize(dim);
                sim.addBox(pretendBox);
                pretendBox.setNMolecules(sim.surfactant, numSurfactants);
                configLattice.initializeCoordinates(pretendBox);
                IMoleculeList surfactants = pretendBox.getMoleculeList(sim.surfactant);
                for (int i=0; i<numSurfactants; i++) {
                    IMolecule surfactant = surfactants.getMolecule(0);
                    pretendBox.removeMolecule(surfactant);
                    double deltaX = 0.55 * dim.getX(0);
                    if (surfactant.getChildList().getAtom(0).getPosition().getX(0) < 0) {
                        deltaX = -deltaX;
                    }
                    for (int j=0; j<2; j++) {
                        Vector pos = surfactant.getChildList().getAtom(j).getPosition();
                        pos.setX(0, pos.getX(0) + deltaX);
                    }
                    sim.box.addMolecule(surfactant);
                }
                sim.removeBox(pretendBox);
                
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).getNeighborManager(sim.box).reset();
                sim.integrator.reset();

                if (numSurfactants > 0) {
                    sim.integrator.getEventManager().addListener(surfactantProfilePump);
                    sim.integrator.getEventManager().addListener(orientationProfilePump);
                }

                nSlider.setEnabled(false);
                nSurfactantSlider.getSlider().setEnabled(false);
                expandButton.getButton().setEnabled(false);
                xSlider.getSlider().setEnabled(false);
                getDisplayBox(sim.box).repaint();
                densityProfileMeter.reset();
                surfactantProfileMeter.reset();
                orientationProfileMeter.reset();
                muProfileMeter.reset();
                
                isExpanded = true;
            }
            Vector dim = space.makeVector();
            ConfigurationLattice configLattice = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        };
        expandButton.setAction(expandAction);
        expandButton.setLabel("Expand");
        ((JPanel)getController().graphic()).add(expandButton.graphic());

        // the reinitialize's preAction halts the integrator 
        final IAction oldPreAction = getController().getReinitButton().getPreAction();
        getController().getReinitButton().setPreAction(new IAction() {
            public void actionPerformed() {
                oldPreAction.actionPerformed();
                if (!isExpanded) return;
                Vector dim = space.makeVector();
                dim.E(sim.box.getBoundary().getBoxSize());
                dim.setX(0, dim.getX(0) / expansionFac);
                sim.box.setNMolecules(sim.surfactant, 0);
                int nMolecules = sim.box.getNMolecules(sim.species);
                // need to delete the molecules so that the NeighborCellManager doesn't try
                // to give them a new cell (which would fail)
                sim.box.setNMolecules(sim.species, 0);
                sim.box.getBoundary().setBoxSize(dim);
                sim.box.setNMolecules(sim.species, nMolecules);
                nSlider.setEnabled(true);
                nSurfactantSlider.getSlider().setEnabled(true);
                expandButton.getButton().setEnabled(true);
                xSlider.getSlider().setEnabled(true);
                isExpanded = false;

                densityProfileMeter.reset();
                surfactantProfileMeter.reset();
                orientationProfileMeter.reset();
                muProfileMeter.reset();
                
                if ((int)nSurfactantSlider.getValue() > 0) {
                    sim.integrator.getEventManager().removeListener(surfactantProfilePump);
                    sim.integrator.getEventManager().removeListener(orientationProfilePump);
                    profilePlot.getDataSet().reset();
                    profilePlot.setDoLegend(false);
                    profilePlot.setDoLegend(true);
                    profilePlot.doUpdate();
                    orientationPlot.getDataSet().reset();
                    orientationPlot.doUpdate();
                }
            }
        });
        
        IAction recenterAction = new IAction() {
            public void actionPerformed() {
                if (!isExpanded) {
                    return;
                }
                double L = sim.box.getBoundary().getBoxSize().getX(0);

                // calculate structure factor and phase angle for lowest-frequency
                // concentration wave (delta rho (x)).
                
                IAtomList leafAtoms = sim.box.getLeafList();
                int nTot = leafAtoms.getAtomCount();
                double sumCos = 0, sumSin = 0;
                double q = 2*Math.PI/L;
                for (int i=0; i<nTot; i++) {
                    Vector pos = leafAtoms.getAtom(i).getPosition();
                    double sinx = Math.sin(q*pos.getX(0));
                    double cosx = Math.cos(q*pos.getX(0));
                    sumCos += cosx;
                    sumSin += sinx;
                }
                double amplitude = (2*Math.sqrt((sumCos*sumCos+sumSin*sumSin))/nTot);
                // concentration wave amplitude must be large enough to correspond
                // to a substantial segregation (phase separation).
                if (amplitude < 0.75) {
//                    System.out.println("not centering "+amplitude);
                    return;
                }
//                System.out.println("centering "+amplitude);
                // phase angle = atan(sumCos/sumSin), where delta rho = 0
                double center = -(Math.atan2(sumCos, sumSin)/q-0.25*L);
                if (center > 1) {
                    center = 1;
                }
                else if (center < -1) {
                    center = -1;
                }
                for (int i=0; i<nTot; i++) {
                    Vector pos = leafAtoms.getAtom(i).getPosition();
                    pos.setX(0, pos.getX(0) - center);
                }
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).getNeighborManager(sim.box).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    // we can cause overlap by increasing tail diameter
                    // if so, that might not have been resolved yet
                }
            }
        };
        IntegratorListenerAction recenterActionListener = new IntegratorListenerAction(recenterAction);
        sim.integrator.getEventManager().addListener(recenterActionListener);
        recenterActionListener.setInterval(100);

        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
        DisplayTimer displayTimer = new DisplayTimer(sim.integrator);
        add(displayTimer);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature(sim.box, space.D());
        DataFork temperatureFork = new DataFork();
        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
        sim.integrator.getEventManager().addListener(temperaturePumpListener);
        temperaturePumpListener.setInterval(20);
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
	    
//	      MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential(), sim.box);
//        AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        energyHistory.setTimeDataSource(timeCounter);
//        DataPump energyPump = new DataPump(eMeter, energyHistory);
//        sim.integrator.addIntervalAction(energyPump);
//        sim.integrator.setActionInterval(energyPump, 60);
//        energyHistory.setPushInterval(5);
//        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(10);
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        sim.integrator.getEventManager().addListener(pePumpListener);
        pePumpListener.setInterval(10);
        peHistory.setPushInterval(1);
        dataStreamPumps.add(pePump);
		
        DisplayPlot ePlot = new DisplayPlot();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());

        ePlot.getPlot().setTitle("Potential Energy History");
		ePlot.setDoLegend(false);
		ePlot.setLabel("Energy");
		
        DataSourceTensorVirialHardProfile pMeter = new DataSourceTensorVirialHardProfile(space);
        pMeter.setIntegrator(sim.integrator);
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
        IntegratorListenerAction pPumpListener = new IntegratorListenerAction(pPump);
        sim.integrator.getEventManager().addListener(pPumpListener);
        pPumpListener.setInterval(20);
        dataStreamPumps.add(pPump);

        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAccumulator);

        DataProcessorInterfacialTension interfacialTension = new DataProcessorInterfacialTension(space);
        interfacialTension.setBox(sim.box);
        virialFork.addDataSink(interfacialTension);
        final AccumulatorAverageCollapsing tensionAvg = new AccumulatorAverageCollapsing();
        interfacialTension.setDataSink(tensionAvg);
        tensionAvg.setPushInterval(10);
        DisplayTextBoxesCAE tensionDisplay = new DisplayTextBoxesCAE();
        tensionDisplay.setAccumulator(tensionAvg);

        DataSourceVirialProfile virialProfileMeter = new DataSourceVirialProfile(pMeter);
        DataFork virialProfileFork = new DataFork();
        DataPump virialProfilePump = new DataPump(virialProfileMeter, virialProfileFork);
        DataGroupSplitter virialSplitter = new DataGroupSplitter();
        virialProfileFork.addDataSink(virialSplitter);
        IntegratorListenerAction virialProfilePumpListener = new IntegratorListenerAction(virialProfilePump);
        sim.integrator.getEventManager().addListener(virialProfilePumpListener);
        virialProfilePumpListener.setInterval(20);
        AccumulatorAverageFixed[] virialProfileAvg = new AccumulatorAverageFixed[space.D()];
        DisplayPlot virialPlot = new DisplayPlot();
        for (int i=0; i<space.D(); i++) {
            virialProfileAvg[i] = new AccumulatorAverageFixed(10);
            virialProfileAvg[i].setPushInterval(10);
            virialSplitter.setDataSink(i, virialProfileAvg[i]);
            virialProfileAvg[i].addDataSink(virialPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{virialProfileAvg[i].AVERAGE});
            virialPlot.setLegend(new DataTag[]{virialProfileAvg[i].getTag()}, comp[i]+" Virial");
        }
        virialPlot.setLabel("Virial Profile");
        add(virialPlot);
        dataStreamPumps.add(virialProfilePump);

        DataProcessorInterfacialTensionProfile interfacialTensionProfile = new DataProcessorInterfacialTensionProfile(space);
        interfacialTensionProfile.setBox(sim.box);
        virialProfileFork.addDataSink(interfacialTensionProfile);
        AccumulatorAverageFixed tensionProfileAvg = new AccumulatorAverageFixed(10);
        interfacialTensionProfile.setDataSink(tensionProfileAvg);
        tensionProfileAvg.setPushInterval(10);
        DisplayPlot tensionPlot = new DisplayPlot();
        tensionPlot.setDoLegend(false);
        tensionProfileAvg.addDataSink(tensionPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{tensionProfileAvg.AVERAGE});
        tensionPlot.setLabel("Tension Profile");
        add(tensionPlot);

        densityProfileMeter = new MeterProfileByVolume(space);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(10);
        DataPump profilePump = new DataPump(densityProfileMeter, densityProfileAvg);
        DataDump profileDump = new DataDump();
        densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        IntegratorListenerAction profilePumpListener = new IntegratorListenerAction(profilePump);
        sim.integrator.getEventManager().addListener(profilePumpListener);
        profilePumpListener.setInterval(10);
        dataStreamPumps.add(profilePump);

        MoleculePositionGeometricCenterPBC positionDefinitionPBC = new MoleculePositionGeometricCenterPBC(space, sim.box.getBoundary());
        surfactantProfileMeter = new MeterProfileByVolume(space);
        surfactantProfileMeter.setBox(sim.box);
        surfactantProfileMeter.setPositionDefinition(positionDefinitionPBC);
        meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.surfactant);
        surfactantProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed surfactantProfileAvg = new AccumulatorAverageFixed(10);
        surfactantProfileAvg.setPushInterval(10);
        surfactantProfilePump = new DataPumpListener(surfactantProfileMeter, surfactantProfileAvg, 10);
        dataStreamPumps.add(surfactantProfilePump);
        
        orientationProfileMeter = new MeterProfileByAtoms(space);
        orientationProfileMeter.setBox(sim.box);
        orientationProfileMeter.setSpecies(sim.surfactant);
        orientationProfileMeter.setPositionDefinition(positionDefinitionPBC);
        MeterOrientation meterOrientation = new MeterOrientation(space);
        meterOrientation.setBox(sim.box);
        orientationProfileMeter.setDataSource(meterOrientation);
        orientationPlot = new DisplayPlot();
        orientationPlot.setDoLegend(false);
        orientationPlot.setLabel("Orientation");
        orientationPlot.getPlot().setTitle("Surfactant Orientation");
        add(orientationPlot);
        orientationProfilePump = new DataPumpListener(orientationProfileMeter, orientationPlot.getDataSet().makeDataSink(), 10);
        dataStreamPumps.add(orientationProfilePump);
        
        final FitTanh fitTanh = new FitTanh();
        densityProfileAvg.addDataSink(fitTanh, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});

        profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        surfactantProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{surfactantProfileAvg.AVERAGE});
        fitTanh.setDataSink(profilePlot.getDataSet().makeDataSink());
        profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag()}, "density");
        profilePlot.setLegend(new DataTag[]{surfactantProfileAvg.getTag()}, "surfactant");
        profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag(), fitTanh.getTag()}, "fit");
        profilePlot.setDoLegend(true);
        profilePlot.setLabel("Density");
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
        JPanel profilePanel = new JPanel(new GridBagLayout());
        getPanel().tabbedPane.add("Density", profilePanel);
        profilePanel.add(profilePlot.graphic(), vertGBC);
        final DisplayTextBox vaporDensityBox = new DisplayTextBox();
        vaporDensityBox.setLabel("Vapor density");
        final DisplayTextBox liquidDensityBox = new DisplayTextBox();
        liquidDensityBox.setLabel("Vapor density");
        final DisplayTextBox interfaceWidthBox = new DisplayTextBox();
        interfaceWidthBox.setLabel("Vapor density");
        final DisplayTextBox interfaceLocationBox = new DisplayTextBox();
        interfaceLocationBox.setLabel("Vapor density");
        JPanel fitParamsPanel = new JPanel();
        fitParamsPanel.add(vaporDensityBox.graphic());
        fitParamsPanel.add(liquidDensityBox.graphic());
        fitParamsPanel.add(interfaceWidthBox.graphic());
        fitParamsPanel.add(interfaceLocationBox.graphic());
        profilePanel.add(fitParamsPanel, vertGBC);
        IAction pullParams = new IAction() {
            public void actionPerformed() {
                double[] params = fitTanh.getLastBestParam();
                data.x = params[0];
                vaporDensityBox.putData(data);
                data.x = params[1];
                liquidDensityBox.putData(data);
                data.x = params[2];
                interfaceWidthBox.putData(data);
                data.x = params[3];
                interfaceLocationBox.putData(data);
            }
            DataDouble data = new DataDouble();
        };
        IntegratorListenerAction pullParamsListener = new IntegratorListenerAction(pullParams);
        sim.integrator.getEventManager().addListener(pullParamsListener);
        pullParamsListener.setInterval(10);
        
        DisplayPlot muPlot = new DisplayPlot();
        muProfileMeter = new MeterProfile(space, sim.getRandom());
        muProfileMeter.setBox(sim.box);
        DataSourcePositionedBoltzmannFactor meterChemicalPotential = new DataSourcePositionedBoltzmannFactor(space);
        meterChemicalPotential.setIntegrator(sim.integrator);
        meterChemicalPotential.setSpecies(sim.species);
        muProfileMeter.setDataSource(meterChemicalPotential);
        AccumulatorAverageFixed chemicalPotentialAverage = new AccumulatorAverageFixed(10);
        chemicalPotentialAverage.setPushInterval(10);
        DataPump muProfilePump = new DataPump(muProfileMeter, chemicalPotentialAverage);
        DataProcessorChemicalPotential dataProcessorChemicalPotential = new DataProcessorChemicalPotential();
        dataProcessorChemicalPotential.setDensityProfileDump(profileDump);
        dataProcessorChemicalPotential.setIntegrator(sim.integrator);
        chemicalPotentialAverage.addDataSink(dataProcessorChemicalPotential, new AccumulatorAverage.StatType[]{chemicalPotentialAverage.AVERAGE});
        dataProcessorChemicalPotential.setDataSink(muPlot.getDataSet().makeDataSink());
        muPlot.setLegend(new DataTag[]{dataProcessorChemicalPotential.getTag()}, "mu");
//        dataProcessorChemicalPotential = new DataProcessorChemicalPotential();
//        dataProcessorChemicalPotential.doE = false;
//        dataProcessorChemicalPotential.setDensityProfileDump(profileDump);
//        dataProcessorChemicalPotential.setIntegrator(sim.integrator);
//        chemicalPotentialAverage.addDataSink(dataProcessorChemicalPotential, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
//        dataProcessorChemicalPotential.setDataSink(muPlot.getDataSet().makeDataSink());
//        muPlot.setLegend(new DataTag[]{dataProcessorChemicalPotential.getTag()}, "muD");
//        dataProcessorChemicalPotential = new DataProcessorChemicalPotential();
//        dataProcessorChemicalPotential.doD = false;
//        dataProcessorChemicalPotential.setDensityProfileDump(profileDump);
//        dataProcessorChemicalPotential.setIntegrator(sim.integrator);
//        chemicalPotentialAverage.addDataSink(dataProcessorChemicalPotential, new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
//        dataProcessorChemicalPotential.setDataSink(muPlot.getDataSet().makeDataSink());
//        muPlot.setLegend(new DataTag[]{dataProcessorChemicalPotential.getTag()}, "muE");
        muPlot.setLabel("Chemical Potential");
        muPlot.setDoLegend(false);
        add(muPlot);
        IntegratorListenerAction muProfilePumpListener = new IntegratorListenerAction(muProfilePump);
        sim.integrator.getEventManager().addListener(muProfilePumpListener);
        muProfilePumpListener.setInterval(50);
        dataStreamPumps.add(muProfilePump);

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
        final ConfigurationLattice config = new ConfigurationLattice((space.D() == 2) ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space);
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                int n = sim.box.getNMolecules(sim.species);
                if(n == 0) {
                	sim.integrator.setThermostatInterval(800);
                }
                else {
                    sim.integrator.setThermostatInterval((800+(n-1))/n);
                }
                
                config.initializeCoordinates(sim.box);
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).getNeighborManager(sim.box).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    if (n == oldN) {
                        throw new RuntimeException(e);
                    }
                    nSlider.setValue(oldN);
                    return;
                }
                oldN = n;

                getController().getSimRestart().actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
            
            int oldN = sim.box.getMoleculeList().getMoleculeCount();
        });

        nSurfactantSlider = new DeviceSlider(null);
        nSurfactantSlider.setMinimum(0);
        nSurfactantSlider.setMaximum(50);
        nSurfactantSlider.setValue(0);
        nSurfactantSlider.setLabel("Number of Surfactants");
        nSurfactantSlider.setShowBorder(true);
        nSurfactantSlider.setShowValues(true);

        DeviceSlider surfactantEpsilon = new DeviceSlider(sim.getController());
        surfactantEpsilon.setShowBorder(true);
        surfactantEpsilon.setModifier(new Modifier() {

            public Dimension getDimension() {
                return Energy.DIMENSION;
            }

            public String getLabel() {
                return "Surfactant well depth";
            }

            public double getValue() {
                return sim.p2HeadHead.getEpsilon();
            }

            public void setValue(double newValue) {
                if (newValue < 0 || newValue > 10) {
                    throw new RuntimeException("Illegal well depth");
                }
                sim.p2HeadHead.setEpsilon(newValue);
                sim.p2Head.setEpsilon(Math.sqrt(newValue));
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
            }
        });
        surfactantEpsilon.setMaximum(5);
        surfactantEpsilon.setShowValues(true);
        surfactantEpsilon.setLabel("head epsilon");
        
        DeviceSlider surfactantSigma = new DeviceSlider(sim.getController());
        surfactantSigma.setShowBorder(true);
        surfactantSigma.setModifier(new Modifier() {

            public Dimension getDimension() {
                return Length.DIMENSION;
            }

            public String getLabel() {
                return "Surfactant tail diameter";
            }

            public double getValue() {
                return sim.p2TailTail.getCollisionDiameter();
            }

            public void setValue(double newValue) {
                if (newValue < 0 || newValue > 10) {
                    throw new RuntimeException("Illegal well depth");
                }
                sim.p2TailTail.setCollisionDiameter(newValue);
                sim.p2Tail.setCollisionDiameter(0.5+0.5*newValue);
                sim.p2HeadTail.setCollisionDiameter(0.5+0.5*newValue);
                sim.p2Bond.setBondLength(0.5+0.5*newValue-0.2);
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    // ignore... if the diameter is too large we can get overlap.  they'll fix themselves over time
                }
                ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.tailType, newValue);
            }
        });
        surfactantSigma.setLabel("tail diameter");
        surfactantSigma.setPrecision(1);
        surfactantSigma.setMinimum(1);
        surfactantSigma.setMaximum(1.6);
        surfactantSigma.setNMajor(3);
        surfactantSigma.setShowValues(true);
        
        IAction reconfig = new IAction() {
            public void actionPerformed() {
                config.initializeCoordinates(sim.box);
                ((PotentialMasterList)sim.integrator.getPotentialMaster()).getNeighborManager(sim.box).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
                getController().getSimRestart().actionPerformed();
                getDisplayBox(sim.box).repaint();
                densityProfileMeter.reset();
                surfactantProfileMeter.reset();
                orientationProfileMeter.reset();
                muProfileMeter.reset();
            }
        };
        xSlider = new DeviceSlider(sim.getController());
        xSlider.setShowValues(true);
        xSlider.setLabel("Box length");
        xSlider.setShowBorder(true);
        xSlider.setMinimum(6);
        xSlider.setMaximum(30);
        xSlider.setModifier(new ModifierBoxSize(space, sim.box, 0, reconfig));
        JPanel systemPanel = new JPanel(new GridBagLayout());
        systemPanel.add(xSlider.graphic(), vertGBC);
        
        //************* Lay out components ****************//


        getDisplayBox(sim.box).setScale(0.7);

        //temperature selector
        temperatureSelect = new DeviceThermoSlider(sim.getController(), sim.integrator);
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
        		
        		orientationProfileMeter.reset();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        systemPanel.add(temperatureSelect.graphic(), vertGBC);
        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.add("System", systemPanel);
        getPanel().controlPanel.add(tabbedPane, vertGBC);
        JPanel numMoleculesPanel = new JPanel(new GridBagLayout());
        numMoleculesPanel.add(nSlider.graphic(), vertGBC);
        numMoleculesPanel.add(nSurfactantSlider.graphic(), vertGBC);
        tabbedPane.add("# of molecules", numMoleculesPanel);
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(surfactantEpsilon.graphic(), vertGBC);
        potentialPanel.add(surfactantSigma.graphic(), vertGBC);
        tabbedPane.add("Surfactant potential", potentialPanel);

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

        InterfacialSW sim = new InterfacialSW(sp);
        InterfacialSWGraphic swGraphic = new InterfacialSWGraphic(sim, sp);
        swGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(9));
		SimulationGraphic.makeAndDisplayFrame
		        (swGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
	        Space sp = Space3D.getInstance();
	        InterfacialSW sim = new InterfacialSW(sp);
            InterfacialSWGraphic ljmdGraphic = new InterfacialSWGraphic(sim, sp);
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
        
        protected IData processData(IData inputData) {
            double[] x = data.getData();
            for (int i=0; i<x.length; i++) {
                x[i] = ((DataTensor)inputData).x.component(i,i);
            }
            return data;
        }

        protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
            dataInfo = new DataDoubleArray.DataInfoDoubleArray(inputDataInfo.getLabel(), inputDataInfo.getDimension(), new int[]{inputDataInfo.getLength()});
            return dataInfo;
        }

        public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
            if (!(inputDataInfo instanceof DataTensor.DataInfoTensor)) {
                throw new IllegalArgumentException("Gotta be a DataInfoTensor");
            }
            return null;
        }

        private static final long serialVersionUID = 1L;
        protected final DataDoubleArray data;
    }
    
    public static class ModifierBoxSize implements Modifier {
        public ModifierBoxSize(Space space, Box box, int dim, IAction reconfig) {
            this.box = box;
            this.dim = dim;
            this.reconfig = reconfig;
            size = space.makeVector();
        }
        
        public Dimension getDimension() {
            return Length.DIMENSION;
        }

        public String getLabel() {
            return "Box size";
        }

        public double getValue() {
            return box.getBoundary().getBoxSize().getX(dim);
        }

        public void setValue(double newValue) {
            if (newValue <= 0) {
                throw new IllegalArgumentException("Gotta be positive");
            }
            //newValue+=0.01;
            size.E(box.getBoundary().getBoxSize());
            double oldValue = size.getX(dim);
            size.setX(dim, newValue);
            if (dim == 1 && size.getD() == 3) {
                size.setX(2, newValue);
            }
            box.getBoundary().setBoxSize(size);
            try {
                reconfig.actionPerformed();
            }
            catch (RuntimeException e) {
                // box is too small.  restore to original size
                size.setX(dim, oldValue);
                box.getBoundary().setBoxSize(size);
                // and reconfig.  this shouldn't throw.
                reconfig.actionPerformed();
            }
        }
        
        protected final Box box;
        protected final int dim;
        protected final IAction reconfig;
        protected final Vector size;
    }
}


