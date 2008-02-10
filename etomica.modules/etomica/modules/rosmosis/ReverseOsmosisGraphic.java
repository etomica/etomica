package etomica.modules.rosmosis;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;

import etomica.action.Action;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.Data;
import etomica.data.DataFork;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterProfile;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.potential.P2LennardJones;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Angstrom;
import etomica.units.Bar;
import etomica.units.CompoundUnit;
import etomica.units.Dalton;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Length;
import etomica.units.Liter;
import etomica.units.Meter;
import etomica.units.Mole;
import etomica.units.Picosecond;
import etomica.units.Pixel;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.systems.MKS;
import etomica.util.Constants.CompassDirection;

/**
 * Graphical and data collection bits for reverse osmosis simulation.
 * 
 * @author Andrew Schultz
 */
public class ReverseOsmosisGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Square-Well Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 2;
    protected DeviceThermoSlider tempSlider;
    protected DeviceSlider solventChamberDensitySlider, soluteChamberDensitySlider;
    protected DeviceSlider soluteMoleFractionSlider;
    protected DeviceBox sigBox, epsBox, massBox, tetherBox;
    protected DeviceSlider membraneThicknessSlider;
    protected double epsilonSolute, massSolute, sigmaSolute;
    protected double epsilonSolvent, massSolvent, sigmaSolvent;
    protected double epsilonMembrane, massMembrane, sigmaMembrane;
    protected Unit eUnit, dUnit, pUnit;
    protected ReverseOsmosis sim;
    
    public ReverseOsmosisGraphic(final ReverseOsmosis simulation) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        ArrayList dataStreamPumps = getController().getDataStreamPumps();

    	this.sim = simulation;

        Unit tUnit = Kelvin.UNIT;

        eUnit = new UnitRatio(Joule.UNIT, Mole.UNIT);
        epsilonSolute = sim.potential11.getEpsilon();
        massSolute = sim.speciesSolute.getLeafType().getMass();
        sigmaSolute = ((AtomTypeSphere)sim.speciesSolute.getLeafType()).getDiameter();
        epsilonSolvent = sim.potential22.getEpsilon();
        massSolvent = sim.speciesSolvent.getLeafType().getMass();
        sigmaSolvent = ((AtomTypeSphere)sim.speciesSolvent.getLeafType()).getDiameter();
        epsilonMembrane = sim.potentialMM.getEpsilon();
        massMembrane = sim.speciesMembrane.getLeafType().getMass();
        sigmaMembrane = ((AtomTypeSphere)sim.speciesMembrane.getLeafType()).getDiameter();

        if (sim.getSpace().D() == 2) {
            dUnit = new UnitRatio(Mole.UNIT, 
                                    new MKS().area());
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        }
        else {
            dUnit = new UnitRatio(Mole.UNIT, Liter.UNIT);
            pUnit = Bar.UNIT;

        }
        

        if (sim.getSpace().D() == 2) {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(400/sim.box.getBoundary().getDimensions().x(1)));
        }
        else {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getDimensions().x(1)));
        }

        sim.activityIntegrate.setSleepPeriod(0);
       
        //combo box to select potentials
        sigBox = new DeviceBox();
        epsBox = new DeviceBox();
        massBox = new DeviceBox();
        tetherBox = new DeviceBox();
        membraneThicknessSlider = new DeviceSlider(sim.getController());

        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
        displayCycles.setPrecision(6);
        DataPump pump = new DataPump(meterCycles,displayCycles);
        sim.integrator.addIntervalAction(pump);
        displayCycles.setLabel("Simulation time");
        
        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController());
        tempSlider.setUnit(Kelvin.UNIT);
//        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(500.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setUnit(tUnit);
        tempSlider.setAdiabatic();
        tempSlider.setIntegrator(sim.integrator);
        
        ModifierGeneral modifier = new ModifierGeneral(sim.configMembrane, "solventChamberDensity");
        solventChamberDensitySlider = new DeviceSlider(sim.getController(), modifier);
        solventChamberDensitySlider.setPrecision(1);
        solventChamberDensitySlider.setMaximum(30);
        solventChamberDensitySlider.setNMajor(3);
//        solventChamberDensitySlider.setValue(sim.configMembrane.getSolventChamberDensity());
        solventChamberDensitySlider.setShowValues(true);
        solventChamberDensitySlider.setUnit(dUnit);
        solventChamberDensitySlider.doUpdate();
        solventChamberDensitySlider.setShowBorder(true);
        solventChamberDensitySlider.setLabel("Solvent Density (mol/L)");
        
        modifier = new ModifierGeneral(sim.configMembrane, "solutionChamberDensity");
        soluteChamberDensitySlider = new DeviceSlider(sim.getController(), modifier);
        soluteChamberDensitySlider.setPrecision(1);
        soluteChamberDensitySlider.setMaximum(30);
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
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);
        
        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.setBorder(new TitledBorder(null, "Potential Selection", TitledBorder.CENTER, TitledBorder.TOP));
        JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(sigBox.graphic());
        parameterPanel.add(epsBox.graphic());
        parameterPanel.add(massBox.graphic());
        parameterPanel.add(tetherBox.graphic());
        potentialPanel.add(parameterPanel,vertGBC);
        potentialPanel.add(membraneThicknessSlider.graphic(), vertGBC);

        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Membrane");
        setupPanel.add(configPanel, "Configuration");


        ModifierAtomDiameter sigMembraneModifier = new ModifierAtomDiameter(sim.speciesMembrane, sim.potentialMM,
                new SpeciesSpheresMono[]{sim.speciesSolute, sim.speciesSolvent}, new P2LennardJones[]{sim.potentialM1, sim.potentialM2});
        ModifierEpsilon epsModifier = new ModifierEpsilon(sim.potentialMM, new P2LennardJones[]{sim.potential11, sim.potential22},
                new P2LennardJones[]{sim.potentialM1, sim.potentialM2});
        ModifierGeneral massModifier = new ModifierGeneral(sim.speciesSolute.getLeafType().getElement(),"mass");
        ModifierGeneral tetherModifier = new ModifierGeneral(sim.potentialTether,"epsilon");
        Modifier membraneThicknessModifier = new Modifier() {
            public Dimension getDimension() { return Quantity.DIMENSION; }
            public String getLabel() { return "membrane thickness"; }
            public double getValue() { return sim.configMembrane.getNumMembraneLayers(); }
            public void setValue(double newValue) { sim.configMembrane.setNumMembraneLayers((int)newValue); }
        };
        sigBox.setModifier(sigMembraneModifier);
        sigBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
        sigBox.doUpdate();
        epsBox.setUnit(eUnit);
        epsBox.setModifier(epsModifier);
        massBox.setModifier(massModifier);
        massBox.setUnit(Dalton.UNIT);
        tetherBox.setUnit(eUnit);
        tetherBox.setModifier(tetherModifier);
        tetherBox.setLabel("Tether Constant ("+eUnit.symbol()+")");
        sigBox.setController(sim.getController());
        epsBox.setController(sim.getController());
        massBox.setController(sim.getController());
        tetherBox.setController(sim.getController());
        membraneThicknessSlider.setLabel("Membrane Thickness");
        membraneThicknessSlider.setShowBorder(true);
        membraneThicknessSlider.setMinimum(1);
        membraneThicknessSlider.setMaximum(4);
        membraneThicknessSlider.setModifier(membraneThicknessModifier);
        
        Action neighborRangeReset = new Action() {
            public void actionPerformed() {
//                ((PotentialMasterList)sim.integrator.getPotential()).reset();
//                double nbrRange = ((PotentialMasterList)sim.integrator.getPotential()).getMaxPotentialRange();
//                nbrRange *= 1.2;
//                ((PotentialMasterList)sim.integrator.getPotential()).setRange(nbrRange);
//                ((PotentialMasterList)sim.integrator.getPotential()).reset();
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {}
                getDisplayBox(sim.box).repaint();
            }
        };
        sigBox.setPostAction(neighborRangeReset);

        //display of box, timer
        ColorSchemeByType colorScheme = (ColorSchemeByType)getDisplayBox(sim.box).getColorScheme();
        colorScheme.setColor(sim.speciesSolute.getLeafType(),java.awt.Color.RED);
        colorScheme.setColor(sim.speciesSolvent.getLeafType(),java.awt.Color.BLUE);
        colorScheme.setColor(sim.speciesMembrane.getLeafType(),java.awt.Color.WHITE);

	    //meters and displays
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);

        //add meter and display for current kinetic temperature

		MeterTemperature thermometer = new MeterTemperature();
        thermometer.setBox(sim.box);
        DataFork temperatureFork = new DataFork();
        final DisplayTextBox tBox = new DisplayTextBox();
        final DataPump temperaturePump = new DataPump(thermometer, tBox);
        sim.integrator.addIntervalAction(temperaturePump);
        sim.integrator.setActionInterval(temperaturePump, 10);
		dataStreamPumps.add(temperaturePump);
        tBox.setUnit(tUnit);
		tBox.setLabel("Measured Temperature");
		tBox.setLabelPosition(CompassDirection.NORTH);

		MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential());
        eMeter.setBox(sim.box);
        final AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setTimeDataSource(timeCounter);
        DataSinkExcludeOverlap foo = new DataSinkExcludeOverlap(sim.box);
        final DataPump energyPump = new DataPump(eMeter, foo);
        foo.setDataSink(energyHistory);
        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.setActionInterval(energyPump, 10);
        dataStreamPumps.add(energyPump);
		
		MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.integrator.getPotential());
		peMeter.setBox(sim.box);
        final AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setTimeDataSource(timeCounter);
        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
        peAccumulator.setPushInterval(2);
        foo = new DataSinkExcludeOverlap(sim.box);
        final DataPump pePump = new DataPump(peMeter, foo);
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        foo.setDataSink(peFork);
        sim.integrator.addIntervalAction(pePump);
        sim.integrator.setActionInterval(pePump, 10);
        dataStreamPumps.add(pePump);

		MeterKineticEnergy keMeter = new MeterKineticEnergy();
		keMeter.setBox(sim.box);
        final AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setTimeDataSource(timeCounter);
        // we do this for the scaling by numAtoms rather than for the overlap exclusion
        foo = new DataSinkExcludeOverlap(sim.box);
        final DataPump kePump = new DataPump(keMeter, foo);
        foo.setDataSink(keHistory);
        sim.integrator.addIntervalAction(kePump);
        sim.integrator.setActionInterval(kePump, 10);
        dataStreamPumps.add(kePump);

        MeterFlux meterFlux = new MeterFlux();
        double xLength = sim.box.getBoundary().getDimensions().x(0);
        meterFlux.setBoundaries(0, new double[]{-0.25*xLength, 0.25*xLength}, new int[]{1, -1});
        meterFlux.setIntegrator(sim.integrator);
        meterFlux.setBox(sim.box);
        meterFlux.setSpecies(new ISpecies[]{sim.speciesSolute, sim.speciesSolvent});
        AccumulatorAverageCollapsing fluxAvg = new AccumulatorAverageCollapsing();
        fluxAvg.setPushInterval(10);
        DataPump fluxPump = new DataPump(meterFlux, fluxAvg);
        sim.integrator.addIntervalAction(fluxPump);
        // has to happen before PBC are applied
        sim.integrator.setIntervalActionPriority(fluxPump, 1);
        dataStreamPumps.add(fluxPump);

        MeterNMolecules densitySolute = new MeterNMolecules();
        MeterNMolecules densitySolvent = new MeterNMolecules();
        densitySolute.setBox(sim.box);
        densitySolvent.setBox(sim.box);
        densitySolute.setSpecies(sim.speciesSolute);
        densitySolvent.setSpecies(sim.speciesSolvent);
        MeterProfile meterProfileSolute = new MeterProfile(sim.getSpace());
        MeterProfile meterProfileSolvent = new MeterProfile(sim.getSpace());
        meterProfileSolute.setBox(sim.box);
        meterProfileSolute.setDataSource(densitySolute);
        meterProfileSolvent.setBox(sim.box);
        meterProfileSolvent.setDataSource(densitySolvent);
        IVector profVector = sim.getSpace().makeVector();
        profVector.setX(0, 1);
        meterProfileSolute.setProfileVector(profVector);
        meterProfileSolvent.setProfileVector(profVector);
        AccumulatorAverageFixed profSoluteAvg = new AccumulatorAverageFixed(10);
        profSoluteAvg.setPushInterval(10);
        DataPump profPumpSolute = new DataPump(meterProfileSolute, profSoluteAvg);
        sim.integrator.addIntervalAction(profPumpSolute);
        sim.integrator.setActionInterval(profPumpSolute, 2);
        AccumulatorAverageFixed profSolventAvg = new AccumulatorAverageFixed(10);
        profSolventAvg.setPushInterval(10);
        DataPump profPumpSolvent = new DataPump(meterProfileSolvent, profSolventAvg);
        sim.integrator.addIntervalAction(profPumpSolvent);
        sim.integrator.setActionInterval(profPumpSolvent, 2);
        dataStreamPumps.add(profPumpSolute);
        dataStreamPumps.add(profPumpSolvent);

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
		ePlot.setXUnit(Picosecond.UNIT);

        final DisplayPlot profPlot = new DisplayPlot();
        profSoluteAvg.addDataSink(profPlot.getDataSet().makeDataSink(),
                new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        profPlot.setLegend(new DataTag[]{meterProfileSolute.getTag()}, "Solute");
        profSolventAvg.addDataSink(profPlot.getDataSet().makeDataSink(),
                new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        profPlot.setLegend(new DataTag[]{meterProfileSolvent.getTag()}, "Solvent");
        profPlot.getPlot().setTitle("DensityProfiles (mol/L)");
        profPlot.setDoLegend(true);
        profPlot.setLabel("Profiles");
        profPlot.setXUnit(Angstrom.UNIT);
        profPlot.setUnit(dUnit);

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

        DisplayTextBoxesCAE fluxDisplay = new DisplayTextBoxesCAE();
        fluxDisplay.setAccumulator(fluxAvg);
        fluxDisplay.setLabel("Flux (Atoms/ps/A^2)");

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


        final Action resetAction = new Action() {
        	public void actionPerformed() {
        	    try {
        	        sim.integrator.reset();
        	    }
        	    catch (ConfigurationOverlapException e){}

        		// Reset temperature (THIS IS NOT WORKING)
                temperaturePump.actionPerformed();
                tBox.repaint();

                // IS THIS WORKING?
//                pPump.actionPerformed();
//                pDisplay.putData(pAccumulator.getData());
//                pDisplay.repaint();
                peDisplay.putData(peAccumulator.getData());
                peDisplay.repaint();

        		getDisplayBox(sim.box).graphic().repaint();
        		
        		displayCycles.putData(meterCycles.getData());
        		displayCycles.repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);
        Action reconfigAction = new Action() {
            public void actionPerformed() {
                sim.configMembrane.initializeCoordinates(sim.box);
                resetAction.actionPerformed();
            }
        };
        solventChamberDensitySlider.setPostAction(reconfigAction);
        soluteChamberDensitySlider.setPostAction(reconfigAction);
        soluteMoleFractionSlider.setPostAction(reconfigAction);
        membraneThicknessSlider.setPostAction(reconfigAction);

        getPanel().controlPanel.add(setupPanel, vertGBC);

    	add(ePlot);
        add(profPlot);
    	add(displayCycles);
    	add(tBox);
//    	add(pDisplay);
    	add(peDisplay);
        add(fluxDisplay);
    }

    protected static class ModifierAtomDiameter implements Modifier {
        
        public ModifierAtomDiameter(SpeciesSpheresMono species, P2LennardJones potential,
                SpeciesSpheresMono[] otherSpecies, P2LennardJones[] otherPotentials) {
            this.species = species;
            this.otherSpecies = otherSpecies;
            this.potential = potential;
            this.otherPotentials = otherPotentials;
        }

        public void setValue(double d) {
            if (d > 4.0) {
                throw new IllegalArgumentException("diameter can't exceed 4.0A");
            }
            //assume one type of atom
            ((AtomTypeSphere)species.getLeafType()).setDiameter(d);
            potential.setSigma(d);
            for (int i=0; i<otherPotentials.length; i++) {
                double otherSigma = ((AtomTypeSphere)otherSpecies[i].getLeafType()).getDiameter();
                otherPotentials[i].setSigma(0.5*(d+otherSigma));
            }
//            new BoxImposePbc(sim.box).actionPerformed();
//            if (sim.potentialWrapper.getWrappedPotential() instanceof P2Ideal) {
//                // 0 is OK, but we need to have a sane number of cells
//                ((PotentialMasterList)sim.integrator.getPotential()).setRange(2.0);
//            }
//            else if (sim.potentialWrapper.getWrappedPotential() instanceof P2HardSphere) {
//                ((PotentialMasterList)sim.integrator.getPotential()).setRange(potentialHS.getRange()*1.4);
//            }
//            else { //SW
//                ((PotentialMasterList)sim.integrator.getPotential()).setRange(potentialSW.getRange()*1.2);
//            }
//            ((PotentialMasterList)sim.integrator.getPotential()).reset();
//            try {
//                integrator.reset();
//            }
//            catch (ConfigurationOverlapException e){
                // can happen when increasing diameter
//            }
//            getDisplayBox(integrator.getBox()).repaint();
//            sim.config.setBoundaryPadding(sigma);
        }

        public double getValue() {
            return ((AtomTypeSphere)species.getLeafType()).getDiameter();

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
        protected final SpeciesSpheresMono[] otherSpecies;
        protected final P2LennardJones potential;
        protected final P2LennardJones[] otherPotentials;
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
    
    public static class DataSinkExcludeOverlap extends DataProcessor {

        public DataSinkExcludeOverlap(Box box) {
            myData = new DataDouble();
            this.box = box;
        }
        
        public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
            return null;
        }
        
        public Data processData(Data data) {
            myData.E(data);
            myData.TE(1.0/box.getLeafList().getAtomCount());
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
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }

        ReverseOsmosisGraphic reverseOsmosisGraphic = new ReverseOsmosisGraphic(new ReverseOsmosis(space));
		SimulationGraphic.makeAndDisplayFrame
		        (reverseOsmosisGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            ReverseOsmosisGraphic reverseOsmosisGraphic = new ReverseOsmosisGraphic(new ReverseOsmosis(Space3D.getInstance()));

		    getContentPane().add(reverseOsmosisGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}


