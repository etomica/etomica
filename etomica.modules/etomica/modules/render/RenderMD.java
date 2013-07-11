package etomica.modules.render;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
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
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceBox;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorHard;
import etomica.listener.IntegratorListenerAction;
import etomica.modules.render.ParseObj.BondInfo;
import etomica.modules.swmd.SwmdGraphic.DataSinkExcludeOverlap;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardBondedList;
import etomica.potential.P2Ideal;
import etomica.potential.P2PenetrableSquareWell;
import etomica.potential.PotentialHard;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Picosecond;
import etomica.util.ParameterBase;
import etomica.util.Constants.CompassDirection;

/**
 * 
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.  
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class RenderMD extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    private static final long serialVersionUID = 1L;
    /**
     * The Box holding the atoms. 
     */
    public final IBox box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;

    public final P2HardBondedList potential;
    public final PotentialHard potentialBonded, potentialNonBonded;
    
    public final IPotentialMaster potentialMaster;
    
    public final ParseObj parser;
    public ActivityIntegrate activityIntegrate;

    
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public RenderMD(ISpace _space) {
        this(_space, new RenderMDParam());
    }
    
    public RenderMD(ISpace _space, RenderMDParam params) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(_space);

        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, 3.0, space) : new PotentialMasterMonatomic(this);
        
        parser = new ParseObj(params.file);

        int numAtoms = parser.nAtoms;
        
        double neighborRangeFac = 1.2;
        double sigma = Double.NaN;//not using neighbor lists because potential is different for each pair
        double lambda = params.lambda;
        if (params.useNeighborLists) {
            ((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*sigma*lambda);
        }

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setTemperature(params.temperature);
        integrator.setIsothermal(true);
        integrator.setTimeStep(params.timeStep);

        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        
        potentialBonded = new P2PenetrableSquareWell(space);
        potentialNonBonded = new P2Ideal(space);
        potential = new P2HardBondedList(this, potentialBonded, potentialNonBonded);
        
        ((P2PenetrableSquareWell)potentialBonded).setEpsilonCore(params.epsilonCore);
        ((P2PenetrableSquareWell)potentialBonded).setEpsilon(params.epsilon);
        ((P2PenetrableSquareWell)potentialBonded).setLambda(params.lambda);
        
        IAtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new IAtomType[]{leafType, leafType});

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
//        inflater.setTargetDensity(params.eta * 2 * space.D() / Math.PI);
        inflater.setTargetDensity(1.0);
        inflater.actionPerformed();
        
        potential.setBox(box);
        IAtomList leafList = box.getLeafList();
        for (int iLeaf=0; iLeaf<numAtoms; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            IVectorMutable pos = a.getPosition();
            pos.E((Vector3D)parser.vertices.get(iLeaf));
        }
        
        int nBonds = parser.bondList.size();
        for(int i=0; i<nBonds; i++) {
            BondInfo bond = parser.bondList.get(i);
            IAtom atom0 = leafList.getAtom(bond.i0);
            IAtom atom1 = leafList.getAtom(bond.i1);
            potential.setBonded(true, atom0, atom1, 0.95*bond.bondLengthSquared);
        }
        
        int bondSum = 0;
        int bondMax = 0;
        int bondMin = 0;
        for (int iLeaf=0; iLeaf<numAtoms; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            int bondCount = potential.getBondedList(a).size();
            bondSum += bondCount;
            bondMax = Math.max(bondMax, bondCount);
            bondMin = Math.min(bondMin, bondCount);
        }
        System.out.println("Avg, max, min bonds per atom: "+((float)bondSum/numAtoms)+" "+bondMax+" "+bondMin);

//        new ConfigurationFile(params.file).initializeCoordinates(box);
//        new ConfigurationCar(space).initializeCoordinates(box);
        
        integrator.setBox(box);

        if (params.useNeighborLists) { 
            NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager(box);
            integrator.getEventManager().addListener(nbrManager);
        }
        else {
            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        }
    }

//    public static void main(String[] args) {
//    	final String APP_NAME = "RenderMD";
//
//    	ISpace sp = Space3D.getInstance();
//        RenderMDParam params = new RenderMDParam();
//        final RenderMD sim = new RenderMD(sp, params);
//                
//        
//        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, sim.space, sim.getController());
//
//        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();
//
//        ArrayList<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
//
//        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
//        nSelector.setResetAction(new SimulationRestart(sim, sp, sim.getController()));
//        nSelector.setSpecies(sim.species);
//        nSelector.setBox(sim.box);
//        ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getAtomType(0),0.1);
//
//        nSelector.setPostAction(simGraphic.getPaintAction(sim.box));
//        simGraphic.add(nSelector);
//
//        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
//        
//        // Simulation Time
//        final DisplayTextBox displayCycles = new DisplayTextBox();
//
//        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integrator);
//        displayCycles.setPrecision(6);
//        DataPump pump= new DataPump(meterCycles,displayCycles);
//        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));
//        displayCycles.setLabel("Simulation time");
//
//        //temperature selector
//        DeviceThermoSlider tempSlider;
//        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
//        //tempSlider.setUnit(tUnit);
////        tempSlider.setPrecision(1);
//        tempSlider.setMinimum(0.0);
//        tempSlider.setMaximum(1500.0);
//        tempSlider.setSliderMajorValues(3);
//        tempSlider.setAdiabatic();
//        simGraphic.add(tempSlider);
//
//        JPanel statePanel = new JPanel(new GridBagLayout());
//        GridBagConstraints gbc2 = new GridBagConstraints();
//        gbc2.gridx = 0;  gbc2.gridy = 0;
//        statePanel.add(tempSlider.graphic(), gbc2);
//
//        DeviceBox sigBox, epsBox, lamBox, massBox;
//
////        JPanel parameterPanel = new JPanel(new GridLayout(0,1));
////        parameterPanel.add(sigBox.graphic());
////        parameterPanel.add(epsBox.graphic());
////        parameterPanel.add(lamBox.graphic());
////        parameterPanel.add(massBox.graphic());
//
//        JTabbedPane setupPanel = new JTabbedPane();
//        setupPanel.add(statePanel, "State");
//
//        
////        ModifierAtomDiameter sigModifier = new ModifierAtomDiameter();
////        sigModifier.setValue(sigma);
////        ModifierGeneral epsModifier = new ModifierGeneral(potentialSW, "epsilon");
////        ModifierGeneral lamModifier = new ModifierGeneral(potentialSW, "lambda");
////        ModifierGeneral massModifier = new ModifierGeneral(sim.species.getLeafType().getElement(),"mass");
////        sigBox.setModifier(sigModifier);
////        sigBox.setLabel("Core Diameter ("+Angstrom.UNIT.symbol()+")");
////        epsBox.setUnit(eUnit);
////        epsBox.setModifier(epsModifier);
////        lamBox.setModifier(lamModifier);
////        massBox.setModifier(massModifier);
////        massBox.setUnit(mUnit);
////        sigBox.setController(sim.getController());
////        epsBox.setController(sim.getController());
////        lamBox.setController(sim.getController());
////        massBox.setController(sim.getController());
//
//        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
//
//        //add meter and display for current kinetic temperature
//
//        MeterTemperature thermometer = new MeterTemperature(sim.box, sim.space.D());
//        DataFork temperatureFork = new DataFork();
//        final DataPump temperaturePump = new DataPump(thermometer,temperatureFork);
//        IntegratorListenerAction temperaturePumpListener = new IntegratorListenerAction(temperaturePump);
//        sim.integrator.getEventManager().addListener(temperaturePumpListener);
//        temperaturePumpListener.setInterval(1);
//        final AccumulatorAverageCollapsing temperatureAverage = new AccumulatorAverageCollapsing();
//        temperatureAverage.setPushInterval(20);
//        final AccumulatorHistory temperatureHistory = new AccumulatorHistory();
//        temperatureHistory.setTimeDataSource(timeCounter);
//        temperatureFork.setDataSinks(new IDataSink[]{temperatureAverage,temperatureHistory});
//        final DisplayTextBoxesCAE tBox = new DisplayTextBoxesCAE();
//        tBox.setAccumulator(temperatureAverage);
//        dataStreamPumps.add(temperaturePump);
//        //tBox.setUnit(tUnit);
//        tBox.setLabel("Measured Temperature (K)");
//        tBox.setLabelPosition(CompassDirection.NORTH);
//
//        // Number density box
//        MeterDensity densityMeter = new MeterDensity(sim.getSpace());
//        densityMeter.setBox(sim.box);
//        final DisplayTextBox densityBox = new DisplayTextBox();
//        //densityBox.setUnit(dUnit);
//        final DataPump densityPump = new DataPump(densityMeter, densityBox);
//        IntegratorListenerAction densityPumpListener = new IntegratorListenerAction(densityPump);
//        sim.integrator.getEventManager().addListener(densityPumpListener);
//        densityPumpListener.setInterval(1);
//        dataStreamPumps.add(densityPump);
//        densityBox.setLabel("Density");
//        
//        MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotentialMaster(), sim.box);
//        final AccumulatorHistory energyHistory = new AccumulatorHistory();
//        energyHistory.setTimeDataSource(timeCounter);
//        final DataSinkExcludeOverlap eExcludeOverlap = new DataSinkExcludeOverlap();
//        eExcludeOverlap.setDataSink(energyHistory);
//        final DataPump energyPump = new DataPump(eMeter, eExcludeOverlap);
//        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
//        sim.integrator.getEventManager().addListener(energyPumpListener);
//        dataStreamPumps.add(energyPump);
//        
//        MeterPotentialEnergyFromIntegrator peMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//        final AccumulatorHistory peHistory = new AccumulatorHistory();
//        peHistory.setTimeDataSource(timeCounter);
//        final AccumulatorAverageCollapsing peAccumulator = new AccumulatorAverageCollapsing();
//        peAccumulator.setPushInterval(2);
//        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peAccumulator});
//        final DataSinkExcludeOverlap peExcludeOverlap = new DataSinkExcludeOverlap();
//        peExcludeOverlap.setDataSink(peFork);
//        final DataPump pePump = new DataPump(peMeter, peExcludeOverlap);
//        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
//        sim.integrator.getEventManager().addListener(pePumpListener);
//        dataStreamPumps.add(pePump);
//
//        MeterKineticEnergyFromIntegrator keMeter = new MeterKineticEnergyFromIntegrator(sim.integrator);
//        final AccumulatorHistory keHistory = new AccumulatorHistory();
//        keHistory.setTimeDataSource(timeCounter);
//        final DataPump kePump = new DataPump(keMeter, keHistory);
//        IntegratorListenerAction kePumpListener = new IntegratorListenerAction(kePump);
//        sim.integrator.getEventManager().addListener(kePumpListener);
//        dataStreamPumps.add(kePump);
//        int numAtoms = sim.box.getLeafList().getAtomCount();
//        energyPumpListener.setInterval(numAtoms > 120 ? 1 : 120/numAtoms);
//        kePumpListener.setInterval(numAtoms > 120 ? 1 : 120/numAtoms);
//        pePumpListener.setInterval(numAtoms > 120 ? 1 : 120/numAtoms);
//        
//        final DisplayPlot ePlot = new DisplayPlot();
//        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        ePlot.setLegend(new DataTag[]{energyHistory.getTag()}, "Total");
//        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        ePlot.setLegend(new DataTag[]{peHistory.getTag()}, "Potential");
//        keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
//        ePlot.setLegend(new DataTag[]{keHistory.getTag()}, "Kinetic");
//
//        ePlot.getPlot().setTitle("Energy History (J/mol)");
//        ePlot.setDoLegend(true);
//        ePlot.setLabel("Energy");
//        //ePlot.setUnit(eUnit);
//        ePlot.setXUnit(Picosecond.UNIT);
//        
//        MeterPressureHard pMeter = new MeterPressureHard(sim.getSpace());
//        pMeter.setIntegrator(sim.integrator);
//        final AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
//        final DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator);
//        sim.integrator.getEventManager().addListener(pPump);
//        pAccumulator.setPushInterval(50);
//        dataStreamPumps.add(pPump);
//
//        final DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
//        peDisplay.setAccumulator(peAccumulator);
//        peDisplay.setLabel("Potential Energy (J/mol)");
////        peDisplay.setUnit(eUnit);
//        
//        //DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
//        
//        simGraphic.getPanel().controlPanel.add(setupPanel, vertGBC);
//        //simGraphic.getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);
//        //simGraphic.add(configButton);
//        //simGraphic.add(velocityButton);
//
//        simGraphic.add(displayCycles);
//        simGraphic.add(densityBox);
//        simGraphic.add(tBox);
//        simGraphic.add(peDisplay);
//        
//        java.awt.Dimension d = ePlot.getPlot().getPreferredSize();
//        d.width -= 50;
//        ePlot.getPlot().setSize(d);
//
//
//        simGraphic.makeAndDisplayFrame(APP_NAME);
//        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
//        colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);
//    }

    public static RenderMDParam getParameters() {
        return new RenderMDParam();
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class RenderMDParam extends ParameterBase {
        public double eta = 0.7;
        public boolean useNeighborLists = false;
        public String file = "/Users/kofke/Documents/workspace/car self-assembly/mustang.txt";
        public double epsilonCore = 10;
        public double lambda = 1.1;
        public double epsilon = 10;
        public double temperature = 1.0;
        public double timeStep = 0.02;
        
    }
}
