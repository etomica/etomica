package etomica.simulation.prototypes;
import etomica.action.Action;
import etomica.action.BoxImposePbc;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.box.Box;
import etomica.potential.P1HardBoundary;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D_noNbr extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public ActivityIntegrate activityIntegrate;
    public AccumulatorAverage pressureAverage;
    public AccumulatorHistory pressureHistory;
    public AccumulatorAverage temperatureAverage;
    public AccumulatorHistory temperatureHistory;
    public Box box;
    public SpeciesSpheresMono species;
    public IntegratorHard integrator;

    public HSMD2D_noNbr() {
    	this(Space2D.getInstance());
    }
    
    public HSMD2D_noNbr(Space2D space) {
        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
	    box = new Box(this);
        addBox(box);
        box.setNMolecules(species, 64);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(box);
	    P2HardSphere potential = new P2HardSphere(space);
	    potentialMaster.addPotential(potential,new Species[]{species,species});
        P1HardBoundary potentialBoundary = new P1HardBoundary(space);
        potentialMaster.addPotential(potentialBoundary, new Species[] {species});
//        potentialBoundary.setActive(0,true,true);
//        potentialBoundary.setActive(1,true,true);
//        potentialBoundary.setActive(0,false,true);
//        potentialBoundary.setActive(1,false,true);
        
        integrator.addIntervalAction(new BoxImposePbc(box));
        integrator.setBox(box);
        integrator.setNullPotential(new P1HardPeriodic(space));
//        integrator.setIsothermal(true);
        
//        MeterPressureHard meterPressure = new MeterPressureHard(integrator);
//        meterPressure.setBox(box);
//        pressureAverage = new AccumulatorAverage();
//        DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
//        IntervalActionAdapter pressureAction = new IntervalActionAdapter(pressurePump, integrator);
//
//        pressureHistory = new AccumulatorHistory();
//        pressureAverage.makeDataPusher(
//          new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE}).
//                                      addDataSink(pressureHistory);
        
        MeterTemperature meterTemperature = new MeterTemperature();
        meterTemperature.setBox(box);
        temperatureAverage = new AccumulatorAverage();
        DataPump temperaturePump = new DataPump(meterTemperature, temperatureAverage);
        integrator.addIntervalAction(temperaturePump);

//        pressureHistory = new AccumulatorHistory();
//        pressureAverage.makeDataPusher(
//          new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE}).
//                                      addDataSink(pressureHistory);
        temperatureHistory = new AccumulatorHistory();
        temperatureAverage.addDataSink(temperatureHistory,new AccumulatorAverage.StatType[] {StatType.AVERAGE});
        DataSourceCountTime timeCounter = new DataSourceCountTime(integrator);
        temperatureHistory.setTimeDataSource(timeCounter);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
    	final String APP_NAME = "HSMD2D no Nbr";

        final HSMD2D_noNbr sim = new HSMD2D_noNbr();
        final SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        sim.activityIntegrate.setDoSleep(true);
//        DisplayTextBoxesCAE pressureDisplay = new DisplayTextBoxesCAE();
//        pressureDisplay.setAccumulator(sim.pressureAverage);
//        DisplayPlot pressurePlot = new DisplayPlot();
//        sim.pressureHistory.addDataSink(pressurePlot.makeDataSink());
        DisplayTextBoxesCAE temperatureDisplay = new DisplayTextBoxesCAE();
        temperatureDisplay.setAccumulator(sim.temperatureAverage);
        DisplayPlot temperaturePlot = new DisplayPlot();
        temperaturePlot.setLabel("Temp");
        sim.temperatureHistory.setDataSink(temperaturePlot.getDataSet().makeDataSink());
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);
        Action repaintAction = graphic.getDisplayBoxPaintAction(sim.box);

        nSelector.setPostAction(repaintAction);
        graphic.getController().getReinitButton().setPostAction(repaintAction);

        sim.integrator.addIntervalAction(repaintAction);

        DeviceThermoSelector thermo = new DeviceThermoSelector(sim, sim.integrator);
        graphic.add(nSelector);
        graphic.add(thermo);
//        graphic.add(pressureDisplay);
//        graphic.add(pressurePlot);
        graphic.add(temperatureDisplay);
        graphic.add(temperaturePlot);
		graphic.makeAndDisplayFrame(APP_NAME);
    }//end of main
    
}