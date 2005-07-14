package etomica.simulations;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.potential.P1HardBoundary;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.space2d.Space2D;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D_noNbr extends Simulation {
    
    public ActivityIntegrate activityIntegrate;
    public AccumulatorAverage pressureAverage;
    public AccumulatorHistory pressureHistory;
    public AccumulatorAverage temperatureAverage;
    public AccumulatorHistory temperatureHistory;
    public Phase phase;
    public SpeciesSpheresMono species;
    public IntegratorHard integrator;

    public HSMD2D_noNbr() {
    	this(Space2D.getInstance());
    }
    
    public HSMD2D_noNbr(Space2D space) {
        super(space);

        integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(64);
	    phase = new Phase(this);
	    P2HardSphere potential = new P2HardSphere(space);
	    potentialMaster.setSpecies(potential,new Species[]{species,species});
        P1HardBoundary potentialBoundary = new P1HardBoundary(space);
        potentialMaster.setSpecies(potentialBoundary, new Species[] {species});
//        potentialBoundary.setActive(0,true,true);
//        potentialBoundary.setActive(1,true,true);
//        potentialBoundary.setActive(0,false,true);
//        potentialBoundary.setActive(1,false,true);
        
        integrator.addListener(new PhaseImposePbc(phase));
        integrator.setNullPotential(new P1HardPeriodic(space));
        integrator.addPhase(phase);
//        integrator.setIsothermal(true);
        
//        MeterPressureHard meterPressure = new MeterPressureHard(integrator);
//        meterPressure.setPhase(phase);
//        pressureAverage = new AccumulatorAverage();
//        DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
//        IntervalActionAdapter pressureAction = new IntervalActionAdapter(pressurePump, integrator);
//
//        pressureHistory = new AccumulatorHistory();
//        pressureAverage.makeDataPusher(
//          new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE}).
//                                      addDataSink(pressureHistory);
        
        MeterTemperature meterTemperature = new MeterTemperature();
        meterTemperature.setPhase(phase);
        temperatureAverage = new AccumulatorAverage(meterTemperature.getDataInfo());
        DataPump temperaturePump = new DataPump(meterTemperature, temperatureAverage);
        new IntervalActionAdapter(temperaturePump, integrator);

//        pressureHistory = new AccumulatorHistory();
//        pressureAverage.makeDataPusher(
//          new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE}).
//                                      addDataSink(pressureHistory);
        temperatureHistory = new AccumulatorHistory(temperatureAverage.getDataInfo());
        temperatureAverage.addDataSink(temperatureHistory,new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE});
}
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D_noNbr sim = new HSMD2D_noNbr();
        final SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.activityIntegrate.setDoSleep(true);
//        DisplayBoxesCAE pressureDisplay = new DisplayBoxesCAE();
//        pressureDisplay.setAccumulator(sim.pressureAverage);
//        DisplayPlot pressurePlot = new DisplayPlot();
//        sim.pressureHistory.addDataSink(pressurePlot.makeDataSink());
        DisplayBoxesCAE temperatureDisplay = new DisplayBoxesCAE();
        temperatureDisplay.setAccumulator(sim.temperatureAverage);
        DisplayPlot temperaturePlot = new DisplayPlot();
        sim.temperatureHistory.addDataSink(temperaturePlot.getDataTable().makeColumn(sim.temperatureHistory.getDataInfo()));
        DeviceNSelector nSelector = new DeviceNSelector(sim, sim.phase.getAgent(sim.species));
        DeviceThermoSelector thermo = new DeviceThermoSelector(sim.getController(), sim.integrator);
        graphic.add(nSelector);
        graphic.add(thermo);
//        graphic.add(pressureDisplay);
//        graphic.add(pressurePlot);
        graphic.add(temperatureDisplay);
        graphic.add(temperaturePlot);
		graphic.makeAndDisplayFrame();
    }//end of main
    
}