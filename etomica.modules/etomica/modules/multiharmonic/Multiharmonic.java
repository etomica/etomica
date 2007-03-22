package etomica.modules.multiharmonic;

import etomica.action.SimulationDataAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.phase.Phase;
import etomica.potential.P1Harmonic;
import etomica.simulation.Simulation;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsing;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class Multiharmonic extends Simulation {

    public Multiharmonic() {
        super(Space1D.getInstance());
        defaults.makeLJDefaults();
        defaults.atomSize = 0.02;
        double x0 = 1;
        species = new SpeciesSpheresMono(this);
        getSpeciesRoot().addSpecies(species);
        phase = new Phase(this);
        phase.getBoundary().setDimensions(new Vector1D(3.0));
        controller = getController();
        integrator = new IntegratorVelocityVerlet(this);
        integrator.setPhase(phase);
        integrator.setTimeStep(0.02);
        integrator.setIsothermal(true);
        integrator.setTemperature(1.0);
        species.getAgent(phase).setNMolecules(20);
        potentialA = new P1Harmonic(space);
        potentialA.setX0(new Vector1D(x0));
        potentialA.setSpringConstant(1.0);
        potentialMaster.addPotential(potentialA, new Species[] {species});
        
        phase.getAgent(species).setNMolecules(20);
        
        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)iterator.nextAtom();
            a.getCoord().getPosition().setX(0,x0);
        }
        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);

        potentialB = new P1Harmonic(space);
        potentialB.setX0(new Vector1D(x0+1));
        potentialB.setSpringConstant(10);
        meter = new MeterFreeEnergy(potentialA, potentialB);
        meter.setPhase(phase);
        accumulator = new AccumulatorAverage(this);
        accumulator.setBlockSize(100);
        dataPump = new DataPump(meter, accumulator);
        new IntervalActionAdapter(dataPump, integrator);
        
        meterEnergy = new MeterEnergy(potentialMaster);
        meterEnergy.setPhase(phase);
        accumulatorEnergy = new AccumulatorAverage(this);
        accumulatorEnergy.setBlockSize(100);
        DataPump dataPumpEnergy = new DataPump(meterEnergy, accumulatorEnergy);
        new IntervalActionAdapter(dataPumpEnergy, integrator);
        
        register(meter,dataPump);
        register(meterEnergy,dataPumpEnergy);
        
        historyEnergy = new AccumulatorHistory(new HistoryCollapsing());
        accumulatorEnergy.addDataSink(historyEnergy, new AccumulatorAverage.StatType[] {AccumulatorAverage.StatType.AVERAGE});

        timeCounter = new DataSourceCountTime();
        integrator.addListener(timeCounter);
        
        historyEnergy.setTimeDataSource(timeCounter);
    }

    private static final long serialVersionUID = 1L;
    MeterEnergy meterEnergy;
    AccumulatorAverage accumulatorEnergy;
    AccumulatorHistory historyEnergy;
    SpeciesSpheresMono species;
    Phase phase;
    Controller controller;
    P1Harmonic potentialA, potentialB;
    IntegratorVelocityVerlet integrator;
    ActivityIntegrate activityIntegrate;
    MeterFreeEnergy meter;
    AccumulatorAverage accumulator;
    DataPump dataPump;
    SimulationDataAction resetAccumulators;
    DataSourceCountTime timeCounter;
}
