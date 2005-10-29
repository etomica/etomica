package etomica.multiharmonic;

import etomica.action.SimulationDataAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceFunction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.phase.Phase;
import etomica.potential.P1Harmonic;
import etomica.simulation.Simulation;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Function;
import etomica.util.HistoryCollapsingAverage;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Oct 28, 2005 by kofke
 */
public class Multiharmonic extends Simulation {

    /**
     * 
     */
    public Multiharmonic() {
        super(Space1D.getInstance());
        defaults.atomSize = 0.2;
        double x0 = 10;
        species = new SpeciesSpheresMono(this);
        phase = new Phase(this);
        controller = getController();
        integrator = new IntegratorVelocityVerlet(this);
        integrator.addPhase(phase);
        integrator.setTimeStep(0.02);
        integrator.setIsothermal(true);
        integrator.setTemperature(1.0);
        species.getAgent(phase).setNMolecules(20);
        System.out.println(integrator.getTimeStep());
        potentialA = new P1Harmonic(space);
        potentialA.setX0(new Vector1D(x0));
        potentialA.setSpringConstant(1.0);
        potentialMaster.setSpecies(potentialA, new Species[] {species});
        
        phase.makeMolecules();
        
        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            a.coord.position().setX(0,x0);
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
        dataPump = new DataPump(meter, accumulator);
        new IntervalActionAdapter(dataPump, integrator);
        
        register(meter,dataPump);
        
        history = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
        accumulator.addDataSink(history, new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE});
         
    }


    SpeciesSpheresMono species;
    Phase phase;
    Controller controller;
    P1Harmonic potentialA, potentialB;
    IntegratorVelocityVerlet integrator;
    ActivityIntegrate activityIntegrate;
    MeterFreeEnergy meter;
    AccumulatorAverage accumulator;
    DataPump dataPump;
    AccumulatorHistory history;
    SimulationDataAction resetAccumulators;
}
