/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.action.SimulationDataAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.IController;
import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P1Harmonic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesSpheresMono;
import etomica.data.history.HistoryCollapsingDiscard;


/**
 * MC version of multi-harmonic simulation.  This version runs much faster.
 *
 * @author Andrew Schultz
 */
public class MultiharmonicMC extends Simulation {

    public MultiharmonicMC() {
        super(Space1D.getInstance());
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        box = new Box(new BoundaryRectangularNonperiodic(space), space);
        addBox(box);
        box.getBoundary().setBoxSize(new Vector1D(3.0));
        controller = getController();
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setBox(box);
        integrator.setTemperature(1.0);
        potentialA = new P1Harmonic(space);
        integrator.getMoveManager().addMCMove(new MCMoveMultiHarmonic(potentialA, random));
        potentialMaster.addPotential(potentialA, new IAtomType[] {species.getLeafType()});
        
        box.setNMolecules(species, 10);
        
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);

        potentialB = new P1Harmonic(space);
        meter = new MeterFreeEnergy(potentialA, potentialB);
        meter.setBox(box);
        accumulator = new AccumulatorAverageCollapsing();
        dataPump = new DataPump(meter, accumulator);
        integrator.getEventManager().addListener(new IntegratorListenerAction(dataPump));
        
        meterEnergy = new MeterPotentialEnergy(potentialMaster);
        meterEnergy.setBox(box);
        accumulatorEnergy = new AccumulatorAverageCollapsing();
        dataPumpEnergy = new DataPump(meterEnergy, accumulatorEnergy);
        integrator.getEventManager().addListener(new IntegratorListenerAction(dataPumpEnergy));
        
        historyEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        accumulatorEnergy.addDataSink(historyEnergy, new AccumulatorAverage.StatType[] {accumulatorEnergy.AVERAGE});

        stepCounter = new DataSourceCountSteps(integrator);
        
        historyEnergy.setTimeDataSource(stepCounter);
    }

    private static final long serialVersionUID = 1L;
    MeterPotentialEnergy meterEnergy;
    AccumulatorAverageCollapsing accumulatorEnergy;
    AccumulatorHistory historyEnergy;
    SpeciesSpheresMono species;
    Box box;
    IController controller;
    P1Harmonic potentialA, potentialB;
    IntegratorMC integrator;
    ActivityIntegrate activityIntegrate;
    MeterFreeEnergy meter;
    AccumulatorAverageCollapsing accumulator;
    DataPump dataPump, dataPumpEnergy;
    SimulationDataAction resetAccumulators;
    DataSourceCountSteps stepCounter;
}
