/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.action.SimulationDataAction;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.P1Harmonic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesSpheresMono;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class Multiharmonic extends Simulation {

    private static final long serialVersionUID = 1L;
    MeterEnergy meterEnergy;
    AccumulatorAverageCollapsing accumulatorEnergy;
    AccumulatorHistory historyEnergy;
    SpeciesSpheresMono species;
    Box box;
    P1Harmonic potentialA, potentialB;
    IntegratorVelocityVerlet integrator;
    MeterFreeEnergy meter;
    AccumulatorAverageCollapsing accumulator;
    DataPump dataPump, dataPumpEnergy;
    SimulationDataAction resetAccumulators;
    DataSourceCountTime timeCounter;
    public Multiharmonic() {
        super(Space1D.getInstance());
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        double x0 = 0;
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.getBoundary().setBoxSize(new Vector1D(6.0));
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.02);
        integrator.setIsothermal(true);
        integrator.setTemperature(1.0);
        potentialA = new P1Harmonic(space);
        potentialA.setX0(new Vector1D(x0));
        potentialA.setSpringConstant(1.0);
        potentialMaster.addPotential(potentialA, new AtomType[]{species.getLeafType()});

        box.setNMolecules(species, 20);

        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
        iterator.setBox(box);
        iterator.reset();
        for (IAtom a = iterator.nextAtom(); a != null;
             a = iterator.nextAtom()) {
            a.getPosition().setX(0, x0);
        }
        getController().setSleepPeriod(1);
        getController().addActivity(new ActivityIntegrate(integrator));

        potentialB = new P1Harmonic(space);
        potentialB.setX0(new Vector1D(x0 + 1));
        potentialB.setSpringConstant(10);
        meter = new MeterFreeEnergy(potentialA, potentialB);
        meter.setBox(box);
        accumulator = new AccumulatorAverageCollapsing();
        dataPump = new DataPump(meter, accumulator);
        integrator.getEventManager().addListener(new IntegratorListenerAction(dataPump));

        meterEnergy = new MeterEnergy(potentialMaster, box);
        accumulatorEnergy = new AccumulatorAverageCollapsing();
        dataPumpEnergy = new DataPump(meterEnergy, accumulatorEnergy);
        integrator.getEventManager().addListener(new IntegratorListenerAction(dataPumpEnergy));

        historyEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        accumulatorEnergy.addDataSink(historyEnergy, new AccumulatorAverage.StatType[]{accumulatorEnergy.AVERAGE});

        timeCounter = new DataSourceCountTime(integrator);

        historyEnergy.setTimeDataSource(timeCounter);
    }
}
