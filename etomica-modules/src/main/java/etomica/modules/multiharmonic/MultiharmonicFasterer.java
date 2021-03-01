/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterEnergyFromIntegratorFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.potential.P1Harmonic;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesGeneral;


public class MultiharmonicFasterer extends Simulation {

    MeterEnergyFromIntegratorFasterer meterEnergy;
    AccumulatorAverageCollapsing accumulatorEnergy;
    AccumulatorHistory historyEnergy;
    SpeciesGeneral species;
    Box box;
    P1Harmonic potentialA, potentialB;
    IntegratorVelocityVerletFasterer integrator;
    MeterFreeEnergy meter;
    AccumulatorAverageCollapsing accumulator;
    DataPumpListener dataPump, dataPumpEnergy;
    DataSourceCountTimeFasterer timeCounter;

    public MultiharmonicFasterer() {
        super(Space1D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        double x0 = 0;
        box.getBoundary().setBoxSize(new Vector1D(6.0));
        PotentialComputeField potentialMaster = new PotentialComputeField(getSpeciesManager(), box);
        integrator = new IntegratorVelocityVerletFasterer(potentialMaster, random, 0.02, 1.0, box);
        integrator.setIsothermal(true);
        potentialA = new P1Harmonic(space);
        potentialA.setX0(new Vector1D(x0));
        potentialA.setSpringConstant(1.0);
        potentialMaster.setFieldPotential(species.getLeafType(), potentialA);

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
        dataPump = new DataPumpListener(meter, accumulator);
        integrator.getEventManager().addListener(dataPump);

        meterEnergy = new MeterEnergyFromIntegratorFasterer(integrator);
        accumulatorEnergy = new AccumulatorAverageCollapsing();
        dataPumpEnergy = new DataPumpListener(meterEnergy, accumulatorEnergy);
        integrator.getEventManager().addListener(dataPumpEnergy);

        historyEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        accumulatorEnergy.addDataSink(historyEnergy, new AccumulatorAverage.StatType[]{accumulatorEnergy.AVERAGE});

        timeCounter = new DataSourceCountTimeFasterer(integrator);

        historyEnergy.setTimeDataSource(timeCounter);
    }
}
