/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.integrator.IntegratorMCFasterer;
import etomica.potential.P1Harmonic;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesGeneral;


/**
 * MC version of multi-harmonic simulation.  This version runs much faster.
 *
 * @author Andrew Schultz
 */
public class MultiharmonicMCFasterer extends Simulation {

    MeterPotentialEnergyFromIntegratorFasterer meterEnergy;
    AccumulatorAverageCollapsing accumulatorEnergy;
    AccumulatorHistory historyEnergy;
    SpeciesGeneral species;
    Box box;
    P1Harmonic potentialA, potentialB;
    IntegratorMCFasterer integrator;
    MeterFreeEnergy meter;
    AccumulatorAverageCollapsing accumulator;
    DataPumpListener dataPump, dataPumpEnergy;
    DataSourceCountSteps stepCounter;

    public MultiharmonicMCFasterer() {
        super(Space1D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.getBoundary().setBoxSize(new Vector1D(6.0));
        PotentialComputeField potentialMaster = new PotentialComputeField(getSpeciesManager(), box);
        integrator = new IntegratorMCFasterer(this, potentialMaster, box);
        integrator.setTemperature(1.0);
        potentialA = new P1Harmonic(space);
        integrator.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integrator, potentialA, random));
        potentialMaster.setFieldPotential(species.getLeafType(), potentialA);

        box.setNMolecules(species, 10);

        getController().setSleepPeriod(1);
        getController().addActivity(new ActivityIntegrate(integrator));

        potentialB = new P1Harmonic(space);
        meter = new MeterFreeEnergy(potentialA, potentialB);
        meter.setBox(box);
        accumulator = new AccumulatorAverageCollapsing();
        dataPump = new DataPumpListener(meter, accumulator);
        integrator.getEventManager().addListener(dataPump);

        meterEnergy = new MeterPotentialEnergyFromIntegratorFasterer(integrator);
        accumulatorEnergy = new AccumulatorAverageCollapsing();
        dataPumpEnergy = new DataPumpListener(meterEnergy, accumulatorEnergy);
        integrator.getEventManager().addListener(dataPumpEnergy);

        historyEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard(102, 3));
        accumulatorEnergy.addDataSink(historyEnergy, new AccumulatorAverage.StatType[]{accumulatorEnergy.AVERAGE});

        stepCounter = new DataSourceCountSteps(integrator);

        historyEnergy.setTimeDataSource(stepCounter);
    }
}
