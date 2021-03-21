/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.mbar;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergyFasterer;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMCFasterer;
import etomica.modules.multiharmonic.MCMoveMultiHarmonic;
import etomica.overlap.IntegratorOverlap;
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

    protected final SpeciesGeneral species;
    protected final Box boxA, boxB;
    protected final P1Harmonic potentialA, potentialB;
    protected final IntegratorMCFasterer integratorA, integratorB;
    protected final IntegratorOverlap integratorOS;

    protected final MeterMBAR meterOverlapA, meterOverlapB;

    public MultiharmonicMCFasterer() {
        super(Space1D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        boxA = this.makeBox(new BoundaryRectangularNonperiodic(space));
        boxA.getBoundary().setBoxSize(new Vector1D(6.0));
        boxA.setNMolecules(species, 10);
        boxB = this.makeBox(new BoundaryRectangularNonperiodic(space));
        boxB.getBoundary().setBoxSize(new Vector1D(6.0));
        boxB.setNMolecules(species, 10);

        PotentialComputeField potentialMasterA = new PotentialComputeField(getSpeciesManager(), boxA);
        PotentialComputeField potentialMasterAinB = new PotentialComputeField(getSpeciesManager(), boxB);
        PotentialComputeField potentialMasterB = new PotentialComputeField(getSpeciesManager(), boxB);
        PotentialComputeField potentialMasterBinA = new PotentialComputeField(getSpeciesManager(), boxA);
        potentialMasterAinB.init();
        potentialMasterBinA.init();

        integratorA = new IntegratorMCFasterer(potentialMasterA, random, 1.0, boxA);
        integratorA.setTemperature(1.0);
        potentialA = new P1Harmonic(space);
        integratorA.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integratorA, potentialA, random));
        potentialMasterA.setFieldPotential(species.getLeafType(), potentialA);
        potentialMasterAinB.setFieldPotential(species.getLeafType(), potentialA);

        integratorB = new IntegratorMCFasterer(potentialMasterB, random, 1.0, boxB);
        integratorB.setTemperature(1.0);
        potentialB = new P1Harmonic(space);
        integratorB.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integratorB, potentialB, random));
        potentialMasterB.setFieldPotential(species.getLeafType(), potentialB);
        potentialMasterBinA.setFieldPotential(species.getLeafType(), potentialB);


        MeterPotentialEnergyFromIntegratorFasterer meterPEAinA = new MeterPotentialEnergyFromIntegratorFasterer(integratorA);
        MeterPotentialEnergyFasterer meterPEAinB = new MeterPotentialEnergyFasterer(potentialMasterAinB);
        MeterPotentialEnergyFasterer meterPEBinA = new MeterPotentialEnergyFasterer(potentialMasterBinA);
        MeterPotentialEnergyFromIntegratorFasterer meterPEBinB = new MeterPotentialEnergyFromIntegratorFasterer(integratorB);
        meterOverlapA = new MeterMBAR(new DataSourceScalar[]{meterPEAinA, meterPEBinA}, 1.0);
        meterOverlapA.setNumAlpha(15);
        meterOverlapA.setAlpha(new double[]{1}, new double[]{5});
        meterOverlapB = new MeterMBAR(new DataSourceScalar[]{meterPEAinB, meterPEBinB}, 1.0);
        meterOverlapB.setNumAlpha(15);
        meterOverlapB.setAlpha(new double[]{1}, new double[]{5});
        meterOverlapB.setBoxIndex(1);

        integratorA.getEventManager().addListener(meterOverlapA);
        integratorB.getEventManager().addListener(meterOverlapB);

        integratorOS = new IntegratorOverlap(new Integrator[]{integratorA, integratorB});
        integratorOS.setRefStepFraction(0.5);
        integratorOS.setAdjustStepFraction(false);

        getController().addActivity(new ActivityIntegrate(integratorOS, true));
    }
}
