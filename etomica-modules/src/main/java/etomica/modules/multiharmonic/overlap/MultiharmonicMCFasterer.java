/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.overlap;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
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
    protected final PotentialComputeField potentialMasterA, potentialMasterB;
    protected final P1Harmonic potentialA, potentialB;
    protected final IntegratorMCFasterer integratorA, integratorB;
    protected final IntegratorOverlap integratorOS;

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

        potentialMasterA = new PotentialComputeField(getSpeciesManager(), boxA);
        potentialMasterB = new PotentialComputeField(getSpeciesManager(), boxB);

        integratorA = new IntegratorMCFasterer(potentialMasterA, random, 1.0, boxA);
        potentialA = new P1Harmonic(space);
        integratorA.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integratorA, potentialA, random));
        potentialMasterA.setFieldPotential(species.getLeafType(), potentialA);

        integratorB = new IntegratorMCFasterer(potentialMasterB, this.getRandom(), 1.0, boxB);
        integratorB.setTemperature(1.0);
        potentialB = new P1Harmonic(space);
        integratorB.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integratorB, potentialB, random));
        potentialMasterB.setFieldPotential(species.getLeafType(), potentialB);

        integratorOS = new IntegratorOverlap(new Integrator[]{integratorA, integratorB});

        integratorOS.setAdjustStepFraction(false);
        integratorOS.setRefStepFraction(0.5);

        getController().addActivity(new ActivityIntegrate(integratorOS, true));
    }
}
