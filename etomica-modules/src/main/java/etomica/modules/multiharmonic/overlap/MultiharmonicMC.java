/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.overlap;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.modules.multiharmonic.MCMoveMultiHarmonic;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P1Harmonic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
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
public class MultiharmonicMC extends Simulation {

    protected final SpeciesGeneral species;
    protected final Box boxA, boxB;
    protected final PotentialMaster potentialMasterA, potentialMasterB;
    protected final P1Harmonic potentialA, potentialB;
    protected final IntegratorMC integratorA, integratorB;
    protected final IntegratorOverlap integratorOS;
    
    public MultiharmonicMC() {
        super(Space1D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        potentialMasterA = new PotentialMasterMonatomic(getSpeciesManager());
        potentialMasterB = new PotentialMasterMonatomic(getSpeciesManager());

        boxA = this.makeBox(new BoundaryRectangularNonperiodic(space));
        boxA.getBoundary().setBoxSize(new Vector1D(6.0));
        boxA.setNMolecules(species, 10);
        boxB = this.makeBox(new BoundaryRectangularNonperiodic(space));
        boxB.getBoundary().setBoxSize(new Vector1D(6.0));
        boxB.setNMolecules(species, 10);

        integratorA = new IntegratorMC(this, potentialMasterA, boxA);
        integratorA.setTemperature(1.0);
        potentialA = new P1Harmonic(space);
        integratorA.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integratorA, potentialA, random));
        potentialMasterA.addPotential(potentialA, new AtomType[]{species.getLeafType()});

        integratorB = new IntegratorMC(this, potentialMasterB, boxB);
        integratorB.setTemperature(1.0);
        potentialB = new P1Harmonic(space);
        integratorB.getMoveManager().addMCMove(new MCMoveMultiHarmonic(integratorB, potentialB, random));
        potentialMasterB.addPotential(potentialB, new AtomType[]{species.getLeafType()});

        integratorOS = new IntegratorOverlap(new Integrator[]{integratorA, integratorB});

        integratorOS.setAdjustStepFraction(false);
        integratorOS.setRefStepFraction(0.5);

        getController().addActivity(new ActivityIntegrate(integratorOS, true));
    }
}
