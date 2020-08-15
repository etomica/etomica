/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.umbrella;


import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.modules.multiharmonic.MCMoveMultiHarmonic;
import etomica.potential.P1Harmonic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesSpheresMono;


/**
 * MC version of multi-harmonic simulation.  This version runs much faster.
 *
 * @author Andrew Schultz
 */
public class MultiharmonicMC extends Simulation {

    private static final long serialVersionUID = 1L;
    protected final SpeciesSpheresMono species;
    protected final Box box;
    protected final P1Harmonic potentialA, potentialB;
    protected final IntegratorMC integrator;
    protected final MCMoveMultiHarmonic moveA, moveB;
    
    protected final MeterUmbrella meterUmbrella;
    protected final AccumulatorRatioAverageCovariance accumulator;
    protected final DataPump dataPumpA;
    public MultiharmonicMC() {
        super(Space1D.getInstance());
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        PotentialMaster potentialMasterA = new PotentialMasterMonatomic(this);
        PotentialMaster potentialMasterB = new PotentialMasterMonatomic(this);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.getBoundary().setBoxSize(new Vector1D(6.0));
        box.setNMolecules(species, 10);

        integrator = new IntegratorMC(this, potentialMasterA, box);
        integrator.setTemperature(1.0);
        potentialA = new P1Harmonic(space);
        moveA = new MCMoveMultiHarmonic(potentialA, random);
        integrator.getMoveManager().addMCMove(moveA);
        potentialMasterA.addPotential(potentialA, new AtomType[]{species.getLeafType()});

        potentialB = new P1Harmonic(space);
        moveB = new MCMoveMultiHarmonic(potentialB, random);
        integrator.getMoveManager().addMCMove(moveB);
        potentialMasterB.addPotential(potentialB, new AtomType[]{species.getLeafType()});

        MeterPotentialEnergy meterPEAinA = new MeterPotentialEnergy(potentialMasterA, box);
        MeterPotentialEnergy meterPEBinA = new MeterPotentialEnergy(potentialMasterB, box);
        meterUmbrella = new MeterUmbrella(meterPEAinA, meterPEBinA, 1.0);

        accumulator = new AccumulatorRatioAverageCovariance(1);
        dataPumpA = new DataPump(meterUmbrella, accumulator);
        integrator.getEventManager().addListener(new IntegratorListenerAction(dataPumpA));


        getController().addActivity(new ActivityIntegrate2(integrator, 0, true));
    }
}
