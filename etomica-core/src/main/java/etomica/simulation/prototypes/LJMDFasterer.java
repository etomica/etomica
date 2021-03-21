/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMDFasterer extends Simulation {

    public IntegratorVelocityVerletFasterer integrator;
    public SpeciesGeneral species;
    public Box box;
    public MeterPotentialEnergyFromIntegratorFasterer energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPumpListener pump;


    public LJMDFasterer(Space space, double density, int numAtoms, boolean useNbrLists) {
        super(space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox();
        double rc = 3;
        PotentialMasterFasterer potentialMaster = useNbrLists ? new PotentialMasterListFasterer(getSpeciesManager(), box, 2, rc + 1, BondingInfo.noBonding())
                : new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());
        double sigma = 1.0;
        integrator = new IntegratorVelocityVerletFasterer(potentialMaster, this.getRandom(), 0.05, 1.0, box);
        integrator.setTimeStep(0.02);
        getController().setSleepPeriod(1);
        getController().addActivity(new ActivityIntegrate(integrator));
        box.setNMolecules(species, numAtoms);
        new BoxInflate(box, space, density).actionPerformed();

        Potential2Soft potential = P2LennardJones.makeTruncated(space, sigma, 1.0, new TruncationFactoryForceShift(space, rc));

        AtomType leafType = species.getLeafType();

        potentialMaster.setPairPotential(leafType, leafType, potential);
        if (!useNbrLists) {
            BoxImposePbc imposepbc = new BoxImposePbc(space);
            imposepbc.setBox(box);
            integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc));
        }

        ConfigurationLattice configuration = new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space);
        configuration.initializeCoordinates(box);
        energy = new MeterPotentialEnergyFromIntegratorFasterer(integrator);
        avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        pump = new DataPumpListener(energy, avgEnergy, 10);
        integrator.getEventManager().addListener(pump);
    }

    public static void main(String[] args) {
        final String APP_NAME = "LJMD3D";
        final LJMDFasterer sim = new LJMDFasterer(Space3D.getInstance(), 0.8, 5000, true);
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.getController().getDataStreamPumps().add(sim.pump);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
        display.setAccumulator(sim.avgEnergy);
        simGraphic.add(display);
    }
}
