/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.DataSourceCountSteps;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

/**
 * Hard sphere Monte Carlo Simulation of two species in 2D.
 * <p>
 * Species have the same sphere diameter but non-additive cross interactions.
 * (sigma12 = 1.2)
 *
 * @author David Kofke
 */

public class HSMC2D extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species, species2;
    public Box box;
    public P2HardGeneric potential11, potential12, potential22;
    public DataSourceCountSteps meterCycles;

    public HSMC2D() {
        super(Space2D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        addSpecies(species2);

        box = this.makeBox();
        PotentialMaster potentialMaster = new PotentialMaster(getSpeciesManager(), box, BondingInfo.noBonding());
        integrator = new IntegratorMC(potentialMaster, this.getRandom(), 1.0, box);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, box);
        getController().addActivity(new ActivityIntegrate(integrator));
        box.setNMolecules(species, 20);
        box.setNMolecules(species2, 20);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        potential11 = P2HardSphere.makePotential(1.0);
        potential12 = P2HardSphere.makePotential(1.2);
        potential22 = P2HardSphere.makePotential(1.0);
        AtomType type1 = species.getLeafType();
        AtomType type2 = species2.getLeafType();
        potentialMaster.setPairPotential(type1, type1, potential11);
        potentialMaster.setPairPotential(type1, type2, potential12);
        potentialMaster.setPairPotential(type2, type2, potential22);
        meterCycles = new DataSourceCountSteps(integrator);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));

    }

    public static void main(String[] args) {

        final HSMC2D sim = new HSMC2D();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        simGraphic.makeAndDisplayFrame();
    }

}
