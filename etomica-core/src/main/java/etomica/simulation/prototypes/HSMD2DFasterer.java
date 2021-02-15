/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHardFasterer;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

/**
 * Hard-sphere molecular dynamics simulation of two species in 2D.
 * <p>
 * Atoms of each species have the same diameter but one has 10 times more mass than the other.
 * </p>
 *
 * @author David Kofke
 */

public class HSMD2DFasterer extends Simulation {

    public IntegratorHardFasterer integrator;
    public SpeciesGeneral species1;
    public SpeciesGeneral species2;
    public Box box;
    public P2HardGeneric potential11;
    public P2HardGeneric potential12;
    public P2HardGeneric potential22;

    public HSMD2DFasterer() {
        super(Space2D.getInstance());

        species1 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species1);
        addSpecies(species2);

        box = this.makeBox();

        double sigma = 1;
        double neighborRangeFac = 1.6;
        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, neighborRangeFac * sigma, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair potentialMaster = new PotentialComputePair(this, box, neighborManager);

        getController().setSleepPeriod(1);
        getController().addActivity(new ActivityIntegrate(integrator));
        AtomType leafType1 = species1.getLeafType();
        AtomType leafType2 = species2.getLeafType();
        ((ElementSimple) leafType2.getElement()).setMass(10);

        potential11 = P2HardSphere.makePotential(sigma);
        potential12 = P2HardSphere.makePotential(sigma);
        potential22 = P2HardSphere.makePotential(sigma);

        potentialMaster.setPairPotential(leafType1, leafType1, potential11);
        potentialMaster.setPairPotential(leafType1, leafType2, potential12);
        potentialMaster.setPairPotential(leafType2, leafType2, potential22);

        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), neighborManager, random, 0.01, 1, box);
        integrator.setIsothermal(false);

        box.setNMolecules(species1, 512);
        box.setNMolecules(species2, 5);

        BoxInflate bi = new BoxInflate(box, space);
        bi.setTargetDensity(0.5);
        bi.actionPerformed();

        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2DFasterer sim = new HSMD2DFasterer();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        simGraphic.makeAndDisplayFrame();
    }

}
