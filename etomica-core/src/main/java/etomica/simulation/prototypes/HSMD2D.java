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
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Hard-sphere molecular dynamics simulation of two species in 2D.
 * <p>
 * Atoms of each species have the same diameter but one has 10 times more mass than the other.
 * </p>
 *
 * @author David Kofke
 */

public class HSMD2D extends Simulation {

    public IntegratorHard integrator;
    public SpeciesSpheresMono species1, species2;
    public Box box;
    public Potential2 potential11;
    public Potential2 potential12;
    public Potential2 potential22;

    public HSMD2D() {
        super(Space2D.getInstance());
        box = new Box(space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        double sigma = 1;

        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac * sigma);

        integrator = new IntegratorHard(this, potentialMaster, space, box);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        potentialMaster.setRange(sigma * 1.6);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species1 = new SpeciesSpheresMono(this, space);
        species1.setIsDynamic(true);
        species2 = new SpeciesSpheresMono(this, space);
        species2.setIsDynamic(true);
        AtomType leafType1 = species1.getLeafType();
        AtomType leafType2 = species2.getLeafType();
        addSpecies(species1);
        addSpecies(species2);
        ((ElementSimple) leafType2.getElement()).setMass(10);

        potential11 = new P2HardSphere(space, sigma, false);
        potential12 = new P2HardSphere(space, sigma, false);
        potential22 = new P2HardSphere(space, sigma, false);

        potentialMaster.addPotential(potential11, new AtomType[]{leafType1, leafType1});

        potentialMaster.addPotential(potential12, new AtomType[]{leafType2, leafType2});

        potentialMaster.addPotential(potential22, new AtomType[]{leafType1, leafType2});

        addBox(box);
        box.setNMolecules(species1, 512);
        box.setNMolecules(species2, 5);

        BoxInflate bi = new BoxInflate(box, space);
        bi.setTargetDensity(0.5);
        bi.actionPerformed();
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        simGraphic.makeAndDisplayFrame();
    }

}
