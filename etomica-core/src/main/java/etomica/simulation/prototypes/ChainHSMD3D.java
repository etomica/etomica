/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.models.ModelChain;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.graphics.BondListener;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardBond;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;

/**
 * Molecular dynamics of chains of hard spheres.
 * <p>
 * Each chain has 4 spheres, joined by a hard tether potential.
 */
public class ChainHSMD3D extends Simulation {

    public Box box;
    public IntegratorHard integrator;
    public SpeciesSpheres species;
    public P2HardSphere potential;
    private ModelChain model;

    public ChainHSMD3D() {
        super(Space3D.getInstance());
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        int numAtoms = 108;
        int chainLength = 4;
        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac);

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator, 1, true);
        getController().addAction(activityIntegrate);

        model = new ModelChain(space, true);
        model.setNumAtoms(chainLength);
        model.setBondingPotential(new P2HardBond(space, 1.0, 0.15, true));

        species = (SpeciesSpheres) model.makeSpecies(this);
        potentialMaster.addModel(model);
        ((ConformationLinear) model.getConformation()).setBondLength(1.0);
        ((ConformationLinear) model.getConformation()).setAngle(1, 0.35);

        box = new Box(space);
        double l = 14.4573 * Math.pow((chainLength * numAtoms / 2020.0), 1.0 / 3.0);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{l, l, l}));
        addBox(box);
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        box.setNMolecules(species, numAtoms);
        config.initializeCoordinates(box);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        potential = new P2HardSphere(space, 1.0, true);
        AtomType leafType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});
        CriterionBondedSimple nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        ((CriterionInterMolecular) potentialMaster.getCriterion(potential)).setIntraMolecularCriterion(nonBondedCriterion);

        integrator.setBox(box);
    }

    public static void main(String[] args) {

        final etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
        BondListener bl = new BondListener(sim.box, (DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas);
        bl.addModel(sim.model);
        ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom());
        simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame();
    }
}
