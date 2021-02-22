/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ensembles;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayTextBox;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.integrator.mcmove.MCMoveInsertDeleteFasterer;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolumeFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.PotentialMasterCellFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.compute.PotentialCompute;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

public class LJMCFasterer extends Simulation {

    public final SpeciesGeneral species;
    public final Box box;

    public final IntegratorMCFasterer integrator;
    public final MCMoveAtomFasterer mcMoveAtom;
    public final MCMoveVolumeFasterer mcMoveVolume;
    public final MCMoveInsertDeleteFasterer mcMoveID;
    protected DisplayTextBox vBox;
    protected DisplayTextBox nBox;
    private final Color gray = new Color(242,242,242);

    public LJMCFasterer(Space _space) {
        super(_space);
        setRandom(new RandomMersenneTwister(2));
        //species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);//index 1
        addSpecies(species);

        box = this.makeBox();
        PotentialMasterCellFasterer potentialMaster = new PotentialMasterCellFasterer(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        int N = 200;  //number of atoms

        integrator = new IntegratorMCFasterer(potentialMaster, random, 1.0, box);
        getController().addActivity(new ActivityIntegrate(integrator));

        //instantiate several potentials for selection in combo-box
        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, 2.5);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2Truncated);

        //construct box
        Vector dim = space.makeVector();
        dim.E(space.D() == 2 ? 15 : 10);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);

        mcMoveAtom = new MCMoveAtomFasterer(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMaxAdjustInterval(50000);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMinAdjustStep(1.05);

        mcMoveVolume = new MCMoveVolumeFasterer(integrator, random, 1) {
            public boolean doTrial() {
                double vOld = box.getBoundary().volume();
                uOld = integrator.getPotentialEnergy();
                hOld = uOld + pressure * vOld;
                biasOld = vBias.f(vOld);
                vScale = (2. * random.nextDouble() - 1.) * stepSize;
                vNew = vOld * Math.exp(vScale); //Step in ln(V)
                double L = Math.pow(vNew, 1.0 / space.D());
                if (L < 6.05 || (box.getLeafList().size() / vNew < 0.00011 && L > 48) || L > 148) {
                    vBox.setBackground(Color.PINK);
                } else {
                    vBox.setBackground(gray);
                }
                if (L < 6 || (box.getLeafList().size() / vNew < 0.0001 && L > 50) || L > 150) return false;
                double rScale = Math.exp(vScale / D);
                inflate.setScale(rScale);
                inflate.actionPerformed();
                PotentialCompute potentialCompute = integrator.getPotentialCompute();
                potentialCompute.init();
                uNew = potentialCompute.computeAll(false);
                hNew = uNew + pressure * vNew;
                return true;
            }
        };

        mcMoveID = new MCMoveInsertDeleteFasterer(potentialMaster, random, space) {
            public void acceptNotify() {
                if (moleculeList.size() > 999 && insert) {
                    rejectNotify();
                    uNew = 0;
                    return;
                }
                if (moleculeList.size() >= 998) {
                    nBox.setBackground(Color.PINK);
                } else {
                    nBox.setBackground(gray);
                }
                super.acceptNotify();
            }
        };
        mcMoveID.setSpecies(species);
    }

    public void setvBox(DisplayTextBox vBox) {
        this.vBox = vBox;
    }
    public void setnBox(DisplayTextBox nBox) {
        this.nBox = nBox;
    }


    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        }

        LJMCFasterer sim = new LJMCFasterer(space);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }
}
