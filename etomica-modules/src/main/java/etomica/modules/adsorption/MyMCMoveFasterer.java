/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.mcmove.MCMoveInsertDeleteFasterer;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeArrayList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * @author kofke
 * <p>
 * Extends insert/delete mcmove class by permitting insertion/deletion in a
 * subregion defined by a range of values of the z coordinate.
 */
public class MyMCMoveFasterer extends MCMoveInsertDeleteFasterer implements IntegratorListener {

    private double zFraction, sigma;
    private final Vector position;
    private final MoleculeArrayList activeAtoms;
    private final AtomActionRandomizeVelocity randomizer;
    private final IntegratorBoxFasterer integrator;
    protected int testMoleculeIndex;
    protected DoubleRange range;
    protected final int dim;

    public MyMCMoveFasterer(IntegratorBoxFasterer integrator, IRandom random,
                            Space space, double zFraction, double sigma, int dim) {
        super(integrator.getPotentialCompute(), random, space);
        position = space.makeVector();
        setZFraction(zFraction, sigma);
        this.dim = dim;
        this.integrator = integrator;
        randomizer = new AtomActionRandomizeVelocity(0, random);
        activeAtoms = new MoleculeArrayList();
    }

    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = (random.nextInt(2) == 0);
//        System.out.println((insert?"insert":"delete")+" trial");
        if (insert) {
            uOld = 0.0;
            if (!reservoir.isEmpty()) {
                testMolecule = reservoir.remove(reservoir.size() - 1);
            } else {
                testMolecule = species.makeMolecule();
            }
            position.E(positionSource.randomPosition());
            double z = position.getX(dim);
            double zBoundary = box.getBoundary().getBoxSize().getX(dim);
            // in theory we insert into the top zFraction of the box
            // but we need to avoid the top sigma as well
            z = z * zFraction + zBoundary * (-zFraction / 2 + 0.5);
//            System.out.println("inserting at "+z);
            position.setX(dim, z);
            atomTranslator.setDestination(position);
            atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
        } else {//delete
            if (activeAtoms.size() == 0) {
                testMolecule = null;//added this line 09/19/02
                return false;
            }
            testMoleculeIndex = random.nextInt(activeAtoms.size());
            testMolecule = activeAtoms.get(testMoleculeIndex);
            uOld = potentialMaster.computeOneOldMolecule(testMolecule);
        }
        uNew = Double.NaN;
        return true;
    }//end of doTrial

    public double getChi(double temperature) {//note that moleculeCount() gives the number of molecules after the trial is attempted
        if (insert) {
            if (testMolecule.getChildList().get(0).getLeafIndex() == 34) {
                System.out.println("energy for 34");
                IAtomList atoms = box.getLeafList();
                System.out.println("29 " + atoms.get(29).getPosition());
                System.out.println("34 " + atoms.get(34).getPosition());
            }
            uNew = potentialMaster.computeOneMolecule(testMolecule);
        } else {
            uNew = 0;
        }
        double b = uOld - uNew;
        if (insert) b += mu;
        else b -= mu;

        double a = insert ? zFraction * box.getBoundary().volume() / (activeAtoms.size() + 1)
                : activeAtoms.size() / zFraction / box.getBoundary().volume();
        return a * Math.exp(b / temperature);
    }

    public void rejectNotify() {
//        if(!insert) {
//            System.out.println("rejected deleting "+testMolecule.getIndex());
//        }
//        else {
//            System.out.println("rejected inserting "+testMolecule.getIndex());
//        }
        super.rejectNotify();
    }

    public void acceptNotify() {
//        System.out.println("accept "+(insert ? "insert" : "delete")+" "+testMolecule.getChildList().get(0).getLeafIndex());
        if (!insert && testMolecule.getChildList().get(0).getLeafIndex() == 2 && box.getLeafList().size() == 35) {
            IAtomList atoms = box.getLeafList();
            System.out.println("2 " + atoms.get(2).getPosition());
            System.out.println("29 " + atoms.get(29).getPosition());
            System.out.println("last(" + (atoms.size() - 1) + ") " + atoms.get(atoms.size() - 1).getPosition());
        }
        super.acceptNotify();
        if (!insert) {
//	        System.out.println("accepted deleting "+testMolecule.getIndex());
            activeAtoms.remove(testMoleculeIndex);
        } else {
//            System.out.println("accepted inserting "+testMolecule.getIndex());
            activeAtoms.add(testMolecule);
            randomizer.setTemperature(integrator.getTemperature());
            randomizer.actionPerformed(testMolecule.getChildList().get(0));
        }
    }

    public void setupActiveAtoms() {
        activeAtoms.clear();
        double zBoundary = box.getBoundary().getBoxSize().getX(dim);
        int nMolecules = moleculeList.size();
        for (int i = 0; i < nMolecules; i++) {
            IMolecule molecule = moleculeList.get(i);
            if (molecule.getType() != species) continue;

            double z = molecule.getChildList().get(0).getPosition().getX(dim);
            if (z < zBoundary * (0.5 - zFraction)) continue;
            activeAtoms.add(molecule);
        }
    }

    /**
     * Sets the zFraction, the fraction of the box volume into which atoms are
     * inserted or deleted.  The volume is on the far left or right side (in
     * the z dimension).  To put the volume on the left side (negative z),
     * specify a negative z fraction (a positive value will result in the
     * volume being on the positive z side.
     */
    public void setZFraction(double zFraction, double zFractionMax) {
        this.zFraction = zFraction;
        this.sigma = zFractionMax;
    }

    public double getZFraction() {
        return zFraction;
    }

    public void setSpecies(ISpecies s) {
        super.setSpecies(s);
        if (box != null) {
            moleculeList = box.getMoleculeList(s);
        }
    }

    @Override
    public void integratorInitialized(IntegratorEvent e) {
        setupActiveAtoms();
    }

    @Override
    public void integratorStepStarted(IntegratorEvent e) {
    }

    @Override
    public void integratorStepFinished(IntegratorEvent e) {
    }
}
