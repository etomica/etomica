/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.action.IAction;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.iterator.MoleculeIteratorAllMolecules;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * 
 * @author cribbin
 * 
 */
public class CheckCBMCHexane implements IAction {

    public CheckCBMCHexane(Box p, PotentialMaster potentialMaster, Space space) {
        box = p;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        energyMeter.setBox(box);

        moleculeIterator = new MoleculeIteratorAllMolecules(box);
        moleculeIterator.reset();

        vex = space.makeVector();
        temp = space.makeVector();
        axial = space.makeVector();
        radial = space.makeVector();
        booink = 0;

        phi = 109.47 * 2.0 * Math.PI / 360.0;
        atom1 = -1;
        atom2 = -2;
    }

    public void actionPerformed() {
        // System.out.println("startcheck #" + booink);
        double tol = 0.000005;
        double test, ang;

        if (energyMeter.getDataAsScalar() != 0.0) {
            throw new RuntimeException("Non-zero potential energy!");
        }
        // System.out.println("NRG chk");

        // Check that bond lengths are 0.4;
        IMoleculeList moleculeList = box.getMoleculeList();
        for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList atomList = molecule.getChildList();
            for (int i = 0; i < atomList.getAtomCount() - 1; i++) {
                // vex.E(((AtomLeaf)atomList.get(i)).getPosition());
                vex.ME(atomList.getAtom(i + 1).getPosition());
                test = Math.sqrt(vex.squared());
                test -= length;
                if (Math.abs(test) > tol) {
                    System.out.println(atomList.getAtom(i));
                    throw new RuntimeException("The bond length is not "
                            + length + "!");
                }
            }
        }
        // System.out.println("length chk");

        // Check that bond angles are 109.47 degrees
        tol = 0.0000005;
        for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList atomList = molecule.getChildList();
            for (int i = 0; i < atomList.getAtomCount() - 2; i++) {
                vex.E(atomList.getAtom(i).getPosition());
                vex.ME(atomList.getAtom(i + 1).getPosition());
                temp.E(atomList.getAtom(i + 2).getPosition());
                temp.ME(atomList.getAtom(i + 1).getPosition());

                ang = Math.acos(vex.dot(temp) / length / length);
                ang -= phi;
                if (Math.abs(ang) > tol) {
                    System.out.println(atomList.getAtom(i));
                    throw new RuntimeException("The bond angle is bad.");
                }
            }
        }
        // System.out.println("bond ang chk");

        // Check that torsional bond angles are between ~108 and ~251
        tol = 0.0000000001;
        double makeGood;
        for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList atomList = molecule.getChildList();
            for (int i = 0; i < atomList.getAtomCount() - 3; i++) {
                vex.E(atomList.getAtom(i).getPosition());
                vex.ME(atomList.getAtom(i + 1).getPosition());
                temp.E(atomList.getAtom(i + 3).getPosition());
                temp.ME(atomList.getAtom(i + 2).getPosition());
                axial.E(atomList.getAtom(i + 2).getPosition());
                axial.ME(atomList.getAtom(i + 1).getPosition());

                // Project each vector onto the axial vector, and subtract the
                // axial portion from the result, leaving the radial portion
                axial.normalize();
                axial.TE(length * Math.cos((Math.PI - phi))); // (Pi - phi is
                                                                // 180 - 109.47,
                                                                // or the acute
                                                                // angle)
                vex.ME(axial);
                axial.TE(-1.0); // we do this because we have assumed that
                // the angle between the two vectors (axial
                // & whatever) is obtuse. The angle between
                // axial and temp is acute, so we reverse
                // axial to make the angle acute.
                temp.ME(axial);

                vex.normalize();
                temp.normalize();

                // Calculate the angle between the projected vectors

                /*
                 * We have to do some stuff to "fix" answers that are just over
                 * 1.0 or under -1.0 because we get NaN, and sometimes we are
                 * just really close (like, into the pico and femto range) and
                 * we can make that happy.
                 */
                makeGood = vex.dot(temp);
                if (makeGood < -1) {
                    if ((makeGood + tol) > -1) {
                        makeGood = -1.0;
                    }
                }
                if (makeGood > 1) {
                    if ((makeGood - tol) < 1) {
                        makeGood = 1.0;
                    }
                }

                ang = Math.acos(makeGood);

                if (ang < lowerTorsLimit) {
                    System.out.println(atomList.getAtom(i));
                    throw new RuntimeException("Torsional angle below limit!");
                }
                if (upperTorsLimit < ang) {
                    System.out.println(atomList.getAtom(i));
                    throw new RuntimeException("Torsional angle above limit!");
                }
            }
        }
        // System.out.println("tors ang chk");

        // System.out.println("end check #" + booink);

        System.out.println(booink);
        booink++;
    }

    MeterPotentialEnergy energyMeter;

    Box box;

    MoleculeIteratorAllMolecules moleculeIterator;

    Vector vex, temp, axial, radial;

    int booink;

    double phi, length, lowerTorsLimit, upperTorsLimit;

    double actualAngle;

    int atom1, atom2;

    public void setLength(double length) {
        this.length = length;
    }

    public void setLowerTorsLimit(double lowerTorsLimit) {
        this.lowerTorsLimit = lowerTorsLimit;
    }

    public void setPhi(double phi) {
        this.phi = phi;
    }

    public void setUpperTorsLimit(double upperTorsLimit) {
        this.upperTorsLimit = upperTorsLimit;
    }
}
