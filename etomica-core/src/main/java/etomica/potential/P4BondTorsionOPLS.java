/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.RandomPositionSourceRectangular;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.util.random.IRandom;
import etomica.util.random.RandomNumberGenerator;

/**
 * OPLS Torsion potential.
 *
 * @author Andrew Schultz
 */
public class P4BondTorsionOPLS extends P4BondTorsion {

    public P4BondTorsionOPLS(Space space, double a0, double a1, double a2, double a3) {
        super(space, a0, a1, a2, a3);
    }

    @Override
    public double u(double costheta) {
        double cos2theta = 2 * costheta * costheta - 1;
        double cos3theta = costheta * (2 * cos2theta - 1);
        double cos4theta = 2 * cos2theta * cos2theta - 1;

        return 0.5 * a0 * (1 + costheta) + 0.5 * a1 * (1 - cos2theta) + 0.5 * a2 * (1 + cos3theta) + 0.5 * a3 * (1 - cos4theta);
    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        throw new RuntimeException("Implement me");
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P4BondTorsionOPLS potential = new P4BondTorsionOPLS(space, 0, 10, 20, 30);
        IRandom random = new RandomNumberGenerator();
        Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
        RandomPositionSourceRectangular positionSource = new RandomPositionSourceRectangular(space, random);
        positionSource.setBox(box);
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        Atom atom3 = new Atom(space);
        AtomArrayList atoms = new AtomArrayList(4);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        atoms.add(atom3);
        int n = 40;
        Vector gradient = space.makeVector();
        Vector dr = space.makeVector();
        for (int i=0; i<n; i++) {
            atom0.getPosition().E(positionSource.randomPosition());
            atom1.getPosition().E(positionSource.randomPosition());
            atom2.getPosition().E(positionSource.randomPosition());
            atom3.getPosition().E(positionSource.randomPosition());
            
            double U = potential.energy(atoms);

            int iRand = random.nextInt(4);
            IAtom atom = atoms.get(iRand);
            gradient.E(potential.gradient(atoms)[iRand]);
            
            dr.setRandomSphere(random);
            dr.TE(0.0001);
            double expectedDeltaU = gradient.dot(dr);
            
            atom.getPosition().PE(dr);
            
            double newU = potential.energy(atoms);
            
            System.out.println(expectedDeltaU+" "+(newU-U)+" "+(expectedDeltaU-newU+U));
        }
    }
}
