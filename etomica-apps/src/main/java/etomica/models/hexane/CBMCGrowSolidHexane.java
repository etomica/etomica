/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.api.*;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Contains specifics for the solid hexane model:
 * <ul>
 * <li>bond length = 0.4
 * <li>rejects moving the whole molecule as a CBMC move
 * </ul>
 * 
 * @author cribbin
 * 
 */
public class CBMCGrowSolidHexane extends CBMCGrowStraightAlkane {

    public CBMCGrowSolidHexane(PotentialMaster p, IRandom random,
                               Space _space,
                               IntegratorMC integrator, Box phs, ISpecies species, int NTrials) {
        super(p, random, integrator, phs, species, _space, 6, NTrials);

        if (_space.D() != 3) {
            throw new IllegalArgumentException("Torsional bond is only "
                    + "used in 3D simulations");
        }

        setBondLength(0.4);

        phi = ((SpeciesHexane)species).getBondAngle();
        lowerTorsLimit = ((SpeciesHexane)species).getLowerLimit();
        upperTorsLimit = ((SpeciesHexane)species).getUpperLimit();
        
//        phi = (180 - 109.47) / 360.0 * 2.0 * Math.PI; // makes sure the vector
        // is pointing in the right direction on the cosine section

//        lowerTorsLimit = 108.6919204 / 360.0 * 2.0 * Math.PI;
//        upperTorsLimit = 251.3080797 / 360.0 * 2.0 * Math.PI;
        limit = (upperTorsLimit - lowerTorsLimit) / 2;

        rotor = _space.makeTensor();
        temp2 = _space.makeVector();
    }

    protected int calcStartIndex() {
        // function will not pick chainlength;
        return random.nextInt(chainlength - 2) + 1;
    }

    // Different because we know the bond angle
    // All moves are accepted
    protected Vector calcRandomBondWithAngle(IAtom a, IAtom b) {
        // temp will be the radial vector
        // vex will be the axial vector

        temp.E(b.getPosition());
        temp.ME(a.getPosition());
        vex.E(temp); // Store it because we need it, and we are going to
        // change temp in the next line of code.

        temp.E(getNormal(temp));

        // Subtract the projection of b-c(temp) onto a-b(vex) from temp,
        // which leaves the radial part of bc as temp
        vex.normalize();
        // This is the projection bit. (equivalent to multiplying by the cosine
        // of the angle between the radial and axial vectors)
        vex.TE(vex.dot(temp));
        vex.TE(1 / getBondlength() / getBondlength());
        // This is the subtraction bit.
        temp.ME(vex);

        // Create the rotation matrix for an arbitrary unit vector and an
        // arbitrary angle, and apply it
        double randomAngle = random.nextDouble() * 2.0 * Math.PI;
        double cosRA = Math.cos(randomAngle);
        double sinRA = Math.sin(randomAngle);

        rotor.E(new double[][] {{ cosRA + (1 - cosRA) * vex.getX(0) * vex.getX(0), // xx
                (1 - cosRA) * vex.getX(0) * vex.getX(1) - sinRA * vex.getX(2), // xy
                (1 - cosRA) * vex.getX(0) * vex.getX(2) + sinRA * vex.getX(1)}, // xz
                {(1 - cosRA) * vex.getX(1) * vex.getX(0) + sinRA * vex.getX(2), // yx
                cosRA + (1 - cosRA) * vex.getX(1) * vex.getX(1), // yy
                (1 - cosRA) * vex.getX(1) * vex.getX(2) - sinRA * vex.getX(0)}, // yz
                {(1 - cosRA) * vex.getX(2) * vex.getX(0) - sinRA * vex.getX(1), // zx
                (1 - cosRA) * vex.getX(2) * vex.getX(1) + sinRA * vex.getX(0), // zy
                cosRA + (1 - cosRA) * vex.getX(2) * vex.getX(2)} // zz
        });

        // Mulitply the rotation tensor by temp to get the rotated vector.
        // We then have a vector perpendicular to the a-b axis, rotated a random
        // angle. This is the radial vector.
        rotor.transform(temp);

        // we normalize our radial vector (temp). Our axial vector (vex) is
        // already normalized.
        temp.normalize();

        // we multiply the axial vector by the axial length and the radial
        // vector by the radial length (worked out by hand, based on a bond
        // length of 0.4 and a bond angle of 109.47
        temp.TE(getBondlength() * Math.cos(phi));
        vex.TE(getBondlength() * Math.sin(phi));

        // Add these together, and return that value:
        temp.PE(vex);
        return temp;

    }

    protected Vector calcRandomBondWithAngleAndTorsion(IAtom a, IAtom b,
                                                       IAtom c) {
        // Get a random number, and place it between the limits on the new
        // atom's placement. The angle must be between lowerTorsLimit,
        // and upperTorsLimit.
        double randomAngle = 2.0 * limit * random.nextDouble() - limit;

        double cosRA = Math.cos(randomAngle);
        double sinRA = Math.sin(randomAngle);

        // get a normal vector to the a-b vector and the b-c vector
        // This vector is, by definition, perpendicular to the a-b vector, which
        // makes it a radius of a circle centered on that axis.
        vex.E(a.getPosition());
        vex.ME(b.getPosition());
        temp.E(b.getPosition());
        temp.ME(c.getPosition());

        // Create the rotation matrix for an arbitrary unit vector
        vex.normalize();
        rotor.E(new double[][] { {cosRA + (1 - cosRA) * vex.getX(0) * vex.getX(0), // xx
                (1 - cosRA) * vex.getX(0) * vex.getX(1) - sinRA * vex.getX(2), // xy
                (1 - cosRA) * vex.getX(0) * vex.getX(2) + sinRA * vex.getX(1)}, // xz
                {(1 - cosRA) * vex.getX(1) * vex.getX(0) + sinRA * vex.getX(2), // yx
                cosRA + (1 - cosRA) * vex.getX(1) * vex.getX(1), // yy
                (1 - cosRA) * vex.getX(1) * vex.getX(2) - sinRA * vex.getX(0)}, // yz
                {(1 - cosRA) * vex.getX(2) * vex.getX(0) - sinRA * vex.getX(1), // zx
                (1 - cosRA) * vex.getX(2) * vex.getX(1) + sinRA * vex.getX(0), // zy
                cosRA + (1 - cosRA) * vex.getX(2) * vex.getX(2)} // zz
        });

        // Mulitply the rotation tensor by temp to get the rotated vector.
        rotor.transform(temp);

        // we normalize our vector (temp) and multiply it by the known
        // bondlength.
        temp.normalize();
        temp.TE(getBondlength());

        return temp;
    }

    protected double calcBondAngleEnergy(double d) {
        throw new RuntimeException("calcBondAngleEnergy should not be called "
                + "in CBMCGrowSolidHexane");
        /*
         * We are using a set bond length, and a rigid bond angle, we don't need
         * to calculate an energy because we wrote this program to only use the
         * proper bond length and angle
         */
    }

    protected double calcBondTorsionalEnergy(Vector v) {
        throw new RuntimeException("calcBondTorsionalEnergy should not be "
                + "called in CBMCGrowSolidHexane");
        /*
         * We are using a set bond length, and a rigid bond angle, we don't need
         * to calculate an energy because we wrote this program to only use the
         * proper bond length and angle
         */
    }

    protected double calcExternalEnergy(IAtom a) {
        // make sure that the proper potential is enabled. Really,
        // we only have the one potential, so this line is unnecessary
        // but I want it in here for reference when I am extending
        // this code.
        // potMast.setEnabled(pots[0], true);
        externalMeter.setBox(box);
        externalMeter.setTarget(a);

        double blind = externalMeter.getDataAsScalar();
        return blind;
    }

    // Since the bondlength is constant, this returns an ordinary getter
    protected double calcBondL() {
        return getBondlength();
    }

    /**
     * Returns a unit normal to the argument vector
     * 
     * @param vect
     * @return a unit normal to the argument vector
     */
    protected Vector getNormal(Vector vect) {
        // Determine the smallest component
        int min = 0;
        if (vect.getX(1) < vect.getX(0)) {
            min = 1;
        }
        if (vect.getX(2) < vect.getX(min)) {
            min = 2;
        }

        // create the unit vector in that direction
        temp2.E(0.0);
        temp2.setX(min, 1.0);
        return temp2;
    }

    public void setBondLength(double d) {
        bondlength = d;
    }

    public double getBondlength() {
        return bondlength;
    }

    public double energyChange() {
        return 0.0;
    }

    private static final long serialVersionUID = 1L;

    private double phi, lowerTorsLimit, upperTorsLimit, limit;

    Tensor rotor;

    Vector temp2;
}
