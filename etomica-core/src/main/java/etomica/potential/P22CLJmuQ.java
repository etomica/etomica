/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;




/**
 * Two-centered Lennard Jones heteronuclear molecule with a mu2 and quadrupole.
 *
 * @author Jayant K. Singh
 */
public class P22CLJmuQ extends PotentialMolecular implements PotentialMolecularSoft{

    private double sigma1, sigma2, sigma12, sigma1Sq, sigma2Sq, sigma12Sq;
    private double epsilon1, epsilon2, epsilon12;
    private double hsdiasq = 1.0 / Math.sqrt(2);
    private double Q2, mu2;
    private Boundary boundary;
    private double siteFraction;

    private final Vector com1;
    private final Vector com2;
    private final Vector v12;
    private final Vector v1;
    private final Vector v2;
    private final Vector dr;

    public P22CLJmuQ(Space space, double sigma1, double sigma2, double epsilon1, double epsilon2, double mu, double Q, double siteFraction) {
        super(2, space);
        this.sigma1 = sigma1;
        sigma1Sq = sigma1 * sigma1;
        this.sigma2 = sigma2;
        sigma2Sq = sigma2 * sigma2;
        sigma12 = (sigma1 + sigma2) / 2;
        sigma12Sq = sigma12 * sigma12;
        this.epsilon1 = epsilon1;
        this.epsilon2 = epsilon2;
        epsilon12 = Math.sqrt(epsilon1 * epsilon2);
        this.mu2 = mu*mu;
        this.Q2 = Q*Q;
        this.siteFraction = siteFraction;
        com1 = space.makeVector();
        com2 = space.makeVector();
        v12 = space.makeVector();
        v1 = space.makeVector();
        v2 = space.makeVector();
        dr = space.makeVector();
    }

    public void setHardCoreDiamterSq(double val) {
        hsdiasq = val;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IMoleculeList pair) {

        IAtomList mol1 = pair.get(0).getChildList();
        IAtomList mol2 = pair.get(1).getChildList();
        IAtom bead11 = mol1.get(0);
        IAtom bead12 = mol1.get(1);

        IAtom bead21 = mol2.get(0);
        IAtom bead22 = mol2.get(1);

        // LJ contribution

        dr.Ev1Mv2(bead11.getPosition(), bead21.getPosition());
        boundary.nearestImage(dr);
        double ener = epsilon1 * calculateEnergy(dr.squared() / sigma1Sq);

        dr.Ev1Mv2(bead11.getPosition(), bead22.getPosition());
        boundary.nearestImage(dr);
        ener += epsilon12 * calculateEnergy(dr.squared() / sigma12Sq);

        dr.Ev1Mv2(bead12.getPosition(), bead21.getPosition());
        boundary.nearestImage(dr);
        ener += epsilon12 * calculateEnergy(dr.squared() / sigma12Sq);

        dr.Ev1Mv2(bead12.getPosition(), bead22.getPosition());
        boundary.nearestImage(dr);
        ener += epsilon2 * calculateEnergy(dr.squared() / sigma2Sq);

        if (Q2 != 0.0 && !Double.isInfinite(ener)) {

            com1.Ea1Tv1(siteFraction, bead12.getPosition());
            com1.PEa1Tv1((1 - siteFraction), bead11.getPosition());

            com2.Ea1Tv1(siteFraction, bead22.getPosition());
            com2.PEa1Tv1((1 - siteFraction), bead21.getPosition());
            v12.Ev1Mv2(com1, com2);
            boundary.nearestImage(v12);
            double r2 = v12.squared();
            double r = Math.sqrt(r2);
            v12.TE(1 / r);

            // axis one
            dr.Ev1Mv2(bead12.getPosition(), bead11.getPosition());
            boundary.nearestImage(dr);

            v1.E(dr);

            v1.normalize();

            //axis two
            dr.Ev1Mv2(bead22.getPosition(), bead21.getPosition());
            boundary.nearestImage(dr);

            v2.E(dr);

            v2.normalize();

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and v12
            double cos1 = v1.dot(v12);
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and v12
            double cos2 = v2.dot(v12);
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to v12
            // between the molecules multiplied by sin1 and sin2

            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double temp = v1.dot(v2) - 5 * cos1 * cos2;

            double uqq = (3.0 * Q2) / (4.0 * r2 * r2 * r);
            double uQuad = (uqq) * (1.0 - 5.0 * (cos1 * cos1 + cos2 * cos2 + 3 * cos1 * cos1 * cos2 * cos2) + 2 * (temp * temp));
            ener += uQuad;


            double uDipole = (v1.dot(v2) - 3.0 * v1.dot(v12) * v2.dot(v12));
            uDipole  = mu2 * uDipole / (r*r*r) ;
            ener += uDipole;


        }
        return ener;
    }

    protected double calculateEnergy(double r2) {
        if (r2 < hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        double s2 = 1 / r2;
        double s6 = s2 * s2 * s2;
        return 4 * s6 * (s6 - 1.0);
    }

    protected double calculateDU(double r2) {
        if (r2 < hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        double s2 = 1 / r2;
        double s6 = s2 * s2 * s2;
        return -4 * s6 * (12*s6 - 6) * s2 ;
    }

    @Override
    public double virial(IMoleculeList pair) {
        IAtomList mol1 = pair.get(0).getChildList();
        IAtomList mol2 = pair.get(1).getChildList();
        IAtom bead11 = mol1.get(0);
        IAtom bead12 = mol1.get(1);

        IAtom bead21 = mol2.get(0);
        IAtom bead22 = mol2.get(1);

        com1.Ea1Tv1(siteFraction, bead12.getPosition());
        com1.PEa1Tv1((1 - siteFraction), bead11.getPosition());

        com2.Ea1Tv1(siteFraction, bead22.getPosition());
        com2.PEa1Tv1((1 - siteFraction), bead21.getPosition());
        v12.Ev1Mv2(com1, com2);
        boundary.nearestImage(v12);

        // LJ contribution

        dr.Ev1Mv2(bead11.getPosition(), bead21.getPosition());
        boundary.nearestImage(dr);
        double vir = epsilon1 * calculateDU(dr.squared() / sigma1Sq) * dr.dot(v12);

        dr.Ev1Mv2(bead11.getPosition(), bead22.getPosition());
        boundary.nearestImage(dr);
        vir += epsilon12 * calculateDU(dr.squared() / sigma12Sq) * dr.dot(v12);

        dr.Ev1Mv2(bead12.getPosition(), bead21.getPosition());
        boundary.nearestImage(dr);
        vir += epsilon12 * calculateDU(dr.squared() / sigma12Sq) * dr.dot(v12);

        dr.Ev1Mv2(bead12.getPosition(), bead22.getPosition());
        boundary.nearestImage(dr);
        vir +=  epsilon2 * calculateDU(dr.squared() / sigma2Sq) * dr.dot(v12);

        if (Q2 != 0.0) {

            double r2 = v12.squared();
            double r = Math.sqrt(r2);
            v12.TE(1 / r);

            // axis one
            dr.Ev1Mv2(bead12.getPosition(), bead11.getPosition());
            boundary.nearestImage(dr);

            v1.E(dr);

            v1.normalize();

            //axis two
            dr.Ev1Mv2(bead22.getPosition(), bead21.getPosition());
            boundary.nearestImage(dr);

            v2.E(dr);

            v2.normalize();

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and v12
            double cos1 = v1.dot(v12);
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and v12
            double cos2 = v2.dot(v12);
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to v12
            // between the molecules multiplied by sin1 and sin2

            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double temp = v1.dot(v2) - 5 * cos1 * cos2;

            double uqq = (3.0 * Q2) / (4.0 * r2 * r2 * r);
            double vQuad = -5 * (uqq) * (1.0 - 5.0 * (cos1 * cos1 + cos2 * cos2 + 3 * cos1 * cos1 * cos2 * cos2) + 2 * (temp * temp));
            vir += vQuad;

            double uDipole = (v1.dot(v2) - 3.0 * v1.dot(v12) * v2.dot(v12));
            double vDipole  = -3*mu2 * uDipole / (r*r2) ;
            vir += vDipole;
        }
        return vir;
    }

    @Override
    public Vector[] gradient(IMoleculeList molecules) {
        throw new RuntimeException("nope");
    }

    @Override
    public Vector[] gradient(IMoleculeList molecules, Tensor pressureTensor) {
        throw new RuntimeException("nope");
    }
}
