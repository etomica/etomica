/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Two-centered Lennard Jones heteronuclear molecule with a mu2 and quadrupole.
 *
 * @author Jayant K. Singh
 */
public class P22CLJmuQAtomic implements Potential2Soft {

    private double sigma12, sigma1Sq, sigma2Sq, sigma12Sq;
    private double epsilon1, epsilon2, epsilon12;
    private double hsdiasq = 1.0 / Math.sqrt(2);
    private double Q2, mu2;
    private double Q, mu;
    private double bondLength1, bondLength2;

    private final Vector dr;

    public P22CLJmuQAtomic(Space space, double sigma1, double sigma2, double epsilon1, double epsilon2, double mu, double Q, double bondLength1, double bondLength2) {
        sigma1Sq = sigma1 * sigma1;
        sigma2Sq = sigma2 * sigma2;
        sigma12 = (sigma1 + sigma2) / 2;
        sigma12Sq = sigma12 * sigma12;
        this.epsilon1 = epsilon1;
        this.epsilon2 = epsilon2;
        epsilon12 = Math.sqrt(epsilon1 * epsilon2);
        this.mu2 = mu * mu;
        this.mu = mu;
        this.Q2 = Q * Q;
        this.Q = Q;
        this.bondLength1 = bondLength1;
        this.bondLength2 = bondLength2;
        dr = space.makeVector();
    }

    public void setHardCoreDiamterSq(double val) {
        hsdiasq = val;
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
        return -4 * s6 * (12 * s6 - 6) * s2;
    }

    @Override
    public double energy(IAtomList atoms) {
        return 0;
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {

        IAtomOriented a1 = (IAtomOriented) atom1;
        IAtomOriented a2 = (IAtomOriented) atom2;

        // LJ contribution
        dr.E(dr12);
        dr.PEa1Tv1(-bondLength1, a1.getOrientation().getDirection());
        dr.PEa1Tv1(bondLength1, a2.getOrientation().getDirection());
        double ener = epsilon1 * calculateEnergy(dr.squared() / sigma1Sq);

        dr.E(dr12);
        dr.PEa1Tv1(-bondLength1, a1.getOrientation().getDirection());
        dr.PEa1Tv1(-bondLength2, a2.getOrientation().getDirection());
        ener += epsilon12 * calculateEnergy(dr.squared() / sigma12Sq);

        dr.E(dr12);
        dr.PEa1Tv1(bondLength2, a1.getOrientation().getDirection());
        dr.PEa1Tv1(bondLength1, a2.getOrientation().getDirection());
        ener += epsilon12 * calculateEnergy(dr.squared() / sigma12Sq);

        dr.E(dr12);
        dr.PEa1Tv1(bondLength2, a1.getOrientation().getDirection());
        dr.PEa1Tv1(-bondLength2, a2.getOrientation().getDirection());
        ener += epsilon2 * calculateEnergy(dr.squared() / sigma2Sq);

        if ((Q2 != 0.0 || mu != 0) && !Double.isInfinite(ener)) {

            double r2 = dr12.squared();
            double r = Math.sqrt(r2);

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and v12
            double cos1 = a1.getOrientation().getDirection().dot(dr12) / r;
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and v12
            double cos2 = a2.getOrientation().getDirection().dot(dr12) / r;
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to v12
            // between the molecules multiplied by sin1 and sin2

            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double dot = a1.getOrientation().getDirection().dot(a2.getOrientation().getDirection());
            double temp = dot - 5 * cos1 * cos2;

            double uqq = (3.0 * Q2) / (4.0 * r2 * r2 * r);
            double uQuad = (uqq) * (1.0 - 5.0 * (cos1 * cos1 + cos2 * cos2 + 3 * cos1 * cos1 * cos2 * cos2) + 2 * (temp * temp));
            ener += uQuad;

            double uDipole = (dot - 3.0 * cos1 * cos2);
            uDipole = mu2 * uDipole / (r * r * r);
            ener += uDipole;

            double temp2 = cos1 * (3 * cos2 * cos2 - 1) - 2 * cos2 * (dot - cos1 * cos2);
            double temp3 = cos2 * (3 * cos1 * cos1 - 1) - 2 * cos1 * (dot - cos1 * cos2);
            double uDipoleQual = (3 / (r2 * r2)) * (mu * Q * temp2 - mu * Q * temp3);
            ener += uDipoleQual;
        }
        return ener;
    }

    @Override
    public double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
        IAtomOriented a1 = (IAtomOriented) atom1;
        IAtomOriented a2 = (IAtomOriented) atom2;

        // LJ contribution

        dr.E(dr12);
        dr.PEa1Tv1(-bondLength1, a1.getOrientation().getDirection());
        dr.PEa1Tv1(bondLength1, a2.getOrientation().getDirection());
        double r2 = dr.squared();
        double ener = epsilon1 * calculateEnergy(r2 / sigma1Sq);
        double du = epsilon1 * calculateDU(r2 / sigma1Sq) / sigma1Sq;
        f1.PEa1Tv1(du, dr);
        f2.PEa1Tv1(-du, dr);

        dr.E(dr12);
        dr.PEa1Tv1(-bondLength1, a1.getOrientation().getDirection());
        dr.PEa1Tv1(-bondLength2, a2.getOrientation().getDirection());
        r2 = dr.squared();
        ener += epsilon12 * calculateEnergy(r2 / sigma12Sq);
        du = epsilon12 * calculateDU(r2 / sigma12Sq) / sigma12Sq;
        f1.PEa1Tv1(du, dr);
        f2.PEa1Tv1(-du, dr);

        dr.E(dr12);
        dr.PEa1Tv1(bondLength2, a1.getOrientation().getDirection());
        dr.PEa1Tv1(bondLength1, a2.getOrientation().getDirection());
        r2 = dr.squared();
        ener += epsilon12 * calculateEnergy(r2 / sigma12Sq);
        du = epsilon12 * calculateDU(r2 / sigma12Sq) / sigma12Sq;
        f1.PEa1Tv1(du, dr);
        f2.PEa1Tv1(-du, dr);

        dr.E(dr12);
        dr.PEa1Tv1(bondLength2, a1.getOrientation().getDirection());
        dr.PEa1Tv1(-bondLength2, a2.getOrientation().getDirection());
        r2 = dr.squared();
        ener += epsilon2 * calculateEnergy(r2 / sigma2Sq);
        du = epsilon2 * calculateDU(r2 / sigma2Sq) / sigma2Sq;
        f1.PEa1Tv1(du, dr);
        f2.PEa1Tv1(-du, dr);

        if (Q2 != 0.0 && !Double.isInfinite(ener)) {

            r2 = dr12.squared();
            double r = Math.sqrt(r2);

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and v12
            double cos1 = a1.getOrientation().getDirection().dot(dr12) / r;
            Vector dcos1dr12 = Vector.d(3);
            dcos1dr12.Ea1Tv1(1 / r, a1.getOrientation().getDirection());
            dcos1dr12.PEa1Tv1(-cos1 / r2, dr12);

            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and v12
            double cos2 = a2.getOrientation().getDirection().dot(dr12) / r;
            Vector dcos2dr12 = Vector.d(3);
            dcos2dr12.Ea1Tv1(1 / r, a2.getOrientation().getDirection());
            dcos2dr12.PEa1Tv1(-cos2 / r2, dr12);

            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to v12
            // between the molecules multiplied by sin1 and sin2

            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double dot = a1.getOrientation().getDirection().dot(a2.getOrientation().getDirection());
            double temp = dot - 5 * cos1 * cos2;
            Vector dtempdr12 = Vector.d(3);
            dtempdr12.Ea1Tv1(-5 * cos1, dcos2dr12);
            dtempdr12.PEa1Tv1(-5 * cos2, dcos1dr12);

            double uqq = (3.0 * Q2) / (4.0 * r2 * r2 * r);
            double rduqqdr = -5 * uqq;
            Vector duqqdr12 = Vector.d(3);
            duqqdr12.Ea1Tv1(rduqqdr / r2, dr12);
            double uQuad = (uqq) * (1.0 - 5.0 * (cos1 * cos1 + cos2 * cos2 + 3 * cos1 * cos1 * cos2 * cos2) + 2 * (temp * temp));
            Vector duQuaddr12 = Vector.d(3);
            duQuaddr12.Ea1Tv1(uQuad / uqq, duqqdr12);
            duQuaddr12.PEa1Tv1(-5 * 2 * uqq * cos1, dcos1dr12);
            duQuaddr12.PEa1Tv1(-5 * 2 * uqq * cos2, dcos2dr12);
            duQuaddr12.PEa1Tv1(-5 * 3 * 2 * uqq * cos1 * cos2 * cos2, dcos1dr12);
            duQuaddr12.PEa1Tv1(-5 * 3 * 2 * uqq * cos1 * cos1 * cos2, dcos2dr12);
            duQuaddr12.PEa1Tv1(2 * 2 * uqq * temp, dtempdr12);
            ener += uQuad;
            f1.PE(duQuaddr12);
            f2.ME(duQuaddr12);


            double uDipole = (dot - 3.0 * cos1 * cos2);
            Vector duDipoledr12 = Vector.d(3);
            duDipoledr12.Ea1Tv1(-3.0 * cos1, dcos2dr12);
            duDipoledr12.PEa1Tv1(-3.0 * cos2, dcos1dr12);
            uDipole = mu2 * uDipole / (r * r * r);
            duDipoledr12.TE(mu2 / (r * r * r));
            duDipoledr12.PEa1Tv1(-3 * uDipole / r2, dr12);
            ener += uDipole;
            f1.PE(duDipoledr12);
            f2.ME(duDipoledr12);

            double temp2 = cos1 * (3 * cos2 * cos2 - 1) - 2 * cos2 * (dot - cos1 * cos2);
            Vector dtemp2dr12 = Vector.d(3);
            dtemp2dr12.Ea1Tv1(3 * cos2 * cos2 - 1, dcos1dr12);
            dtemp2dr12.PEa1Tv1(2 * 3 * cos1 * cos2, dcos2dr12);
            dtemp2dr12.PEa1Tv1(-2 * (dot - cos1 * cos2), dcos2dr12);
            dtemp2dr12.PEa1Tv1(2 * cos2 * cos2, dcos1dr12);
            dtemp2dr12.PEa1Tv1(2 * cos2 * cos1, dcos2dr12);
            double temp3 = cos2 * (3 * cos1 * cos1 - 1) - 2 * cos1 * (dot - cos1 * cos2);
            Vector dtemp3dr12 = Vector.d(3);
            dtemp3dr12.Ea1Tv1(3 * cos1 * cos1 - 1, dcos2dr12);
            dtemp3dr12.PEa1Tv1(2 * 3 * cos2 * cos1, dcos1dr12);
            dtemp3dr12.PEa1Tv1(-2 * (dot - cos1 * cos2), dcos1dr12);
            dtemp3dr12.PEa1Tv1(2 * cos1 * cos1, dcos2dr12);
            dtemp3dr12.PEa1Tv1(2 * cos2 * cos1, dcos1dr12);
            double uDipoleQuad = (3 / (r2 * r2)) * (mu * Q * temp2 - mu * Q * temp3);
            Vector duDipoleQuad = Vector.d(3);
            duDipoleQuad.Ea1Tv1(-4 * uDipoleQuad / r2, dr12);
            duDipoleQuad.PEa1Tv1(3 / (r2 * r2) * mu * Q, dtemp2dr12);
            duDipoleQuad.PEa1Tv1(-3 / (r2 * r2) * mu * Q, dtemp3dr12);
            ener += uDipoleQuad;
            f1.PE(duDipoleQuad);
            f2.ME(duDipoleQuad);
        }
        return ener;
    }

    @Override
    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        return 0;
    }

}
