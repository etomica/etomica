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
 * Hard sphere molecule with a dipole sitting at the center.
 *
 * @author Weisong & Shu
 * date:May 2015
 *
 */
public class P2HSDQPlus implements IPotential2 {

    private double sigma , sigma2, epsilon;
    private double dipole, theta, kappa, alpha, D;
    protected final Vector[] a;

    public P2HSDQPlus(Space space, double sigma, double epsilon, double dipole, double theta, double alpha, double kappa, double D) {
        setSigma(sigma);
        this.epsilon = epsilon;
        a = new Vector[3];
        a[0] = space.makeVector();
        a[1] = space.makeVector();
        a[2] = space.makeVector();

        setDipole(dipole);
        setQuadrupole(theta);
        setAlpha(alpha);
        setKappa(kappa);
        this.D = D;
    }

    @Override
    public double du(double r2) {
        return 0;
    }

    @Override
    public double d2u(double r2) {
        return 0;
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        double r2 = dr12.squared();

        if(r2 < 0.1*sigma2) { // hard core
            return Double.POSITIVE_INFINITY;
        }
        // normalize dr, the vector between the molecules
        double r = Math.sqrt(dr12.squared());

        // v1 (unit vector) is the orientation of molecule 1: dipole1 direction
        Vector v1 = ((IAtomOriented)atom1).getOrientation().getDirection();
        // v2 (unit vector) is the orientation of molecule 2: dipole1 direction
        Vector v2 = ((IAtomOriented)atom2).getOrientation().getDirection();

        // cos(dipole 1 and dipole 2)=cos(v1 and v2)
        double cos_D1_D2 = v1.dot(v2);
        //cos(dipole 1 and r12)
        double c1 = v1.dot(dr12)/r;
        //cos(r12 and dipole 2)
        double c2 = dr12.dot(v2)/r;
        double s1 = Math.sqrt(1 - c1*c1);
        double s2 = Math.sqrt(1 - c2*c2);
        // 0.5 si sj Cij - ci cj = cos_D1_D2 - 3 ci cj
        // Cij = 2 (cos_D1_D2 - 2 ci cj) / si sj
        double Cij = 2*(cos_D1_D2 - 2*c1*c2) / (s1 * s2);
        double umumu = dipole * dipole *  (cos_D1_D2 - 3.0  * c1* c2) / (r*r2);
        double umuQ = 3*dipole*theta/(2*r2*r2) * (c2 - c1) * (3*c1*c2 - 2*s1*s2*Cij + 1);
        double uQQ = 3*theta*theta/(4*r*r2*r2) * (1 - 5*c1*c1 - 5*c2*c2 + 17*c1*c2
                - 16*c1*c2*Cij - 2*s1*s1*s2*s2*Cij*Cij);
        double uDi = -dipole*dipole*alpha/(2*r2*r2*r2)*(3*c1*c1 + 3*c2*c2 + 2);
        double uQi = -9*theta*theta*alpha/(8*r2*r2*r2*r2)*(s1*s1*s1*s1 + 4*c1*c1*c1*c1 + s2*s2*s2*s2 + 4*c2*c2*c2*c2);
        double s = sigma/r;
        double s6 = s*s*s*s*s*s;
        double uRepAni = 4*D*epsilon*s6*s6*(3*c1*c1 + 3*c2*c2 - 2);
        double uAttAni = 4*epsilon*s6*(kappa - 1.5*kappa*(1-kappa)*(c1*c1+c2*c2) - 1.5*kappa*kappa*(s1*c1*Cij - 2*c1*c2)*(s1*c1*Cij - 2*c1*c2));

        return umumu + umuQ + uQQ + uDi + uQi + uRepAni + uAttAni;
    }

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s * s;
    }

    public void setDipole(double moment){
        dipole=moment;
    }

    public void setQuadrupole(double theta){
        this.theta = theta;
    }

    public void setAlpha(double alpha){
        this.alpha = alpha;
    }

    public void setKappa(double kappa){
        this.kappa = kappa;
    }

    @Override
    public double energy(IAtomList atoms) {
        return 0;
    }

    @Override
    public double u(double r2) {
        return 0;
    }
}
