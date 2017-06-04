/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Lennard Jones molecule with a quadrupole.
 *
 * @author Jayant K. Singh
 */
public class P2LJQ extends Potential2 implements Potential2Soft {

    public P2LJQ(Space space) {
        this(space, 1, 1, 1);
    }

    public P2LJQ(Space space, double sigma, double epsilon, double momentSquared) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setQuadrupolarMomentSquare(momentSquared);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        drunit = space.makeVector();
        dcos1dr = space.makeVector();
        dcos2dr = space.makeVector();
    }

    public void setHardCoreDiamterSq(double val){
        hsdiasq=val;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IAtomList pair){
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)pair.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if(r2<hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        double s2 = sigma2/(r2);
        double s6 = s2*s2*s2;
        double ener = epsilon4*s6*(s6 - 1.0);

        if(Q2!=0.0){
            // normalize dr, the vector between the molecules
            dr.normalize();

            // v1 is the orientation of molecule 1
            Vector v1 = atom1.getOrientation().getDirection();

            // v2 is the orientation of molecule 2
            Vector v2 = atom2.getOrientation().getDirection();

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and dr
            double cos1 = v1.dot(dr);
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and dr
            double cos2 = v2.dot(dr);
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to dr
            // between the molecules multiplied by sin1 and sin2
            
            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double temp = v1.dot(v2) - 5*cos1*cos2;
            
            double uqq = (3.0*Q2)/(4.0*r2*r2*Math.sqrt(r2));
            double uQuad = (uqq)*(1.0-5.0*(cos1*cos1+cos2*cos2+3*cos1*cos1*cos2*cos2) +2*(temp*temp));
            ener += uQuad;
        }
        return ener;
    }

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }

    public double getEpsilon() {return epsilon;}

    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4*epsilon;
        epsilon48 = 48*epsilon;
    }

    public void setQuadrupolarMomentSquare(double moment){
        Q2=moment;
    }
    
    public double getQuadrupolarMomentSquare() {
        return Q2;
    }
    
    public Vector[] gradient(IAtomList pair, Tensor pressureTensor) {
        return gradient(pair);
    }

    public Vector[] gradient(IAtomList pair) {
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)pair.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        double rdudr = epsilon48*s6*(s6 - 0.5);
        gradient[0].Ea1Tv1(rdudr/r2,dr);

        if(Q2!=0.0){
            // normalize dr, the vector between the molecules
            drunit.E(dr);
            double r12 = Math.sqrt(r2);
            drunit.TE(1/r12);

            // v1 is the orientation of molecule 1
            Vector v1 = atom1.getOrientation().getDirection();

            // v2 is the orientation of molecule 2
            Vector v2 = atom2.getOrientation().getDirection();

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and dr
            double cos1 = v1.dot(drunit);
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and dr
            double cos2 = v2.dot(drunit);
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to dr
            // between the molecules multiplied by sin1 and sin2
            
            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double temp = v1.dot(v2) - 5*cos1*cos2;
            
            double uqq = (3.0*Q2)/(4.0*r2*r2*r12);
            double uQuad = 5/r12*(uqq)*(1.0-5.0*(cos1*cos1+cos2*cos2+3*cos1*cos1*cos2*cos2) +2*(temp*temp));
            gradient[0].PEa1Tv1(uQuad, drunit);

            // calculate dcos1 / dr1  (gradient of cos1)
            dcos1dr.Ea1Tv1(-r2, v1);
            dcos1dr.PEa1Tv1(v1.dot(dr), dr);
            dcos1dr.TE(-1/(r12*r2));

            // calculate dcos2 / dr1  (gradient of cos2)
            dcos2dr.Ea1Tv1(-r2, v2);
            dcos2dr.PEa1Tv1(v2.dot(dr), dr);
            dcos2dr.TE(-1/(r12*r2));

            // now actually add terms to the gradients
            gradient[0].PEa1Tv1(uqq*10*cos1, dcos1dr);
            gradient[0].PEa1Tv1(uqq*10*cos2, dcos2dr);
            gradient[0].PEa1Tv1(uqq*20*v1.dot(v2)*cos2, dcos1dr);
            gradient[0].PEa1Tv1(uqq*20*v1.dot(v2)*cos1, dcos2dr);
            gradient[0].PEa1Tv1(-uqq*70*cos1*cos2*cos2, dcos1dr);
            gradient[0].PEa1Tv1(-uqq*70*cos1*cos1*cos2, dcos2dr);
        }

        // pairwise additive, so
        gradient[1].Ea1Tv1(-1,gradient[0]);
        
        return gradient;
    }

    public double virial(IAtomList atoms) {
        gradient(atoms);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)atoms.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        boundary.nearestImage(dr);
        double v = gradient[1].dot(dr);
        if (Double.isInfinite(v)) {
            throw new RuntimeException("oops "+v);
        }
        return gradient[1].dot(dr);
    }

    public double hyperVirial(IAtomList pair) {
        throw new MethodNotImplementedException();
    }

    /**
     * Sets the temperature used for Boltzmann-weighting of the orientational
     * average energy used in u(double) and integral(double)
     */
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    // return a Boltzmann-weighted orientational average of the energy
    // for the given distance
    public double u(double r2) {
        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        double r4 = r2*r2;
        return epsilon4*s6*(s6 - 1.0) - (14.0/(5.0*temperature))*Q2*Q2/(r4*r4*r2);
    }
    
    public double du(double r2) {
        // should return a Boltzmann-weighted orientational average of du/dr
        throw new RuntimeException("please implement me");
    }
    
    public double integral(double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        double rc = sigma/rC;
        double sigmaD = space.powerD(sigma);
        double rc3 = rc*rc*rc;
        double rc6 = rc3*rc3;
        double rc12 = rc6*rc6;
        // LJ part
        double uInt = 4.0*epsilon*sigmaD*A*(rc12/9 - rc6/3)/rc3;  //complete LRC is obtained by multiplying by N1*N2/V
        double r2 = rC;
        double r4 = r2*r2;
        return uInt - (7.0/(25.0*temperature))*Q2*Q2/(r4*r4*r2*rC);
    }

    private static final long serialVersionUID = 1L;
    private double sigma , sigma2;
    private double epsilon, epsilon4, epsilon48;
    private double hsdiasq=1.0/Math.sqrt(2);
    private double Q2;
    private Boundary boundary;
    private final Vector dr, drunit, dcos1dr, dcos2dr;
    private final Vector[] gradient;
    protected double temperature;
}
