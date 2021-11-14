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
public class P2HSDipoleAtomic implements Potential2Soft {

    public P2HSDipoleAtomic(Space space, double sigma, double dipole, double rCut) {
        this.space = space;
        setSigma(sigma);
        a = new Vector[3];
        a[0] = space.makeVector();
        a[1] = space.makeVector();
        a[2] = space.makeVector();

        dr = space.makeVector();
        drunit = space.makeVector();
        runit = space.makeVector();
        work = space.makeVector();

        setDipole(dipole);
        this.rCut = rCut;
        hessian = new Hessian(space.makeTensor(), space.makeTensor(), space.makeTensor(),
                space.makeTensor(), space.makeTensor(), space.makeTensor());
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
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

        if(r2 > rCut*rCut){
            return 0;
        }

        if(r2 < sigma2) { // hard core
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
        double cos_D1_r = v1.dot(dr12)/r;
        //cos(r12 and dipole 2)
        double cos_r_D2=dr12.dot(v2)/r;
        double r12Magnitude = Math.sqrt(r2);
        double ener = dipole * dipole *  (cos_D1_D2 - 3.0  * cos_D1_r* cos_r_D2);
        ener = ener/r12Magnitude /r2 ;
        return ener;
    }

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s * s;
    }

    public void setDipole(double moment){
        dipole=moment;
    }

    @Override
    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        double r2 = dr12.squared();
        if (r2 > rCut*rCut) {
            return 0;
        }
        if (r2 < sigma2) {
            return Double.POSITIVE_INFINITY;
        }

        double momentSq = dipole*dipole;
        double u = 0;
        if(momentSq!=0.0){
            double s2 = 1/r2;
            double s1 = Math.sqrt(s2);
            // normalize dr, the vector between the molecules
            drunit.E(dr12);
            drunit.TE(s1);

            Vector v1 = ((IAtomOriented)atom1).getOrientation().getDirection();

            Vector v2 = ((IAtomOriented)atom2).getOrientation().getDirection();

            double cos_D1_D2 = v1.dot(v2);
            double cos_D1_r = v1.dot(drunit);
            double cos_D2_r = v2.dot(drunit);

            double fac = momentSq * (s2*s1);
            double dfac = momentSq * (3*s2*s2*s1);
            double udd = cos_D1_D2 - 3.0*cos_D1_r*cos_D2_r;
            work.Ea1Tv1(cos_D2_r, v1);
            work.PEa1Tv1(cos_D1_r, v2);
            work.PEa1Tv1(-2*cos_D1_r*cos_D2_r*s1, dr12);
            work.TE(3.0*s1*fac);
            f1.ME(work);
            f1.PEa1Tv1(-dfac*udd, dr12);
            f2.PE(work);
            f2.PEa1Tv1(+dfac*udd, dr12);

            work.E(v1);
            work.XE(v2);
            Vector t = Vector.d(t1.getD());
            t.E(v1);
            t.XE(drunit);
            t.TE(3.0*cos_D2_r);
            t.ME(work);
            t.TE(fac);
            t1.PE(t);

            t.E(v2);
            t.XE(drunit);
            t.TE(3.0*cos_D1_r);
            t.PE(work);
            t.TE(fac);
            t2.PE(t);

            double r12Magnitude = Math.sqrt(r2);
            u = dipole * dipole *  (cos_D1_D2 - 3.0  * cos_D1_r* cos_D2_r);
            u = u/r12Magnitude /r2 ;
        }


        // pairwise additive, so
        return u;
    }

    public Hessian d2u(Vector dr12, IAtom a1, IAtom a2){
        IAtomOriented atom1 = (IAtomOriented) a1;
        IAtomOriented atom2 = (IAtomOriented) a2;
        Vector ei =  atom1.getOrientation().getDirection();
        Vector ej =  atom2.getOrientation().getDirection();

        double exi = ei.getX(0);//ei and ej is the dipole orientation
        double eyi = ei.getX(1);
        double ezi = ei.getX(2);
        double exj = ej.getX(0);
        double eyj = ej.getX(1);
        double ezj = ej.getX(2);

        Vector deidxi = space.makeVector();
        Vector deidyi = space.makeVector();
        Vector deidzi = space.makeVector();
        Vector dejdxj = space.makeVector();
        Vector dejdyj = space.makeVector();
        Vector dejdzj = space.makeVector();

        double [] dejdxjD = {0,-ezj,eyj};
        double [] dejdyjD = {ezj,0,-exj};
        double [] dejdzjD = {-eyj,exj,0};
        dejdxj.E(dejdxjD);
        dejdyj.E(dejdyjD);
        dejdzj.E(dejdzjD);

        double r2 = dr12.squared();
        hessian.o1o1.E(0);
        hessian.o1o2.E(0);
        hessian.o2o2.E(0);
        if (r2 > rCut*rCut) {
            return hessian;
        }
        double r = Math.sqrt(r2);
        runit.Ea1Tv1(1.0/r, dr12);
        double coeff = dipole*dipole/r2/r;

        //ei cross dejdxj - 3.0*(dejdxj.runit)*(ei cross runite)
        dr.E(ei);
        dr.XE(dejdxj);
        a[0].E(ei);
        a[0].XE(runit);
        a[0].TE(-3.0*dejdxj.dot(runit));
        a[0].PE(dr);

        dr.E(ei);
        dr.XE(dejdyj);
        a[1].E(ei);
        a[1].XE(runit);
        a[1].TE(-3.0*dejdyj.dot(runit));
        a[1].PE(dr);

        dr.E(ei);
        dr.XE(dejdzj);
        a[2].E(ei);
        a[2].XE(runit);
        a[2].TE(-3.0*dejdzj.dot(runit));
        a[2].PE(dr);

        hessian.o1o2.E(a);
        hessian.o1o2.TE(coeff);//ij

        double [] deidxiD = {0,-ezi,eyi};
        double [] deidyiD = {ezi,0,-exi};
        double [] deidziD = {-eyi,exi,0};
        deidxi.E(deidxiD);
        deidyi.E(deidyiD);
        deidzi.E(deidziD);

//		deidxi cross ej - 3.0*(deidxi cross runite)*(ej.runit)
        dr.E(deidxi);
        dr.XE(ej);
        a[0].E(deidxi);
        a[0].XE(runit);
        a[0].TE(-3.0*ej.dot(runit));
        a[0].PE(dr);

        dr.E(deidyi);
        dr.XE(ej);
        a[1].E(deidyi);
        a[1].XE(runit);
        a[1].TE(-3.0*ej.dot(runit));
        a[1].PE(dr);

        dr.E(deidzi);
        dr.XE(ej);
        a[2].E(deidzi);
        a[2].XE(runit);
        a[2].TE(-3.0*ej.dot(runit));
        a[2].PE(dr);


        hessian.o1o1.E(a);
        hessian.o1o1.TE(coeff);//ii

        dr.E(dejdxj);
        dr.XE(ei);
        a[0].E(dejdxj);
        a[0].XE(runit);
        a[0].TE(-3.0*ei.dot(runit));
        a[0].PE(dr);

        dr.E(dejdyj);
        dr.XE(ei);
        a[1].E(dejdyj);
        a[1].XE(runit);
        a[1].TE(-3.0*ei.dot(runit));
        a[1].PE(dr);

        dr.E(dejdzj);
        dr.XE(ei);
        a[2].E(dejdzj);
        a[2].XE(runit);
        a[2].TE(-3.0*ei.dot(runit));
        a[2].PE(dr);

        hessian.o2o2.E(a);
        hessian.o2o2.TE(coeff);//jj

        return hessian;

    }


    private final Space space;
    private double sigma , sigma2;
    private final Vector dr,drunit,work,runit;
    private double dipole;
    private double rCut;
    protected final Vector[] a;
    protected final Hessian hessian;

    @Override
    public double energy(IAtomList atoms) {
        return 0;
    }

    @Override
    public double u(double r2) {
        return 0;
    }
}
