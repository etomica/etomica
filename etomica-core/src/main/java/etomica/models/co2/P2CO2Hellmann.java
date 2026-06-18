/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.co2;

import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.potential.IPotential2;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.CompoundUnit;
import etomica.units.Kelvin;
import etomica.units.Unit;
import etomica.units.dimensions.Null;
import etomica.util.Constants;

/**
 * Ab initio potential for CO2-CO2 developed by R. Hellmann.
 *
 * http://dx.doi.org/10.1016/j.cplett.2014.08.057
 * 
 * @author Andrew Schultz
 */
public class P2CO2Hellmann implements IPotential2 {

    public static final double[] posA = new double[]{-1.28815171291, -1.17769797231, -0.18133162098, 0.00000000000, 0.18133162098, 1.17769797231, 1.28815171291};
    public static final double[] posB = new double[]{-1.28741781626, -1.18192825424, -0.18607849166, 0.00000000000, 0.18607849166, 1.18192825424, 1.28741781626};
    protected static final double[] qA = new double[]{-191.214602050, 163.370217205, -2747.80131789, 5551.29140547, -2747.80131789, 163.370217205, -191.214602050};
    protected static final double[] qB = new double[]{-197.417207828, 168.070083318, -2559.64083227, 5177.97591356, -2559.64083227, 168.070083318, -197.417207828};
    protected static final int[] siteID = new int[]{0,1,2,3,2,1,0};
    protected static final double[] AA = new double[]{-0.206097356066E+07, 0.636857652220E+07, -0.194431965659E+08, 0.379428245681E+08, -0.112075555946E+07, 0.460535048972E+08, -0.118847782748E+09, -0.130239550952E+08, 0.425605973419E+08, -0.114999594389E+09};
    protected static final double[] AB = new double[]{-0.247910365353E+07, 0.659160470472E+07, -0.197776308389E+08, 0.384165630648E+08, -0.124570324466E+07, 0.451317323034E+08, -0.116048612008E+09, -0.103079402689E+08, 0.340824968085E+08, -0.915027698701E+08};
    protected static final double[] alphaA = new double[]{2.15463968686, 3.20330234880, 2.44546409876, 2.47207342773, 1.67720522533, 2.65642378772, 2.76521425413, 2.99122464510, 2.82606075114, 2.94148736356};
    protected static final double[] alphaB = new double[]{2.08319218048, 3.16681447768, 2.46163539534, 2.48589087370, 1.67813668662, 2.65969570294, 2.77169644514, 2.98535796569, 2.75870881239, 2.87267355769};
    protected static final double[] bA = new double[]{3.10392728450, 2.52351588006, 1.60598918542, 1.91489841285, 2.06812316859, 1.49177611494, 4.13165202187, 2.77767103339, 2.48867537323, 2.30146553200};
    protected static final double[] bB = new double[]{3.14980106637, 2.46903752251, 1.57103563097, 1.89845841233, 2.14451960163, 1.46843191121, 4.14021127755, 2.72634741238, 2.44815795987, 2.27614875317};
    protected static final double[] C6A = new double[]{-0.259827667988E+08, 0.606040735312E+08, -0.128063652046E+09, 0.106278254027E+09, -0.957062410858E+08, 0.729091325685E+08, 0.259781327401E+08, 0.146735186656E+09, -0.345229777876E+09, 0.709212843263E+09};
    protected static final double[] C6B = new double[]{-0.306747626563E+08, 0.698469835305E+08, -0.143806191593E+09, 0.121226824365E+09, -0.109398472925E+09, 0.811702881095E+08, 0.263241896284E+08, 0.126349448908E+09, -0.285769208067E+09, 0.551179708953E+09};
    protected static final double[] C8A = new double[]{0.189027654711E+09, -0.687384314073E+09, 0.312617758706E+10, -0.253453395293E+10, 0.930718903055E+09, -0.183474802396E+10, -0.175211033214E+09, -0.698068389813E+09, 0.191961111413E+10, -0.337571242220E+10};
    protected static final double[] C8B = new double[]{0.211522217149E+09, -0.810638994730E+09, 0.355929066714E+10, -0.286891373977E+10, 0.114677667224E+10, -0.210805303525E+10, -0.173569859005E+09, -0.496759975158E+09, 0.122323855871E+10, -0.131218053988E+10};
    protected static final int nsites = 7;
    
    protected final Vector ri, rj, drij, torque;
    protected final double[] pos, q;
    protected final double[][] A, alpha, b, C6, C8;

    protected final Space space;
    protected final Vector[][] gradientAndTorque;
    protected final static double mass;
    static {
        mass = Carbon.INSTANCE.getMass() + 2*Oxygen.INSTANCE.getMass();
    }
    
    public enum Parameters {
        A, B
    }
    
    public P2CO2Hellmann(Space space, Parameters param) {
        this.space = space;
        A = new double[4][4];
        alpha = new double[4][4];
        b = new double[4][4];
        C6 = new double[4][4];
        C8 = new double[4][4];
        q = new double[4];
        Unit sqrtK = new CompoundUnit(new Unit[]{Kelvin.UNIT}, new double[]{0.5});
        if (param == Parameters.A) {
            pos = posA;
            for (int i=0; i<q.length; i++) {
                q[i] = sqrtK.toSim(qA[i]);
            }
            ijInit(AA, A, Kelvin.UNIT);
            ijInit(alphaA, alpha, Null.UNIT);
            ijInit(bA, b, Null.UNIT);
            ijInit(C6A, C6, Kelvin.UNIT);
            ijInit(C8A, C8, Kelvin.UNIT);
        }
        else {
            pos = posB;
            for (int i=0; i<q.length; i++) {
                q[i] = sqrtK.toSim(qB[i]);
            }
            ijInit(AB, A, Kelvin.UNIT);
            ijInit(alphaB, alpha, Null.UNIT);
            ijInit(bB, b, Null.UNIT);
            ijInit(C6B, C6, Kelvin.UNIT);
            ijInit(C8B, C8, Kelvin.UNIT);
        }
        ri = space.makeVector();
        rj = space.makeVector();
        drij = space.makeVector();
        torque = space.makeVector();
        gradientAndTorque = new Vector[2][2];
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                gradientAndTorque[i][j] = space.makeVector();
            }
        }
    }

    public double getPos(int idx) {
        return pos[idx];
    }

    public double getQ(int idx) {
        return q[siteID[idx]];
    }

    protected void ijInit(double[] x, double[][] xx, Unit unit) {
        int k = 0;
        for (int i=0; i<4; i++) {
            xx[i][i] = unit.toSim(x[k]);
            k++;
            for (int j=i+1; j<4; j++) {
                xx[i][j] = xx[j][i] = unit.toSim(x[k]);
                k++;
            }
        }
    }

    protected boolean debug = false;

    public double uduTorque(Vector dr12, IAtom a1, IAtom a2, Vector f1, Vector f2, Vector t1, Vector t2) {

        IAtomOriented atom1 = (IAtomOriented)a1;
        IAtomOriented atom2 = (IAtomOriented)a2;
        Vector cm0 = atom1.getPosition();
        Vector cm1 = atom2.getPosition();
        Vector or0 = atom1.getOrientation().getDirection();
        Vector or1 = atom2.getOrientation().getDirection();
        double u = 0;
        for (int i=0; i<7; i++) {
            int ii = siteID[i];
            ri.E(cm0);
            ri.PEa1Tv1(pos[i], or0);
            for (int j=0; j<7; j++) {
                int jj = siteID[j];
                rj.E(cm1);
                rj.PEa1Tv1(pos[j], or1);
                drij.Ev1Mv2(rj, ri);
                double rij2 = drij.squared();
                double rij = Math.sqrt(rij2);
                if (rij < 0.9) {
                    f1.E(Double.NaN);
                    f2.E(Double.NaN);
                    return Double.POSITIVE_INFINITY;
                }
                double ar = alpha[ii][jj]*rij;
                double rduExpdr = -A[ii][jj]*ar*Math.exp(-ar);

                double sum = 1;
                double dsum = 0;
                double br = b[ii][jj]*rij;
                double term = 1;
                double dterm = 1;
                for (int k=1; k<=6; k++) {
                    term *= br/k;
                    if (k==1) {
                        dterm = br;
                    }
                    else {
                        dterm *= br/(k-1);
                    }
                    sum += term;
                    dsum += dterm;
                }
                if (sum==1) {
                    f1.E(Double.NaN);
                    f2.E(Double.NaN);
                    return Double.POSITIVE_INFINITY;
                }
                double uExp = A[ii][jj]*Math.exp(-alpha[ii][jj]*rij);
                double rij6 = rij2*rij2*rij2;
                double expbr = Math.exp(-br);
                double u6 = -(1-expbr*sum)*C6[ii][jj]/rij6;
                double rdu6dr = (-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum))*C6[ii][jj]/rij6;

                for (int k=7; k<=8; k++) {
                    term *= br/k;
                    dterm *= br/(k-1);
                    sum += term;
                    dsum += dterm;
                }
                double rij8 = rij6*rij2;
                double u8 = -(1-expbr*sum)*C8[ii][jj]/rij8;
                double rdu8dr = (-expbr*br*sum + expbr*dsum + 8*(1-expbr*sum))*C8[ii][jj]/rij8;
//                rdu8dr = (expbr*dsum);

                double uCharge = q[ii]*q[jj]/rij;
                double rduChargedr = -q[ii]*q[jj]/rij;

                u += uExp + u6 + u8 + uCharge;
                if (Double.isNaN(u)) throw new RuntimeException("oops");
                double rdudr = rduExpdr + rdu6dr + rdu8dr + rduChargedr;
                // we're done with drij.  replace it with gradient on 1 (force on 0)
//                System.out.println("g drij "+drij+" "+rdudr/rij2);
                drij.TE(rdudr/rij2);
                f1.PE(drij);
//                System.out.println("grad 1 "+gradientAndTorque[0][1]);
                f2.ME(drij);
                torque.Ea1Tv1(pos[j], or1);
                torque.XE(drij);
                t2.ME(torque);
                torque.Ea1Tv1(pos[i], or0);
                torque.XE(drij);
                t1.PE(torque);
            }
        }
        return u;
    }

    public double u(Vector dr12, IAtom a1, IAtom a2) {
        IAtomOriented atom1 = (IAtomOriented)a1;
        IAtomOriented atom2 = (IAtomOriented)a2;
        Vector cm0 = atom1.getPosition();
        Vector cm1 = atom2.getPosition();
        Vector or0 = atom1.getOrientation().getDirection();
        Vector or1 = atom2.getOrientation().getDirection();
        double u = 0;
        boolean checkme = false;
        for (int i=0; i<7; i++) {
            int ii = siteID[i];
            ri.E(cm0);
            ri.PEa1Tv1(pos[i], or0);
            for (int j=0; j<7; j++) {
                int jj = siteID[j];
                rj.E(cm1);
                rj.PEa1Tv1(pos[j], or1);
                double rij2 = rj.Mv1Squared(ri);
                double rij = Math.sqrt(rij2);
//                System.out.println("rij "+rj+" "+ri+" "+rij);
                if (rij < 0.9) return Double.POSITIVE_INFINITY;
                if (rij < 1.2) checkme = true;
                double uExp = A[ii][jj]*Math.exp(-alpha[ii][jj]*rij);

                double sum = 1;
                double br = b[ii][jj]*rij;
                double term = 1;
                for (int k=1; k<=6; k++) {
                    term *= br/k;
                    sum += term;
                }
                if (sum==1) return Double.POSITIVE_INFINITY;
                double rij6 = rij2*rij2*rij2;
                double expbr = Math.exp(-br);
                double u6 = -(1-expbr*sum)*C6[ii][jj]/rij6;

                for (int k=7; k<=8; k++) {
                    term *= br/k;
                    sum += term;
                }
                double rij8 = rij6*rij2;
                double u8 = -(1-expbr*sum)*C8[ii][jj]/rij8;

                double uCharge = q[ii]*q[jj]/rij;
                u += uExp + u6 + u8 + uCharge;
                if (debug) {
                    System.out.println(rij+" "+(uExp+u6+u8+uCharge));
                }
                if (Double.isNaN(u)) throw new RuntimeException("oops");
            }
        }
        if (debug) return u;
        if (u<-1000 || Double.isNaN(u)) {
            System.out.println(u);
            debug = true;
            u(dr12, a1, a2);
            throw new RuntimeException("oops, too much");
        }
        if (checkme && u<10000) {
            System.out.println(u);
        }
        return u;
    }

    public P2CO2SC makeSemiclassical(double temperature) {
        return new P2CO2SC(temperature);
    }

    public class P2CO2SC implements IPotential2 {

        protected final Vector[][] gi;
        protected final Tensor tt0Tensor, tt1Tensor, rr0Tensor, rr1Tensor;
        protected final Tensor ijTensor, rTensor0, rTensor1, identity;
        protected final Tensor ijRTensor;
        protected final Tensor rot0, rot1;
        protected final Vector or01, or11, or02, or12;
        protected final Vector[] allOr0, allOr1;
        protected final Vector drijRot;
        protected final double moment;
        public double[][] d2tot = new double[2][6];
        protected final double temperature, fac;
        protected final double bondL = 1.1625;

        public P2CO2SC(double temperature) {
            ijTensor = space.makeTensor();
            identity = space.makeTensor();
            tt0Tensor = space.makeTensor();
            tt1Tensor = space.makeTensor();
            rr0Tensor = space.makeTensor();
            rr1Tensor = space.makeTensor();
            rTensor0 = space.makeTensor();
            rTensor1 = space.makeTensor();
            ijRTensor = space.makeTensor();
            identity.E(new double[][]{{1,0,0},{0,1,0},{0,0,1}});
            gi = new Vector[2][7];
            for (int i=0; i<7; i++) {
                gi[0][i] = space.makeVector();
                gi[1][i] = space.makeVector();
            }
            or01 = space.makeVector();
            or11 = space.makeVector();
            or02 = space.makeVector();
            or12 = space.makeVector();
            allOr0 = new Vector[]{null, or01, or02};
            allOr1 = new Vector[]{null, or11, or12};
            drijRot = space.makeVector();
            rot0 = space.makeTensor();
            rot1 = space.makeTensor();
            moment = 2*Oxygen.INSTANCE.getMass()*bondL*bondL;
            
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24/2)/temperature;
        }

        protected void getPerp(Vector or, Vector perp1, Vector perp2) {
            int max = 0;
            if (Math.abs(or.getX(1)) > Math.abs(or.getX(0))) max=1;
            if (Math.abs(or.getX(2)) > Math.abs(or.getX(max))) max=2;
            int min = 0;
            if (Math.abs(or.getX(1)) < Math.abs(or.getX(0))) min=1;
            if (Math.abs(or.getX(2)) < Math.abs(or.getX(min))) min=2;
            perp1.E(0);
            perp1.setX(min, or.getX(max));
            perp1.setX(max, -or.getX(min));
            perp1.normalize();
            perp1.PEa1Tv1(-perp1.dot(or), or);
            perp1.normalize();
            perp2.E(or);
            perp2.XE(perp1);
        }

        public double u(Vector dr12, IAtom a1, IAtom a2) {
            IAtomOriented atom1 = (IAtomOriented) a1;
            IAtomOriented atom2 = (IAtomOriented) a2;
            Vector cm0 = atom1.getPosition();
            Vector cm1 = atom2.getPosition();
            Vector or0 = atom1.getOrientation().getDirection();
            getPerp(or0, or01, or02);
            allOr0[0] = or0;
            rot0.E(allOr0);
            rot0.invert();

            Vector or1 = atom2.getOrientation().getDirection();
            getPerp(or1, or11, or12);
            allOr1[0] = or1;
            rot1.E(allOr1);
            rot1.invert();

            tt0Tensor.E(0);
            tt1Tensor.E(0);
            rr0Tensor.E(0);
            rr1Tensor.E(0);
            for (int i=0; i<7; i++) {
                gi[0][i].E(0);
                gi[1][i].E(0);
            }
            for (int i=0; i<6; i++ ){
                d2tot[0][i] = 0;
                d2tot[1][i] = 0;
            }
            double u = 0;
            for (int i=0; i<7; i++) {
                int ii = siteID[i];
                ri.Ea1Tv1(pos[i], or0);
//                rTensor0.setComponent(0, 1, 0);
//                rTensor0.setComponent(1, 0, -0);
//                rTensor0.setComponent(0, 2, -0);
//                rTensor0.setComponent(2, 0, 0);
                rTensor0.setComponent(1, 2, pos[i]);
                rTensor0.setComponent(2, 1, -pos[i]);
                ri.PE(cm0);
                for (int j=0; j<7; j++) {
                    int jj = siteID[j];
                    rj.Ea1Tv1(pos[j], or1);
//                    rTensor1.setComponent(0, 1, 0);
//                    rTensor1.setComponent(1, 0, -0);
//                    rTensor1.setComponent(0, 2, -0);
//                    rTensor1.setComponent(2, 0, 0);
                    rTensor1.setComponent(1, 2, pos[j]);
                    rTensor1.setComponent(2, 1, -pos[j]);
                    rj.PE(cm1);
                    drij.Ev1Mv2(rj, ri);
                    double rij2 = drij.squared();
                    double rij = Math.sqrt(rij2);
                    if (rij < 0.9) {
                        return Double.POSITIVE_INFINITY;
                    }
                    double ar = alpha[ii][jj]*rij;
                    double uExp = A[ii][jj]*Math.exp(-alpha[ii][jj]*rij);
                    double rduExpdr = -A[ii][jj]*ar*Math.exp(-ar);
                    double r2du2Expdr2 = A[ii][jj]*ar*ar*Math.exp(-ar);

                    double sum = 1;
                    double dsum = 0;
                    double d2sum = 0;
                    double br = b[ii][jj]*rij;
                    double term = 1;
                    double dterm = 1;
                    double d2term = 0;
                    for (int k=1; k<=6; k++) {
                        term *= br/k;
                        if (k==1) {
                            dterm = br;
                        }
                        else {
                            dterm *= br/(k-1);
                        }
                        if (k==2) {
                            d2term = br*br;
                        }
                        else if (k>2) {
                            d2term *= br/(k-2);
                        }
                        sum += term;
                        dsum += dterm;
                        d2sum += d2term;
                    }
                    if (sum==1) {
                        return Double.POSITIVE_INFINITY;
                    }
                    double rij6 = rij2*rij2*rij2;
                    double expbr = Math.exp(-br);
                    double u6 = -(1-expbr*sum)*C6[ii][jj]/rij6;
                    double rdu6dr = (-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum))*C6[ii][jj]/rij6;
                    double r2du26dr2 = (expbr*br*br*sum-expbr*br*dsum -br*expbr*dsum+expbr*d2sum + 6*(-1 + br*expbr*sum-expbr*dsum+expbr*sum) + -6*(-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum)))*C6[ii][jj]/rij6;

                    for (int k=7; k<=8; k++) {
                        term *= br/k;
                        dterm *= br/(k-1);
                        d2term *= br/(k-2);
                        sum += term;
                        dsum += dterm;
                        d2sum += d2term;
                    }
                    double rij8 = rij6*rij2;
                    double u8 = -(1-expbr*sum)*C8[ii][jj]/rij8;
                    double rdu8dr = (-expbr*br*sum + expbr*dsum + 8*(1-expbr*sum))*C8[ii][jj]/rij8;
//                    rdu8dr = (expbr*dsum);
                    double r2du28dr2 = (expbr*br*br*sum-expbr*br*dsum -br*expbr*dsum+expbr*d2sum + 8*(-1 + br*expbr*sum-expbr*dsum+expbr*sum) + -8*(-expbr*br*sum + expbr*dsum + 8*(1-expbr*sum)))*C8[ii][jj]/rij8;
//                    r2du28dr2 = (-br*expbr*dsum+expbr*d2sum);

                    double uCharge = q[ii]*q[jj]/rij;
                    double rduChargedr = -q[ii]*q[jj]/rij;
                    double r2d2uChargedr2 = 2*q[ii]*q[jj]/rij;

                    u += uExp + u6 + u8 + uCharge;
                    double rdudr = rduExpdr + rdu6dr + rdu8dr + rduChargedr;
                    double r2d2udr2 = r2du2Expdr2 + r2du26dr2 + r2du28dr2 + r2d2uChargedr2;
//                    System.out.println(rij+" "+foo8/rij+" "+foo82/rij2);
//                    if (i==0) System.out.println(i+" "+j+" "+rij+" "+u+" "+rdudr/rij+" "+r2d2udr2/rij2);

                    // molecule 0
                    drijRot.E(drij);
                    rot0.transform(drijRot);
                    ijTensor.Ev1v2(drijRot, drijRot);
                    ijTensor.TE((rdudr - r2d2udr2)/(rij2*rij2));
                    ijTensor.PEa1Tt1(-rdudr/rij2, identity);

                    tt0Tensor.ME(ijTensor);

                    rTensor0.transpose();
                    ijRTensor.E(rTensor0);
                    ijRTensor.TE(ijTensor);
                    rTensor0.transpose();
                    ijRTensor.TE(rTensor0);
                    rr0Tensor.ME(ijRTensor);

                    drijRot.TE(rdudr/rij2);
                    gi[0][i].ME(drijRot);


                    // molecule 1
                    drijRot.E(drij);
                    rot1.transform(drijRot);
                    ijTensor.Ev1v2(drijRot, drijRot);
//                    System.out.println("r2 "+rij2+" "+ijTensor.component(0,0));
                    ijTensor.TE((rdudr - r2d2udr2)/(rij2*rij2));
                    ijTensor.PEa1Tt1(-rdudr/rij2, identity);

                    // we only need tt for 0, both have the same contribution
                    tt1Tensor.ME(ijTensor);
//                    System.out.println("d2 "+r2d2udr2/rij2+" "+ijTensor.component(0,0)+" "+tt1Tensor.component(0,0));

                    rTensor1.transpose();
                    ijRTensor.E(rTensor1);
                    ijRTensor.TE(ijTensor);
                    rTensor1.transpose();
                    ijRTensor.TE(rTensor1);
                    rr1Tensor.ME(ijRTensor);

                    drijRot.TE(rdudr/rij2);
                    gi[1][j].PE(drijRot);
                }
            }

            for (int i=0; i<7; i++) {
                ri.E(0);
                ri.setX(0, pos[i]);
                rr0Tensor.PEv1v2(ri, gi[0][i]);
                rr0Tensor.PEa1Tt1(-pos[i]*gi[0][i].getX(0), identity);

                rj.E(0);
                rj.setX(0, pos[i]);
                rr1Tensor.PEv1v2(rj, gi[1][i]);
                rr1Tensor.PEa1Tt1(-pos[i]*gi[1][i].getX(0), identity);

            }
            double sum = 0;
            for (int i=0; i<3; i++){
                d2tot[0][i] += tt0Tensor.component(i,i);
                d2tot[1][i] += tt1Tensor.component(i,i);
                sum += tt0Tensor.component(i,i)/mass;
            }
            for (int i=1; i<3; i++){
                d2tot[0][3+i] += rr0Tensor.component(i,i);
                d2tot[1][3+i] += rr1Tensor.component(i,i);
                sum += (rr0Tensor.component(i,i) + rr1Tensor.component(i,i))/(2*moment);
            }
            return u + fac*sum;
        }
    }
}
