/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.co2;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Oxygen;
import etomica.potential.IPotentialTorque;
import etomica.potential.P2SemiclassicalAtomic;
import etomica.potential.P2SemiclassicalAtomic.AtomInfo;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.IOrientation3D;
import etomica.space3d.Orientation3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.CompoundUnit;
import etomica.units.Kelvin;
import etomica.units.Null;
import etomica.units.Unit;
import etomica.util.Constants;
import etomica.util.RandomMersenneTwister;

/**
 * Ab initio potential for CO2-CO2 developed by R. Hellmann.
 *
 * http://dx.doi.org/10.1016/j.cplett.2014.08.057
 * 
 * @author Andrew Schultz
 */
public class P2CO2Hellmann implements IPotentialTorque {

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

    public double virial(IAtomList atoms) {
        return 0;
    }

    public Vector[] gradient(IAtomList atoms) {
        return gradientAndTorque(atoms)[0];
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradientAndTorque(atoms)[0];
    }
    protected boolean debug = false;

    public Vector[][] gradientAndTorque(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
        Vector cm0 = atom0.getPosition();
        Vector cm1 = atom1.getPosition();
        Vector or0 = atom0.getOrientation().getDirection();
        Vector or1 = atom1.getOrientation().getDirection();
        gradientAndTorque[0][0].E(0);
        gradientAndTorque[0][1].E(0);
        gradientAndTorque[1][0].E(0);
        gradientAndTorque[1][1].E(0);
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
                    gradientAndTorque[0][0].E(Double.NaN);
                    gradientAndTorque[1][0].E(Double.NaN);
                    gradientAndTorque[0][1].E(Double.NaN);
                    gradientAndTorque[1][1].E(Double.NaN);
                    return gradientAndTorque;
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
                    gradientAndTorque[0][0].E(Double.NaN);
                    gradientAndTorque[1][0].E(Double.NaN);
                    gradientAndTorque[0][1].E(Double.NaN);
                    gradientAndTorque[1][1].E(Double.NaN);
                    return gradientAndTorque;
                }
                double rij6 = rij2*rij2*rij2;
                double expbr = Math.exp(-br);
                double rdu6dr = (-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum))*C6[ii][jj]/rij6; 

                for (int k=7; k<=8; k++) {
                    term *= br/k;
                    dterm *= br/(k-1);
                    sum += term;
                    dsum += dterm;
                }
                double rij8 = rij6*rij2;
                double rdu8dr = (-expbr*br*sum + expbr*dsum + 8*(1-expbr*sum))*C8[ii][jj]/rij8;
//                rdu8dr = (expbr*dsum);

                double rduChargedr = -q[ii]*q[jj]/rij;
                
                double rdudr = rduExpdr + rdu6dr + rdu8dr + rduChargedr;
                // we're done with drij.  replace it with gradient on 1 (force on 0)
//                System.out.println("g drij "+drij+" "+rdudr/rij2);
                drij.TE(rdudr/rij2);
                gradientAndTorque[0][1].PE(drij);
//                System.out.println("grad 1 "+gradientAndTorque[0][1]);
                gradientAndTorque[0][0].ME(drij);
                torque.Ea1Tv1(pos[j], or1);
                torque.XE(drij);
                gradientAndTorque[1][1].ME(torque);
                torque.Ea1Tv1(pos[i], or0);
                torque.XE(drij);
                gradientAndTorque[1][0].PE(torque);
            }
        }
        return gradientAndTorque;
    }
    
    public double energy(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
        Vector cm0 = atom0.getPosition();
        Vector cm1 = atom1.getPosition();
        Vector or0 = atom0.getOrientation().getDirection();
        Vector or1 = atom1.getOrientation().getDirection();
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
            energy(atoms);
            throw new RuntimeException("oops, too much");
        }
        if (checkme && u<10000) {
            System.out.println(u);
        }
        return u;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        
    }

    public int nBody() {
        return 2;
    }
    
    public P2CO2SC makeSemiclassical(double temperature) {
        return new P2CO2SC(temperature);
    }

    public class P2CO2SC implements IPotentialAtomic {

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
        
        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public void setBox(Box box) {
            
        }

        public int nBody() {
            return 2;
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
        
        public double energy(IAtomList atoms) {
            IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
            IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
            Vector cm0 = atom0.getPosition();
            Vector cm1 = atom1.getPosition();
            Vector or0 = atom0.getOrientation().getDirection();
            getPerp(or0, or01, or02);
            allOr0[0] = or0;
            rot0.E(allOr0);
            rot0.invert();
            
            Vector or1 = atom1.getOrientation().getDirection();
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

    public static void main1(String[] args) {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(200);
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating species = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(species);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(species, 2);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
//        ((IAtomOriented)pair.getAtom(0)).getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
//        atom1.getOrientation().setDirection(space.makeVector(new double[]{1/Math.sqrt(2),0.5,0.5}));
        P2CO2Hellmann p2 = new P2CO2Hellmann(space, Parameters.B);
        P2SemiclassicalAtomic p2TI = new P2SemiclassicalAtomic(space, p2, temperature);
        final Vector[] rv = new Vector[4];
        for (int i=0; i<4; i++) {
            rv[i] = space.makeVector();
        }
        double om = Oxygen.INSTANCE.getMass();
        double bondLength = p2.pos[5] - p2.pos[1];
        rv[0].setX(0, om*bondLength*bondLength*0.25);
        rv[0].setX(1, om*bondLength*bondLength*0.25);
        p2TI.setAtomInfo(species.getLeafType(), new AtomInfo() {
            public Vector[] getMomentAndAxes(IAtomOriented molecule) {
                // rv[0,2] = 0
                // rv[3] is the orientation
                rv[3].E(molecule.getOrientation().getDirection());
                // rv[1] is an axis perpendicular to rv[3]
                rv[1].E(0);
                if (Math.abs(rv[3].getX(0)) < 0.5) {
                    rv[1].setX(0, 1);
                }
                else if (Math.abs(rv[3].getX(1)) < 0.5) {
                    rv[1].setX(1, 1);
                }
                else {
                    rv[1].setX(2, 1);
                }
                rv[2].Ea1Tv1(rv[1].dot(rv[3]), rv[3]);
                rv[1].ME(rv[2]);
                // rv[2] is an axis perpendicular to rv[3] and rv[1]
                rv[2].E(rv[1]);
                rv[2].XE(rv[3]);
                return rv;
            }
        });
        p2TI.setTemperature(temperature);
        System.out.println("or: "+((IAtomOriented)pair.getAtom(0)).getOrientation().getDirection()+" "+atom1.getOrientation().getDirection());
        Vector y = space.makeVector(new double[]{0.0,1.0,0.0});
        Vector z = space.makeVector(new double[]{0.0,0.0,1.0});
        double lg = 0;
        double lu = 0;
        atom1.getPosition().setX(0, 5);
        double dx = 0.1;
        P2CO2SC p2SC = p2.makeSemiclassical(temperature);

        for (int i=0; i<=4; i++) {
//            atom1.getPosition().setX(0, 5+i*dx);
            double u = p2.energy(pair);
            double uti = p2TI.energy(pair);
            double usc = p2SC.energy(pair);
            System.out.println(String.format("%+5.4f  %+18.10e  %+18.10e  %+18.10e   ", i*dx, u, uti, usc));
//            System.out.print(String.format("%+5.4f  %+18.10e  %+10.4e   ", i*dx, p2.energy(pair), (u-lu)/dx));
//            Vector[][] gradientAndTorque = p2.gradientAndTorque(pair);
//            double g = gradientAndTorque[0][1].getX(0);
//            System.out.print(String.format("%+10.4e  %+10.4e  %+10.4e\n", g, (g-lg)/dx, d2[1].component(0,0)));
//            lg = g;
//            lu = u;
            ((Orientation3D)atom1.getOrientation()).rotateBy(dx, y);
            ((Orientation3D)atom0.getOrientation()).rotateBy(dx, z);
        }

//        atom1.getPosition().setX(0, 5);
//        
//        for (int i=0; i<=4; i++) {
////            atom1.getPosition().setX(1, i*0.001);
//            double u = p2.energy(pair);
//            System.out.print(String.format("%+5.4f  %+18.10e  %+10.4e   ", i*dx, p2.energy(pair), (u-lu)/dx));
//            Vector[][] gradientAndTorque = p2.gradientAndTorque(pair);
//            Tensor[] d2 = p2.secondDerivative(pair);
//            double g = -gradientAndTorque[1][1].getX(1);
//            System.out.print(String.format("%+10.4e  %+10.4e  %+10.4e\n", g, (g-lg)/dx, d2[1].component(4,4)));
//            lg = g;
//            lu = u;
//            ((Orientation3D)atom1.getOrientation()).rotateBy(dx, y);
//        }
    }
    
    /*
     * Randomly moves molcules all over, computing and checking 2nd derivative
     */
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(200);
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating speciesCO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(speciesCO2);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 2);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
//        ((IAtomOriented)pair.getAtom(0)).getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
        Vector o1 = space.makeVector(new double[]{-1,0,0});
        atom1.getOrientation().setDirection(o1);
        P2CO2Hellmann p2 = new P2CO2Hellmann(space, Parameters.B);
        P2CO2SC p2SC = p2.makeSemiclassical(temperature);
        System.out.println("or: "+atom0.getOrientation().getDirection()+" "+atom1.getOrientation().getDirection());
        double lu = 0, lg = 0;
        double dx = 0.0001;
        Vector x = null;
        Vector y = space.makeVector(new double[]{0,1,0});
        Vector z = space.makeVector(new double[]{0,0,1});
        atom1.getPosition().setX(0, 5);
//        ((OrientationFull3D)atom1.getOrientation()).rotateBy(-Math.atan2(p2.sitesHH*0.5,p2.sitesOH-cmx), z);
//        ((OrientationFull3D)atom1.getOrientation()).rotateBy(Math.PI/2, y);
        Vector vdx = space.makeVector();
//        ((OrientationFull3D)atom1.getOrientation()).rotateBy(0.28, z);

        RandomMersenneTwister random = new RandomMersenneTwister(5);
        Vector[] xyzAxes = new Vector[]{x,y,z};
        double[] u = new double[3];
        for (int j=0; j<100; j++) {
            int imol = random.nextInt(2);
//            imol = 1;
            IAtomOriented iAtom = imol == 0 ? atom0 : atom1;
            x = iAtom.getOrientation().getDirection();
            xyzAxes[0] = x;
            p2SC.getPerp(x, y, z);
//            System.out.println(x+" "+y+" "+z);
            boolean rot = random.nextInt(2) == 0;
//            rot = true;
            double d2 = 0;
            int xyz = random.nextInt(3);
//            xyz = 1;
            vdx.E(xyzAxes[xyz]);
            if (random.nextInt(2)==0) {
                if (rot) {
                    ((IOrientation3D)iAtom.getOrientation()).rotateBy(0.5, vdx);
                }
                else {
                    iAtom.getPosition().PEa1Tv1(0.5, vdx);
                }
                continue;
            }
            if (rot) {
//                System.out.println("vdx "+vdx);
//                if (imol==1 && xyz==0) dx *= 0.1;
                for (int i=0; i<3; i++) {
                    ((IOrientation3D)iAtom.getOrientation()).rotateBy(dx, vdx);
                    u[i] = p2.energy(pair);
                    
//                    if (imol==1 && xyz==0) System.out.println(u[i]);
                    if (i==1) {
                        p2SC.energy(pair);
                        d2 = p2SC.d2tot[imol][3+xyz];
                    }
                }
            }
            else {
                for (int i=0; i<3; i++) {
                    iAtom.getPosition().PEa1Tv1(dx, vdx);
                    u[i] = p2.energy(pair);
                    if (i==1) {
                        p2SC.energy(pair);
                        d2 = p2SC.d2tot[imol][xyz];
                    }
                }
            }
            double d2fd = (u[0] - 2*u[1] + u[2])/(dx*dx);
            if (d2==0 && d2fd==0) continue;
            double check = (d2fd-d2)/(0.5*(Math.abs(d2)+Math.abs(d2fd)));
            System.out.print(String.format("%d %d %d %d %+10.4e %+10.4e %10.4e", j, imol, rot?1:0, xyz, d2, d2fd, check));
            if (Math.abs(check) > 0.01) {
                System.out.println("  oops");
                if (!rot) throw new RuntimeException("oops");
            }
            else {
                System.out.println();
            }
        }
    }
}
