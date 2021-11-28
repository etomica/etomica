/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.models.water;

import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.potential.IPotential2;
import etomica.potential.Potential3Soft;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.units.Hartree;
import etomica.util.Constants;

/**
 * Semiclassical potential class for Szalewicz water model.  When invoked as a 3-body
 * potential, only the pair contributions are handled semiclassically.
 */
public class P2H2OSzalewiczSC extends P2WaterSzalewicz {

    protected final double[] drvDamp;
    protected final Vector[][] gi;
    protected final Tensor tt0Tensor, tt1Tensor, rr0Tensor, rr1Tensor;
    protected final Tensor ijTensor, rTensor0, rTensor1, identity;
    protected final Tensor ijRTensor;
    protected final Tensor[] rot;
    protected final Vector drijRot;
    public double[][] d2tot;
    protected final double temperature, fac;
    protected final Vector drij;
    protected final Vector moment;
    protected final Vector ri, rj;

    public P2H2OSzalewiczSC(double temperature, int nBody) {
        super(nBody);
        drvDamp = new double[6];
        drij = space.makeVector();
        ijTensor = space.makeTensor();
        identity = space.makeTensor();
        tt0Tensor = space.makeTensor();
        tt1Tensor = space.makeTensor();
        rr0Tensor = space.makeTensor();
        rr1Tensor = space.makeTensor();
        rTensor0 = space.makeTensor();
        rTensor1 = space.makeTensor();
        ijRTensor = space.makeTensor();
        identity.E(new double[][]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
        gi = new Vector[2][sites.length];
        for (int i = 0; i < sites.length; i++) {
            gi[0][i] = space.makeVector();
            gi[1][i] = space.makeVector();
        }
        drijRot = space.makeVector();
        rot = new Tensor[sitePos.length];
        for (int i = 0; i < rot.length; i++) {
            rot[i] = space.makeTensor();
        }

        this.temperature = temperature;
        double hbar = Constants.PLANCK_H / (2 * Math.PI);
        fac = hbar * hbar / (24 / 2) / temperature;

        moment = space.makeVector();
        moment.setX(0, 2 * Hydrogen.INSTANCE.getMass() * siteDoubles[1][2] * siteDoubles[1][2] + Oxygen.INSTANCE.getMass() * siteDoubles[0][2] * siteDoubles[0][2]);
        moment.setX(1, 2 * Hydrogen.INSTANCE.getMass() * sites[1].squared() + Oxygen.INSTANCE.getMass() * sites[0].squared());
        moment.setX(2, 2 * Hydrogen.INSTANCE.getMass() * siteDoubles[1][0] * siteDoubles[1][0]);

        ri = space.makeVector();
        rj = space.makeVector();
        d2tot = new double[sitePos.length][6];
    }

    public static P2H2OSzalewiczSC.P2H2OSzalewiczSC2 make2BodySC(double temperature) {
        return make2BodySC(temperature, Component.ALL);
    }

    public static P2H2OSzalewiczSC.P2H2OSzalewiczSC2 make2BodySC(double temperature, Component comp) {
        return new P2H2OSzalewiczSC.P2H2OSzalewiczSC2(temperature, comp);
    }

    public static class P2H2OSzalewiczSC2 extends P2H2OSzalewiczSC implements IPotential2 {
        public P2H2OSzalewiczSC2(double temperature, Component comp) {
            super(temperature, 2);
            setComponent(comp);
        }

        @Override
        public double u(double r2) {
            throw new RuntimeException("Not a spherical potential");
        }

        @Override
        public double u(Vector dr12, IAtom atom1, IAtom atom2) {
            Vector zero = space.makeVector();
            return energy(new IAtom[]{atom1, atom2}, new Vector[]{zero, dr12});
        }
    }

    public static P2H2OSzalewiczSC.P2H2OSzalewiczSC3 make3BodySC(double temperature) {
        return make3BodySC(temperature, Component.ALL);
    }

    public static P2H2OSzalewiczSC.P2H2OSzalewiczSC3 make3BodySC(double temperature, Component comp) {
        return new P2H2OSzalewiczSC.P2H2OSzalewiczSC3(temperature, comp);
    }

    public static class P2H2OSzalewiczSC3 extends P2H2OSzalewiczSC implements Potential3Soft {
        public P2H2OSzalewiczSC3(double temperature, Component comp) {
            super(temperature, 3);
            setComponent(comp);
        }

        @Override
        public double u(double r212, double r213, double r223) {
            throw new RuntimeException("Not a spherical potential");
        }

        @Override
        public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3) {
            Vector zero = space.makeVector();
            return energy(new IAtom[]{atom1, atom2, atom3}, new Vector[]{zero, dr12, dr13});
        }
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IAtom[] atoms, Vector[] pos) {
        for (int i = 0; i < atoms.length; i++) {
            rTmp.Ea1Tv1(bohrConv, atoms[i].getPosition());
            // everything after this is in Bohr
            OrientationFull3D ori = (OrientationFull3D) ((IAtomOriented) atoms[i]).getOrientation();
            allOr[0] = ori.getDirection();
            allOr[1] = ori.getSecondaryDirection();
            or2.E(allOr[0]);
            or2.XE(allOr[1]);
            for (int j = 0; j < sites.length; j++) {
                sitePos[i][j].E(rTmp);
                for (int k = 0; k < 3; k++) {
                    sitePos[i][j].PEa1Tv1(sites[j].getX(k), allOr[k]);
                }
            }

            rot[i].E(allOr);
            rot[i].invert();

            com0[i].E(rTmp);
            // according to our atomic masses, the COM is not at 0,0,0
            // we shifted the positions from siteDoubles to sitePos so that the COM is at 0,0,0
            // however, the original COM is used for 3-body interactions.  compute that now
            com0[i].PEa1Tv1(cmzFix, allOr[2]);
        }

        double u0 = 0;
        double d2tsum = 0, d2rsum = 0;
        if (component != Component.INDUCTION && component != Component.NON_PAIR && component != Component.THREE_BODY) {
            for (int iA = 0; iA < atoms.length - 1; iA++) {
                for (int iB = iA + 1; iB < atoms.length; iB++) {
                    tt0Tensor.E(0);
                    tt1Tensor.E(0);
                    rr0Tensor.E(0);
                    rr1Tensor.E(0);
                    for (int i = 0; i < sites.length; i++) {
                        gi[0][i].E(0);
                        gi[1][i].E(0);
                    }
                    for (int nsA = 0; nsA < sites.length; nsA++) {
                        rTensor0.setComponent(0, 1, siteDoubles[nsA][2]);
                        rTensor0.setComponent(1, 0, -siteDoubles[nsA][2]);
                        rTensor0.setComponent(0, 2, -siteDoubles[nsA][1]);
                        rTensor0.setComponent(2, 0, siteDoubles[nsA][1]);
                        rTensor0.setComponent(1, 2, siteDoubles[nsA][0]);
                        rTensor0.setComponent(2, 1, -siteDoubles[nsA][0]);
                        for (int nsB = 0; nsB < sites.length; nsB++) {
                            rTensor1.setComponent(0, 1, siteDoubles[nsB][2]);
                            rTensor1.setComponent(1, 0, -siteDoubles[nsB][2]);
                            rTensor1.setComponent(0, 2, -siteDoubles[nsB][1]);
                            rTensor1.setComponent(2, 0, siteDoubles[nsB][1]);
                            rTensor1.setComponent(1, 2, siteDoubles[nsB][0]);
                            rTensor1.setComponent(2, 1, -siteDoubles[nsB][0]);

                            drij.Ev1Mv2(sitePos[iB][nsB], sitePos[iA][nsA]);
//                                if (nsA==1&&nsB==2) System.out.println(nsA+" "+nsB+" "+drij);
                            double R2 = drij.squared();
                            double R = Math.sqrt(R2);
                            if (R < core) return Double.POSITIVE_INFINITY;
                            double u = 0;
                            double rdudr = 0;
                            double r2d2udr2 = 0;

                            int ib = ind_beta[nsA][nsB] - 1;
                            if (ib != 0) {  // exponential
                                double beta = params[ib];
                                double eks = Math.exp(-beta * R);
                                int indlin = ib - 98;
                                if (indlin < 0) indlin = indlin + 65;
                                int ind0 = indlin;
                                int ind1 = ind0 + 36;
                                int ind2 = ind1 + 36;
                                int ind3 = ind2 + 36;
//                                   aj[ind0] = aj[ind0] +eks;
//                                   aj[ind1] = aj[ind1] +eks*R;
//                                   aj[ind2] = aj[ind2] +eks*R*R;
//                                   aj[ind3] = aj[ind3] +eks*R*R*R;
                                // damp the site-site u_exp function if u_exp is negative
                                double u_exp = eks * (c[ind0] + c[ind1] * R + c[ind2] * R2 + c[ind3] * R2 * R);
                                if (u_exp < 0) {
                                    double fdamp = 1 / (1 + Math.exp(-gamma * (R - R0)));
//                                       aj[ind0] = aj[ind0] *fdamp;
//                                       aj[ind1] = aj[ind1] *fdamp;
//                                       aj[ind2] = aj[ind2] *fdamp;
//                                       aj[ind3] = aj[ind3] *fdamp;
//                                       System.out.println("usc "+nsA+" "+nsB+" "+u_exp+" "+fdamp+" "+u_exp*fdamp);
                                    u += u_exp * fdamp;
                                    double blah = -gamma * R * Math.exp(-gamma * (R - R0));
                                    rdudr += (eks * (c[ind1] * R + 2 * c[ind2] * R2 + 3 * c[ind3] * R2 * R) - beta * R * u_exp) * (fdamp) +
                                            -u_exp * fdamp * fdamp * blah;
                                    r2d2udr2 += (eks * (2 * c[ind2] * R2 + 6 * c[ind3] * R2 * R) - 2 * beta * R * eks * (c[ind1] * R + 2 * c[ind2] * R2 + 3 * c[ind3] * R2 * R) + beta * beta * R2 * u_exp) * fdamp
                                            - 2 * (eks * (c[ind1] * R + 2 * c[ind2] * R2 + 3 * c[ind3] * R2 * R) - beta * R * u_exp) * fdamp * fdamp * blah
                                            + 2 * u_exp * fdamp * fdamp * fdamp * blah * blah + u_exp * fdamp * fdamp * blah * gamma * R;
//                                       System.out.println(R/bohrConv+" "+u+" "+rdudr/R*bohrConv+" "+r2d2udr2/R2*bohrConv*bohrConv);
                                } else {
                                    u += u_exp;
                                    rdudr += eks * (c[ind1] * R + 2 * c[ind2] * R2 + 3 * c[ind3] * R2 * R) - beta * R * u_exp;
                                    r2d2udr2 += eks * (2 * c[ind2] * R2 + 6 * c[ind3] * R2 * R) - 2 * beta * R * eks * (c[ind1] * R + 2 * c[ind2] * R2 + 3 * c[ind3] * R2 * R) + beta * beta * R2 * u_exp;
//                                       System.out.println(R/bohrConv+" "+u+" "+rdudr/R*bohrConv+" "+r2d2udr2/R2*bohrConv*bohrConv);
                                }
                                if (Double.isNaN(u)) {
                                    throw new RuntimeException("oops");
                                }
                            }

                            if (nsA < 5 && nsB < 5) { // elst
                                double qA = params[ind_charge[nsA]];
                                double qB = params[ind_charge[nsB]];
                                double d1 = params[ind_d1[nsA][nsB]];
                                double[] f = ddamp(1, d1, R);
                                double br = d1 * R;
                                double expbr = Math.exp(-br);
                                u += f[3] * qA * qB / R;
                                rdudr += -(-expbr * f[4] + 1 * f[3]) * qA * qB / R;
                                r2d2udr2 += -(expbr * br * f[4] - expbr * f[5] + 1 * (-f[3] + expbr * f[4]) + -1 * (-expbr * f[4] + 1 * f[3])) * qA * qB / R;
//                                   if (nsA+nsB==0) System.out.println(R/bohrConv+" "+Hartree.UNIT.toSim(u)+" "+Hartree.UNIT.toSim(rdudr)/R*bohrConv+" "+Hartree.UNIT.toSim(r2d2udr2)/R2*bohrConv*bohrConv);
//                                   if (nsA+nsB==0) System.out.println(R+" "+u+" "+rdudr+" "+r2d2udr2);

//                                   double rdu6dr = (-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum))*C6[ii][jj]/rij6;
//                                   double r2du26dr2 = (expbr*br*br*sum-expbr*br*dsum -br*expbr*dsum+expbr*d2sum + 6*(-1 + br*expbr*sum-expbr*dsum+expbr*sum) + -6*(-expbr*br*sum + expbr*dsum + 6*(1-expbr*sum)))*C6[ii][jj]/rij6;

                                if (Double.isNaN(u)) {
                                    throw new RuntimeException("oops");
                                }
                            }


                            if (nsA < 3 && nsB < 3) {    // ind-disp
                                double d6 = params[ind_d6[nsA][nsB]];
                                double d8 = params[ind_d8[nsA][nsB]];
                                double d10 = params[ind_d10[nsA][nsB]];
                                double C6 = params[ind_C6[nsA][nsB]];
                                double C8 = params[ind_C8[nsA][nsB]];
                                double C10 = params[ind_C10[nsA][nsB]];
                                double R6 = R2 * R2 * R2;
                                double R8 = R6 * R2;
                                double R10 = R8 * R2;

                                double[] f = ddamp(6, d6, R);
                                // f3 = 1-expbr*sum
                                // f4 = f0*br - f1
                                // f5 = f1*br - f2
                                u += -f[3] * C6 / R6;
                                double br = d6 * R;
                                double expbr = Math.exp(-br);
                                rdudr += (-expbr * f[4] + 6 * f[3]) * C6 / R6;
                                r2d2udr2 += (expbr * br * f[4] - expbr * f[5] + 6 * (-f[3] + expbr * f[4]) + -6 * (-expbr * f[4] + 6 * f[3])) * C6 / R6;

                                f = ddamp(8, d8, R);
                                br = d8 * R;
                                expbr = Math.exp(-br);
                                u += -f[3] * C8 / R8;
                                rdudr += (-expbr * f[4] + 8 * f[3]) * C8 / R8;
                                r2d2udr2 += (expbr * br * f[4] - expbr * f[5] + 8 * (-f[3] + expbr * f[4]) + -8 * (-expbr * f[4] + 8 * f[3])) * C8 / R8;

                                f = ddamp(10, d10, R);
                                u += -f[3] * C10 / R10;
                                br = d10 * R;
                                expbr = Math.exp(-br);
//                                    System.out.println(d10+" "+R+" "+expbr*br*f[0]+" "+expbr*f[1]+" "+10*f[3]+" "+(-expbr*br*f[0] + expbr*f[1] + 10*f[3]));
                                rdudr += (-expbr * f[4] + 10 * f[3]) * C10 / R10;
                                r2d2udr2 += (expbr * br * f[4] - expbr * f[5] + 10 * (-f[3] + expbr * f[4]) + -10 * (-expbr * f[4] + 10 * f[3])) * C10 / R10;
//                                    System.out.println(R/bohrConv+" "+Hartree.UNIT.toSim(u)+" "+Hartree.UNIT.toSim(rdudr)/R*bohrConv+" "+Hartree.UNIT.toSim(r2d2udr2)/R2*bohrConv*bohrConv);
                                if (Double.isNaN(u)) {
                                    throw new RuntimeException("oops");
                                }
                            }

                            u0 += u;
                            if (Double.isNaN(u0)) {
                                throw new RuntimeException("oops");
                            }
//                                System.out.println("u0sc "+u0);

                            // molecule 0
                            drijRot.E(drij);
                            rot[iA].transform(drijRot);
                            ijTensor.Ev1v2(drijRot, drijRot);
                            ijTensor.TE((rdudr - r2d2udr2) / (R2 * R2));
                            ijTensor.PEa1Tt1(-rdudr / R2, identity);

                            tt0Tensor.ME(ijTensor);

                            rTensor0.transpose();
                            ijRTensor.E(rTensor0);
                            ijRTensor.TE(ijTensor);
                            rTensor0.transpose();
                            ijRTensor.TE(rTensor0);
                            rr0Tensor.ME(ijRTensor);

                            drijRot.TE(rdudr / R2);
                            gi[0][nsA].ME(drijRot);
                            if (gi[0][nsB].isNaN()) {
                                throw new RuntimeException("oops");
                            }


                            // molecule 1
                            drijRot.E(drij);
                            rot[iB].transform(drijRot);
                            ijTensor.Ev1v2(drijRot, drijRot);
//                                System.out.println("r2 "+rij2+" "+ijTensor.component(0,0));
                            ijTensor.TE((rdudr - r2d2udr2) / (R2 * R2));
                            ijTensor.PEa1Tt1(-rdudr / R2, identity);

                            // we only need tt for 0, both have the same contribution
                            tt1Tensor.ME(ijTensor);
//                                System.out.println("d2 "+r2d2udr2/rij2+" "+ijTensor.component(0,0)+" "+tt1Tensor.component(0,0));

                            rTensor1.transpose();
                            ijRTensor.E(rTensor1);
                            ijRTensor.TE(ijTensor);
                            rTensor1.transpose();
                            ijRTensor.TE(rTensor1);
                            rr1Tensor.ME(ijRTensor);

                            drijRot.TE(rdudr / R2);
                            gi[1][nsB].PE(drijRot);
                            if (gi[1][nsB].isNaN()) {
                                throw new RuntimeException("oops");
                            }

                        }
                    }
                    for (int i = 0; i < sites.length; i++) {
                        rr0Tensor.PEv1v2(sites[i], gi[0][i]);
                        // we really just want to modify the diagonal, but we don't care
                        // abuot the off-diagonal
                        rr0Tensor.PE(-sites[i].dot(gi[0][i]));

                        rr1Tensor.PEv1v2(sites[i], gi[1][i]);
                        rr1Tensor.PE(-sites[i].dot(gi[1][i]));
                    }
                    for (int i = 0; i < 3; i++) {
                        d2tot[iA][i] = tt0Tensor.component(i, i);
                        d2tot[iB][i] = tt1Tensor.component(i, i);
                        d2tsum += tt0Tensor.component(i, i) / mass;
                    }
                    for (int i = 0; i < 3; i++) {
                        d2tot[iA][3 + i] = rr0Tensor.component(i, i);
                        d2tot[iB][3 + i] = rr1Tensor.component(i, i);
                        d2rsum += (rr0Tensor.component(i, i) + rr1Tensor.component(i, i)) / (2 * moment.getX(i));
                    }
                    for (int j = 0; j < 3; j++) {
                        d2tot[iA][j] = Hartree.UNIT.toSim(d2tot[iA][j] * bohrConv * bohrConv);
                        d2tot[iA][3 + j] = Hartree.UNIT.toSim(d2tot[iA][3 + j]);
                        d2tot[iB][j] = Hartree.UNIT.toSim(d2tot[iB][j] * bohrConv * bohrConv);
                        d2tot[iB][3 + j] = Hartree.UNIT.toSim(d2tot[iB][3 + j]);
                    }
                }
            }
        }
        if (Hartree.UNIT.toSim(u0) > 10000) {
            return Double.POSITIVE_INFINITY;
        }
        double uInd = 0;
        if (component != Component.TWO_BODY && component != Component.THREE_BODY) {
            uInd = indN_iter();
        }
        double u3 = 0;
        if (component != Component.TWO_BODY && component != Component.INDUCTION) {
            for (int i = 0; i < atoms.length - 2; i++) {
                com3[0] = com0[i];
                sitePos3[0] = sitePos[i];
                for (int j = i + 1; j < atoms.length - 1; j++) {
                    com3[1] = com0[j];
                    sitePos3[1] = sitePos[j];
                    for (int k = j + 1; k < atoms.length; k++) {
                        com3[2] = com0[k];
                        sitePos3[2] = sitePos[k];
                        u3 += fit3b();
                    }
                }
            }
        }
        // semiclassical correction only accounts for u0
        if (Math.random() < 0.001 && false) {
            System.out.println(String.format("% 10.4e  % 10.4e  % 10.4e  % 10.4e", Hartree.UNIT.toSim(u0), Hartree.UNIT.toSim(uInd), Hartree.UNIT.toSim(u3), Hartree.UNIT.toSim(fac * (d2tsum * bohrConv * bohrConv + d2rsum))));
        }
        return Hartree.UNIT.toSim(u0 + uInd + u3 + fac * (d2tsum + d2rsum) * bohrConv * bohrConv);
    }

    protected double[] ddamp(int n, double beta, double r) {
        double br = beta * r;
        drvDamp[0] = 1;
        drvDamp[1] = drvDamp[2] = 0;
        drvDamp[3] = 0;

        if (br == 0) {
            return drvDamp;
        }
        double sum = 1;
        double dsum = 0;
        double d2sum = 0;
        double term = 1;
        double dterm = 1;
        double d2term = 0;
        for (int i = 1; i <= n; i++) {
            term *= br / i;
            if (i == 1) {
                dterm = br;
            } else {
                dterm *= br / (i - 1);
            }
            if (i == 2) {
                d2term = br * br;
            } else if (i > 2) {
                d2term *= br / (i - 2);
            }
            sum += term;
            dsum += dterm;
            d2sum += d2term;
        }

        drvDamp[0] = sum;
        drvDamp[1] = dsum;
        drvDamp[2] = d2sum;
        drvDamp[3] = 1.0 - Math.exp(-br) * sum;
        drvDamp[4] = term * br;  // f0*br - f1
        drvDamp[5] = dterm * br;  // f1*br - f2
        // in case of d --> 0 use
        // d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
        if (Math.abs(drvDamp[3]) < 1.0e-8) {
            drvDamp[3] = 0.0;
            for (int i = n + 1; i <= 1000; i++) {
                term = term * br / i;
                drvDamp[3] += term;
                if (term / drvDamp[3] < 1.0e-8) break;
            }
            if (term / drvDamp[3] > 1e-8) throw new RuntimeException("No convergence in damp" + " " + r);
            drvDamp[3] *= Math.exp(-br);
        }
        //     write(6,'(i4,2f10.5,e20.10)') n,beta,r,d
        return drvDamp;
    }

}
