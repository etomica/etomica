/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.units.*;
import etomica.util.collections.IntArrayList;

public class PotentialMoleculePairAmoebaMPole implements IPotentialMolecular {
    protected final Space space;
    protected final IntArrayList[][] bonding;
    protected final Axis[][] axis;
    protected final Multipole[][] multipoles;

    public PotentialMoleculePairAmoebaMPole(Space space, IntArrayList[][] bonding, Axis[][] axis, Multipole[][] multipoles) {
        this.space = space;
        this.bonding = bonding;
        this.axis = axis;
        this.multipoles = multipoles;
    }

    public static Tensor getRotmat(IAtom a, Axis myAxisAtoms) {
        IAtomList atoms = a.getParentGroup().getChildList();
//        System.out.println(a.getLeafIndex()+" "+myAxisAtoms[0]+" "+myAxisAtoms[1]);
        Vector3D xAxis = new Vector3D();
        Vector3D yAxis = new Vector3D();
        Vector3D zAxis = new Vector3D();
        switch (myAxisAtoms.scheme) {
            case ZONLY:
                IAtom zAtom = atoms.get(myAxisAtoms.atoms[0]);
                zAxis.Ev1Mv2(zAtom.getPosition(), a.getPosition());
                zAxis.normalize();
                xAxis.E(1,0,0);
                double dot = zAxis.getX(0);
                if (Math.abs(dot) > 0.707) {
                    xAxis.E(0,1,0);
                    dot = zAxis.getX(1);
                }
                xAxis.PEa1Tv1(-dot, zAxis);
                xAxis.normalize();
                break;
            case ZTHENX:
                zAtom = atoms.get(myAxisAtoms.atoms[0]);
                zAxis.Ev1Mv2(zAtom.getPosition(), a.getPosition());
                zAxis.normalize();
                IAtom xAtom = atoms.get(myAxisAtoms.atoms[1]);
                xAxis.Ev1Mv2(xAtom.getPosition(), a.getPosition());
                xAxis.PEa1Tv1(-xAxis.dot(zAxis), zAxis);
                xAxis.normalize();
                break;
            case ZBISECT:
                // ? define z axis as normal; define x axis with bisector of 2nd and 3rd atoms (perp to z)
                zAtom = atoms.get(myAxisAtoms.atoms[0]);
                zAxis.Ev1Mv2(zAtom.getPosition(), a.getPosition());
                zAxis.normalize();
                xAtom = atoms.get(myAxisAtoms.atoms[1]);
                xAxis.Ev1Mv2(xAtom.getPosition(), a.getPosition());
                xAxis.normalize();
                IAtom yAtom = atoms.get(myAxisAtoms.atoms[2]);
                yAxis.Ev1Mv2(yAtom.getPosition(), a.getPosition());
                yAxis.normalize();
                xAxis.PE(yAxis);
                xAxis.normalize();
                dot = zAxis.dot(xAxis);
                if (false) {
                    // ???? fortran checks 180 - 180/pi acos(dot) < 0
                    // trying to check for z and x axes pointing in exactly opposite directions????
                }
                xAxis.PEa1Tv1(-dot, zAxis);
                xAxis.normalize();
                break;
            case THREEFOLD:
                // ? z points in x+y+z direction, x points toward original z, perp to x+y+z
                zAxis.Ev1Mv2(atoms.get(myAxisAtoms.atoms[0]).getPosition(), a.getPosition());
                zAxis.normalize();
                Vector3D zz = new Vector3D(zAxis);
                xAxis.Ev1Mv2(atoms.get(myAxisAtoms.atoms[1]).getPosition(), a.getPosition());
                xAxis.normalize();
                zAxis.PE(xAxis);
                xAxis.Ev1Mv2(atoms.get(myAxisAtoms.atoms[2]).getPosition(), a.getPosition());
                xAxis.normalize();
                zAxis.PE(xAxis);
                zAxis.normalize();
                dot = zz.dot(zAxis);
                xAxis.E(zz);
                xAxis.PEa1Tv1(-dot, zAxis);
                xAxis.normalize();
                break;
            default:
                throw new RuntimeException("Unknown axis scheme");
        }
        yAxis.E(xAxis);
        yAxis.XE(zAxis);
//        System.out.println("z "+zAtom.getPosition()+" "+zAxis);
        Tensor t = new Tensor3D();
        t.E(new Vector[]{xAxis, yAxis, zAxis});
        t.setComponent(0, 2, zAxis.getX(0));
        t.setComponent(1, 2, zAxis.getX(1));
        t.setComponent(2, 2, zAxis.getX(2));

        t.setComponent(0, 0, xAxis.getX(0));
        t.setComponent(1, 0, xAxis.getX(1));
        t.setComponent(2, 0, xAxis.getX(2));

        t.setComponent(0, 1, yAxis.getX(0));
        t.setComponent(1, 1, yAxis.getX(1));
        t.setComponent(2, 1, yAxis.getX(2));

//        if (a.getLeafIndex()==7) System.out.println("rotmat"+a.getLeafIndex()+" "+myAxisAtoms.scheme+"\n"+t);

        return t;
    }

    protected Tensor transformQ(Tensor Q, Tensor rot) {
        Tensor A = new Tensor3D();
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                double sum = 0;
                for (int k=0; k<3; k++) {
                    for (int m=0; m<3; m++) {
                        sum += rot.component(i,k)*rot.component(j,m)*Q.component(m,k);
                    }
                }
                A.setComponent(i,j,sum);
            }
        }
        return A;
    }

    protected double doubledot(Tensor A, Tensor B) {
        Tensor x = new Tensor3D();
        x.E(A);
        x.TE(B);
        return x.trace();
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {
        IAtomList atoms1 = molecule1.getChildList();
        IAtomList atoms2 = molecule2.getChildList();
        double[] charges2 = new double[atoms2.size()];
        Vector[] dipoles2 = new Vector[atoms2.size()];
        double[] muMag2 = new double[atoms1.size()];
        Tensor[] quadrupoles2 = new Tensor[atoms2.size()];
        double u = 0;

        for (IAtom a2 : atoms2) {
            Axis myAxis2 = axis[a2.getParentGroup().getType().getIndex()][a2.getIndex()];
            Tensor rot = getRotmat(a2, myAxis2);
//            System.out.println("rot "+a2.getLeafIndex()+"\n"+rot);
            int i2 = a2.getIndex();
            Multipole myMultipole = multipoles[molecule2.getType().getIndex()][i2];
            charges2[i2] = myMultipole.q;
            muMag2[i2] = Math.sqrt(myMultipole.mu.squared());
            dipoles2[i2] = new Vector3D();
//            System.out.println("mu "+a2.getLeafIndex()+" "+Electron.UNIT.fromSim(myMultipole.mu.getX(2)));
            dipoles2[i2].E(myMultipole.mu);
            rot.transform(dipoles2[i2]);
//            System.out.println("mu "+a2.getLeafIndex()+" "+Electron.UNIT.fromSim(dipoles2[i2].getX(2)));
//            System.out.println("Q "+a2.getLeafIndex()+"\n"+myMultipole.Q);
            quadrupoles2[i2] = transformQ(myMultipole.Q, rot);
//            System.out.println("Qxx "+a2.getLeafIndex()+" "+Electron.UNIT.fromSim(quadrupoles2[i2].component(0,0)));

        }
        for (IAtom a1 : atoms1) {
            int i1 = a1.getIndex();
            Axis myAxis1 = axis[a1.getParentGroup().getType().getIndex()][i1];
            Tensor rot = getRotmat(a1, myAxis1);
            Multipole myMultipole = multipoles[molecule2.getType().getIndex()][i1];
            double q1 = myMultipole.q;
            double muMag1 = Math.sqrt(myMultipole.mu.squared());
            Vector mu1 = new Vector3D();
            mu1.E(myMultipole.mu);
            rot.transform(mu1);
            Tensor Q1 = transformQ(myMultipole.Q, rot);

//            System.out.println("Q1 "+ a1.getLeafIndex());
//            Tensor3D QQ1 = new Tensor3D();
//            QQ1.E(myMultipole.Q);
//            QQ1.TE(Electron.UNIT.fromSim(1));
//            System.out.println(QQ1);
//            System.out.println("rot");
//            System.out.println(rot);
//            System.out.println("Q1 => "+ a1.getLeafIndex());
//            QQ1.E(Q1);
//            QQ1.TE(Electron.UNIT.fromSim(1));
//            System.out.println(QQ1);

            for (IAtom a2 : atoms2) {
                double q2 = charges2[a2.getIndex()];
                Vector mu2 = dipoles2[a2.getIndex()];
                Tensor Q2 = quadrupoles2[a2.getIndex()];
                Vector dr = space.makeVector();
                dr.Ev1Mv2(a2.getPosition(), a1.getPosition());
                double r2 = dr.squared();
                double r = Math.sqrt(r2);

                double uqq = q1*q2/r;

                double mmcos_D1_D2 = mu1.dot(mu2);
                double cos_D1_D2 = mmcos_D1_D2 / muMag1 / muMag2[a2.getIndex()];
                //cos(dipole 1 and r12)
                double c1 = mu1.dot(dr)/r;  // dir/r
                //cos(r12 and dipole 2)
                double c2 = dr.dot(mu2)/r;  // dkr/r
                double uqmu = -q1*c2/r2;
                double umuq = q2*c1/r2;
                double umumu = (mmcos_D1_D2 - 3.0 * c1* c2) / (r*r2);
                Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);

                Vector Q1dr = new Vector3D();
                Q1dr.E(dr);
                Q1.transform(Q1dr);
                double Q1drdr = Q1dr.dot(dr); // qir
                double Q1drmu2 = Q1dr.dot(mu2); // dkqi

                Vector Q2dr = new Vector3D();
                Q2dr.E(dr);
                Q2.transform(Q2dr);
                double Q2drdr = Q2dr.dot(dr); // qkr
                double Q2drmu1 = Q2dr.dot(mu1); // diqk

                double Q1drQ2dr = Q1dr.dot(Q2dr); // qik

                double uqQ = 3*q1*Q2drdr/(r*r2*r2);
                double uQq = 3*Q1drdr*q2/(r*r2*r2);

                double umuQ = -2*3*Q2drmu1/(r*r2*r2) + 3*5*c1*Q2drdr/(r2*r2*r2);
                double uQmu = 2*3*Q1drmu2/(r*r2*r2) - 3*5*c2*Q1drdr/(r2*r2*r2);

//                System.out.println("Q1 "+ a1.getLeafIndex());
//                Tensor3D QQ1 = new Tensor3D();
//                QQ1.E(Q1);
//                QQ1.TE(Electron.UNIT.fromSim(1));
//                System.out.println(QQ1);
//                Tensor3D QQ2 = new Tensor3D();
//                QQ2.E(Q2);
//                QQ2.TE(Electron.UNIT.fromSim(1));
//                System.out.println("Q2 "+ a2.getLeafIndex());
//                System.out.println(QQ2);

                double uQQ = 2*3*doubledot(Q1, Q2)/(r*r2*r2)
                           - 4*3*5*Q1drQ2dr/(r*r2*r2*r2) + 3*5*7*Q1drdr*Q2drdr/(r*r2*r2*r2*r2);
//                if (i1==0 && a2.getLeafIndex() == 7) {
//                    System.out.println("QQ " + kcalpmole.fromSim(2 * 3 * doubledot(Q1, Q2) / (r * r2 * r2)) + " " +
//                            kcalpmole.fromSim(-4 * 3 * 5 * Q1drQ2dr / (r * r2 * r2 * r2)) +
//                            kcalpmole.fromSim(3 * 5 * 7 * Q1drdr * Q2drdr / (r * r2 * r2 * r2 * r2)));
//                }

                u += uqq + uqmu + umuq + umumu + uqQ + uQq + umuQ + uQmu + uQQ;
//                double uu = uqq + uqmu + umuq + umumu + uqQ + uQq + umuQ + uQmu + uQQ;

//                if (i1==0 && a2.getLeafIndex() == 7) {
//                    System.out.printf("mp %d %d  % 10.5e  % 10.5e % 10.5e  % 10.5e  % 10.5e % 10.5e  % 10.5e % 10.5e  % 10.5e  % 10.5e\n", a1.getLeafIndex(), a2.getLeafIndex(),
//                            kcalpmole.fromSim(uqq), kcalpmole.fromSim(umuq), kcalpmole.fromSim(uqmu),
//                            kcalpmole.fromSim(umumu), kcalpmole.fromSim(uqQ), kcalpmole.fromSim(uQq),
//                            kcalpmole.fromSim(umuQ), kcalpmole.fromSim(uQmu), kcalpmole.fromSim(uQQ), kcalpmole.fromSim(uu));
//                }

                if (r < 0.1 && u != 0) return Double.POSITIVE_INFINITY;
//                System.out.printf("mp %d %d  % 10.5e\n", a1.getLeafIndex(), a2.getLeafIndex(), kcalpmole.fromSim(uqq + uqmu + umuq + umumu + uqQ + uQq + umuQ + uQmu + uQQ));
            }
        }
//        Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
//        System.out.println("multipole "+kcalpmole.fromSim(u));
        if (Double.isNaN(u)) throw new RuntimeException();
        return u;
    }

    public static class Multipole {
        public double q;
        public Vector mu;
        public Tensor Q;
    }

    public enum AxisScheme {ZONLY, ZTHENX, BISECTOR, ZBISECT, THREEFOLD};

    public static class Axis {
        public Axis(int[] atoms, AxisScheme scheme) {
            this.atoms = atoms;
            this.scheme = scheme;
        }
        public int[] atoms;
        public AxisScheme scheme;
    }
}
