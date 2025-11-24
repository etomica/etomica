package etomica.potential.mbnrg;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.potential.amoeba.PotentialMoleculePairAmoebaMPole;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.units.Electron;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;

public class PNTTMPolarization implements IPotentialMolecular {

    protected final Space space;
    protected final IntArrayList[][] bonding;
    protected final PotentialMoleculePairAmoebaMPole.Axis[][] axis;
    protected final PotentialMoleculePairAmoebaMPole.Multipole[][] multipoles;
    protected final double[] polarizability;
    protected final double[] damping, pdamp;
    protected double tol;

    public PNTTMPolarization(Space space, IntArrayList[][] bonding, PotentialMoleculePairAmoebaMPole.Axis[][] axis, PotentialMoleculePairAmoebaMPole.Multipole[][] multipoles, double[] polarizability, double[] damping, double tol) {
        this.space = space;
        this.bonding = bonding;
        this.axis = axis;
        this.multipoles = multipoles;
        this.polarizability = polarizability;
        this.damping = damping;
        pdamp = new double[polarizability.length];
        for (int i=0; i<pdamp.length; i++) {
            pdamp[i] = Math.pow(polarizability[i], 1.0/6.0);
        }
        this.tol = tol;
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

    @Override
    public double energy(IMoleculeList molecules) {
        IMolecule[] m = new IMolecule[molecules.size()];
        for (int i=0; i<m.length; i++) m[i] = molecules.get(i);
        return energy(m);
    }

    public double[] dampthole(int i, int k, int order, double r) {
        double[] rv = new double[order+1];
        Arrays.fill(rv, 1);
        double damp = pdamp[i] * pdamp[k];
        if (damp == 0) return rv;
        double pgamma = Math.min(damping[i], damping[k]);
        if (pgamma == 0) pgamma = damping[i] + damping[k];
        if (pgamma == 0) return rv;
        double x = r/damp;
        damp = pgamma * x*x*x;
        if (damp < 50) {
            double expdamp = Math.exp(-damp);
            rv[3] = 1 - expdamp;
            rv[5] = 1 - expdamp*(1+damp);
            if (order >= 7) {
                double damp2 = damp*damp;
                rv[7] = 1 - expdamp*(1 + damp + 0.6*damp2);
                if (order >= 9) {
                    double damp3 = damp*damp2;
                    rv[9] = 1 - expdamp*(1 + damp + 18.0/35.0*damp2 + 9.0/35.0*damp3);
                }
            }
        }
        return rv;
    }

    protected Vector[][] computeField(IMolecule[] molecules, PotentialMoleculePairAmoebaMPole.Multipole[][] mpoles) {
        Vector[][] field = new Vector[mpoles.length][];
        double rmin = Double.POSITIVE_INFINITY;
        for (int i=0; i<molecules.length; i++) {
            IAtomList atoms1 = molecules[i].getChildList();
            field[i] = new Vector[atoms1.size()];

            for (int j=0; j<atoms1.size(); j++) {
                IAtom a1 = atoms1.get(j);
                PotentialMoleculePairAmoebaMPole.Multipole ijmpole = mpoles[i][j];
//                System.out.println(i+" "+j+" "+ijmpole.mu+" "+inducedDipoles[i][j]+" "+mu1);
                field[i][j] = new Vector3D();

                for (int k=0; k<i; k++) {
                    IAtomList atoms2 = molecules[k].getChildList();
                    for (int l=0; l<atoms2.size(); l++) {
                        PotentialMoleculePairAmoebaMPole.Multipole klmpole = mpoles[k][l];
                        IAtom a2 = atoms2.get(l);

                        Vector dr = space.makeVector();
                        dr.Ev1Mv2(a2.getPosition(), a1.getPosition());
                        double r2 = dr.squared();
                        double r = Math.sqrt(r2);

                        //cos(dipole 1 and r12)
                        double c1 = ijmpole.mu.dot(dr);  // dir
                        //cos(r12 and dipole 2)
                        double c2 = dr.dot(klmpole.mu);  // dkr

                        Vector Q1dr = new Vector3D();
                        Q1dr.E(dr);
                        ijmpole.Q.transform(Q1dr);
                        // Q2dr now has components qix, qiy, qiz
                        double Q1drdr = Q1dr.dot(dr); // qir

                        Vector Q2dr = new Vector3D();
                        Q2dr.E(dr);
                        klmpole.Q.transform(Q2dr);
                        // Q2dr now has components qkx, qky, qkz
                        double Q2drdr = Q2dr.dot(dr); // qkr

                        double[] dmpik = dampthole(a1.getType().getIndex(), a2.getType().getIndex(), 7, r);
                        double rr3 = dmpik[3] / (r*r2);
                        double rr5 = 3 * dmpik[5] / (r*r2*r2);
                        double rr7 = 15 * dmpik[7] / (r*r2*r2*r2);
                        double a = rr3*klmpole.q - rr5*c2 + rr7*Q2drdr;
                        Vector fid = new Vector3D();
                        fid.Ea1Tv1(-a, dr);
                        fid.PEa1Tv1(-rr3, klmpole.mu);
                        fid.PEa1Tv1(2*rr5, Q2dr);
                        a = rr3*ijmpole.q + rr5*c1 + rr7*Q1drdr;
                        Vector fkd = new Vector3D();
                        fkd.Ea1Tv1(a, dr);
                        fkd.PEa1Tv1(-rr3, ijmpole.mu);
                        fkd.PEa1Tv1(-2*rr5, Q1dr);
                        Vector fid2 = new Vector3D();
                        fid2.Ea1Tv1(Electron.UNIT.fromSim(1), fid);
                        Vector fkd2 = new Vector3D();
                        fkd2.Ea1Tv1(Electron.UNIT.fromSim(1), fkd);

                        if (r < 0.6) {
                            if (polarizability[a1.getType().getIndex()] != 0 && fid.squared() != 0) return null;
                            if (polarizability[a2.getType().getIndex()] != 0 && fkd.squared() != 0) return null;
                            if (polarizability[a1.getType().getIndex()] * polarizability[a2.getType().getIndex()] != 0) return null;
                        }
                        if (r < 0.2) {
                            // very close but apparently not interacting???
                            System.err.println("close ! "+a1.getLeafIndex()+" "+a2.getLeafIndex()+" "+r);
                        }
                        rmin = Math.min(r, rmin);
                        //TODO incorporate dscale for larger molecules
                        field[i][j].PE(fid);
                        field[k][l].PE(fkd);
                    }

                }
            }
        }
//        System.out.println("rmin "+rmin);
        return field;
    }

    protected Vector[][] computeField2(IMolecule[] molecules, Vector[][] inducedDipoles) {
        Vector[][] field = new Vector[inducedDipoles.length][];
        for (int i=0; i<molecules.length; i++) {
            IAtomList atoms1 = molecules[i].getChildList();
            field[i] = new Vector[atoms1.size()];

            for (int j=0; j<atoms1.size(); j++) {
                IAtom a1 = atoms1.get(j);
//                System.out.println(i+" "+j+" "+ijmpole.mu+" "+inducedDipoles[i][j]+" "+mu1);
                field[i][j] = new Vector3D();

                for (int k=0; k<=i; k++) {
                    IAtomList atoms2 = molecules[k].getChildList();
                    for (int l=0; l<atoms2.size(); l++) {
                        if  (k==i && l>=j) break;
                        IAtom a2 = atoms2.get(l);

                        Vector dr = space.makeVector();
                        dr.Ev1Mv2(a2.getPosition(), a1.getPosition());
                        double r2 = dr.squared();
                        double r = Math.sqrt(r2);

                        //cos(dipole 1 and r12)
                        double c1 = inducedDipoles[i][j].dot(dr);  // dir
                        //cos(r12 and dipole 2)
                        double c2 = dr.dot(inducedDipoles[k][l]);  // dkr

                        double[] dmpik = dampthole(a1.getType().getIndex(), a2.getType().getIndex(), 5, r);
                        double rr3 = dmpik[3] / (r*r2);
                        double rr5 = 3 * dmpik[5] / (r*r2*r2);
                        double a = - rr5*c2;
                        Vector fid = new Vector3D();
                        fid.Ea1Tv1(-a, dr);
                        fid.PEa1Tv1(-rr3, inducedDipoles[k][l]);
                        a = rr5*c1;
                        Vector fkd = new Vector3D();
                        fkd.Ea1Tv1(a, dr);
                        fkd.PEa1Tv1(-rr3, inducedDipoles[i][j]);
                        Vector fid2 = new Vector3D();
                        fid2.Ea1Tv1(Electron.UNIT.fromSim(1), fid);
                        Vector fkd2 = new Vector3D();
                        fkd2.Ea1Tv1(Electron.UNIT.fromSim(1), fkd);
//                        if (a2.getLeafIndex()==0) {
//                            System.out.println("field "+a1.getLeafIndex()+" "+a2.getLeafIndex()+" "+fid2+" "+fkd2);
//                        }

                        //TODO incorporate dscale for larger molecules
                        field[i][j].PE(fid);
                        field[k][l].PE(fkd);
                    }

                }
            }
        }
        return field;
    }


    public double energy(IMolecule[] molecules) {
        PotentialMoleculePairAmoebaMPole.Multipole[][] mpoles = new PotentialMoleculePairAmoebaMPole.Multipole[molecules.length][];
        for (int i=0; i<molecules.length; i++) {
            IAtomList atoms = molecules[i].getChildList();
            mpoles[i] = new PotentialMoleculePairAmoebaMPole.Multipole[atoms.size()];

            for (IAtom a : atoms) {
                PotentialMoleculePairAmoebaMPole.Axis myAxis = axis[a.getParentGroup().getType().getIndex()][a.getIndex()];
                Tensor rot = PotentialMoleculePairAmoebaMPole.getRotmat(a, myAxis);
                int ia = a.getIndex();
                mpoles[i][ia] = new PotentialMoleculePairAmoebaMPole.Multipole();
                PotentialMoleculePairAmoebaMPole.Multipole myMultipole = multipoles[molecules[i].getType().getIndex()][ia];
                mpoles[i][ia].q = myMultipole.q;
                mpoles[i][ia].mu = new Vector3D();
                mpoles[i][ia].mu.E(myMultipole.mu);
                rot.transform(mpoles[i][ia].mu);
                mpoles[i][ia].Q = transformQ(myMultipole.Q, rot);
            }
        }

        Vector[][] uind = new Vector[molecules.length][];
        for (int i=0; i<uind.length; i++) {
            IAtomList atoms = molecules[i].getChildList();
            uind[i] = new Vector[atoms.size()];
            for (int j=0; j<uind[i].length; j++) {
                uind[i][j] = new Vector3D();
            }
        }

        double lastRes = Double.POSITIVE_INFINITY;
        Vector[][] field0 = computeField(molecules, mpoles);
        if (field0 == null) return Double.POSITIVE_INFINITY;
        for (int iter=0; ; iter++) {
            Vector[][] ifield = computeField2(molecules, uind);
            Vector f = new Vector3D();
            f.Ea1Tv1(Electron.UNIT.fromSim(1), field0[0][0]);
            f.PEa1Tv1(Electron.UNIT.fromSim(1), ifield[0][0]);
//            System.out.println("field00 "+field0[0][0]+" "+ifield[0][0]+" "+f);
//            f.Ea1Tv1(Electron.UNIT.fromSim(1), field[0][1]);
//            System.out.println("field01 "+field[0][1]+" "+f);
//            f.Ea1Tv1(Electron.UNIT.fromSim(1), field[1][0]);
//            System.out.println("field10 "+field[1][0]+" "+f);
//            f.Ea1Tv1(Electron.UNIT.fromSim(1), field[1][1]);
//            System.out.println("field11 "+field[1][1]+" "+f);
            double res = 0;
            for (int i = 0; i < uind.length; i++) {
                IAtomList atoms = molecules[i].getChildList();
                for (int j = 0; j < uind[i].length; j++) {
                    IAtom a = atoms.get(j);
                    Vector uindNew = new Vector3D();
                    uindNew.Ea1Tv1(polarizability[a.getType().getIndex()], field0[i][j]);
                    uindNew.PEa1Tv1(polarizability[a.getType().getIndex()], ifield[i][j]);
                    double b = 0.8;
                    uindNew.TE(b);
                    uindNew.PEa1Tv1(1-b, uind[i][j]);
                    res += uindNew.Mv1Squared(uind[i][j]);
                    uind[i][j].E(uindNew);

                    if (i==0 && j==0 && false) {
//                        System.out.println("induced "+polarizability[a.getType().getIndex()]+" "+field0[i][j]+" "+ifield[i][j]+" "+uind[0][0]);
                        Vector foo = new Vector3D(), bar = new Vector3D(), blah = new Vector3D();
                        foo.Ea1Tv1(Electron.UNIT.fromSim(1), field0[i][j]);
                        blah.Ea1Tv1(Electron.UNIT.fromSim(1), ifield[i][j]);
                        bar.Ea1Tv1(Electron.UNIT.fromSim(1), uind[0][0]);
//                        System.out.println("        "+foo+" "+blah+" "+bar);
                        System.out.println(iter+" "+ifield[i][j].getX(0)+" "+ifield[0][0].getX(1)+" "+ifield[0][0].getX(2));
                    }
                }
            }
//            System.out.println(iter+" "+Math.sqrt(res));
            if (res < tol || res >= lastRes) break;
            lastRes = res;
//            Vector foo = new Vector3D();
//            foo.Ea1Tv1(Debye.UNIT.fromSim(1), uind[0][0]);
//            System.out.println("induced dipole "+foo);
        }

        double u = 0;
        for (int i=0; i<molecules.length; i++) {
            for (int j=i+1; j<molecules.length; j++) {
                u += energy(molecules[i], molecules[j], mpoles[i], mpoles[j], uind[i], uind[j]);
            }
        }

        return u;
    }

    public double energy(IMolecule molecule1, IMolecule molecule2, PotentialMoleculePairAmoebaMPole.Multipole[] mpole1, PotentialMoleculePairAmoebaMPole.Multipole[] mpole2, Vector[] uind1, Vector[] uind2) {
        IAtomList atoms1 = molecule1.getChildList();
        IAtomList atoms2 = molecule2.getChildList();
        double u = 0;

        for (IAtom a1 : atoms1) {
            int i1 = a1.getIndex();
            PotentialMoleculePairAmoebaMPole.Multipole ampole1 = mpole1[i1];
            for (IAtom a2 : atoms2) {
                int i2 = a2.getIndex();
                PotentialMoleculePairAmoebaMPole.Multipole ampole2 = mpole2[i2];
                Vector dr = space.makeVector();
                dr.Ev1Mv2(a2.getPosition(), a1.getPosition());
                double r2 = dr.squared();
                double r = Math.sqrt(r2);
                double[] dmpik = dampthole(a1.getType().getIndex(), a2.getType().getIndex(), 7, r);
                double rr3 = dmpik[3] / (r*r2);
                double rr5 = 3 * dmpik[5] / (r*r2*r2);
                double rr7 = 15 * dmpik[7] / (r*r2*r2*r2);

                double mmcos_D1_u2 = ampole1.mu.dot(uind2[i2]);
                double mmcos_u1_D2 = uind1[i1].dot(ampole2.mu);
                //cos(dipole 1 and r12)
                double uir = uind1[i1].dot(dr);  // uir/r
                double dir = ampole1.mu.dot(dr);  // uir/r
                //cos(r12 and dipole 2)
                double ukr = dr.dot(uind2[i2]);  // ukr/r
                double dkr = dr.dot(ampole2.mu);  // ukr/r
                double uqmu = -ampole1.q*ukr*rr3;
                double umuq =  ampole2.q*uir*rr3;
                double uumu = mmcos_u1_D2*rr3 - uir * dkr * rr5;
                double umuu = mmcos_D1_u2*rr3 - dir * ukr * rr5;

                Vector Q1dr = new Vector3D();
                Q1dr.E(dr);
                ampole1.Q.transform(Q1dr);
                double Q1drdr = Q1dr.dot(dr); // qir
                double Q1drmu2 = Q1dr.dot(uind2[i2]); // dkqi

                Vector Q2dr = new Vector3D();
                Q2dr.E(dr);
                ampole2.Q.transform(Q2dr);
                double Q2drdr = Q2dr.dot(dr); // qkr
                double Q2drmu1 = Q2dr.dot(uind1[i1]); // diqk

                double umuQ = -2*Q2drmu1 * rr5 + uir*Q2drdr * rr7;
                double uQmu = 2*Q1drmu2 * rr5 - ukr*Q1drdr * rr7;

//                System.out.printf("mp %d %d  % 10.5e % 10.5e  % 10.5e % 10.5e  % 10.5e % 10.5e\n", a1.getLeafIndex(), a2.getLeafIndex(),
//                        0.5*kcalpmole.fromSim(umuq), 0.5*kcalpmole.fromSim(uqmu),
//                        0.5*kcalpmole.fromSim(umuu), 0.5*kcalpmole.fromSim(uumu),
//                        0.5*kcalpmole.fromSim(uQmu), 0.5*kcalpmole.fromSim(umuQ));

                u += 0.5*(uqmu + umuq + uumu + umuu+ umuQ + uQmu);
//                System.out.printf("mp %d %d  % 10.7e\n", a1.getLeafIndex(), a2.getLeafIndex(), 0.5*kcalpmole.fromSim(uqmu + umuq + uumu + umuu + umuQ + uQmu));
            }
        }
//        Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
//        System.out.println("polarization tot "+kcalpmole.fromSim(u));
        return u;
    }
}
