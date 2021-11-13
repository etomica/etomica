/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.DipoleSourceAtomic;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Potential to implement the reaction field treatment of electrostatic interactions. In this approach, the long
 * range electrostatic interactions are modeled via a dielectric continuum with dielectric constant epsilon.
 * The net effect is that molecules within a specified cutoff distance interact with an effective pair potential
 * that is otherwise independent of their separation distance.
 * <p>
 * In conjunction with this two-body potential, a zero-body potential must also be added to the PotentialMaster.
 * The zero-body potential is provided by the method makeP0.
 */

public class P2ReactionFieldDipole implements Potential2Soft {
    protected final Vector iDipole, cavityDipole;
    protected final Vector dr;
    protected final Vector[][] gradientAndTorque;
    protected final Vector3D[] a;
    protected DipoleSourceAtomic dipoleSource;
    protected double cutoff2, cutoff;
    protected double epsilon;
    protected double fac;
    protected final Hessian h;


    /**
     * @param space              defines the dimensionality of the simulation
     */
    public P2ReactionFieldDipole(Space space, DipoleSourceAtomic dipoleSource) {
        this.dipoleSource = dipoleSource;
        iDipole = space.makeVector();
        cavityDipole = space.makeVector();
        dr = space.makeVector();
        gradientAndTorque = new Vector[2][2];
        gradientAndTorque[0][0] = space.makeVector();
        gradientAndTorque[0][1] = space.makeVector();
        gradientAndTorque[1][0] = space.makeVector();
        gradientAndTorque[1][1] = space.makeVector();
        a = new Vector3D[3];
        a[0] = (Vector3D) space.makeVector();
        a[1] = (Vector3D) space.makeVector();
        a[2] = (Vector3D) space.makeVector();
        h = new Hessian(space.makeTensor(), space.makeTensor(), space.makeTensor(),
                space.makeTensor(), space.makeTensor(), space.makeTensor());
    }

    /**
     * The DipoleSource returns the dipole vector of a given molecule.
     *
     * @return the dipole source used by this object.
     */
    public DipoleSourceAtomic getDipoleSource() {
        return dipoleSource;
    }

    /**
     * @param newDipoleSource the DipoleSource to set
     */
    public void setDipoleSource(DipoleSourceAtomic newDipoleSource) {
        dipoleSource = newDipoleSource;
    }

    /**
     * The radius of the sphere within which molecules contribute to the reaction
     * field of a given molecule. This value should match the truncation radius
     * for electrostatic interactions.
     *
     * @return the cutoff radius of the reactionField.
     */
    public double getRange() {
        return cutoff;
    }

    /**
     * @param newRange the new cutoff radius of the reactionField.
     */
    public void setRange(double newRange) {
        cutoff = newRange;
        cutoff2 = newRange * newRange;
        if (epsilon < Double.POSITIVE_INFINITY) {
            fac = 2 * (epsilon - 1) / (2 * epsilon + 1) / (cutoff2 * cutoff);
        } else {
            fac = 1 / (cutoff2 * cutoff);
        }
    }

    /**
     * @return the dielectric constant of the medium surrounding the cavity.
     */
    public double getDielectric() {
        return epsilon;
    }

    /**
     * @param newDielectric the new value of the dielectric constant of
     *                      the medium surrounding the cavity.
     */
    public void setDielectric(double newDielectric) {
        epsilon = newDielectric;
        if (cutoff > 0) {
            if (epsilon < Double.POSITIVE_INFINITY) {
                fac = 2 * (epsilon - 1) / (2 * epsilon + 1) / (cutoff2 * cutoff);
            } else {
                fac = 1 / (cutoff2 * cutoff);
            }
        }
    }

    @Override
    public int nBody() {
        return 2;
    }

    @Override
    public Hessian d2u(Vector dr12, IAtom atom1, IAtom atom2) {
        double r2 = dr12.squared();
        h.o1o1.E(0);
        h.o2o2.E(0);
        h.o1o2.E(0);
        if (r2 > cutoff * cutoff) return h;
        iDipole.E(dipoleSource.getDipole(atom1));
        Vector jDipole = Vector.d(dr.getD());
        jDipole.E(dipoleSource.getDipole(atom2));

        double exi = iDipole.getX(0);//ei and ej is the dipole orientation with mu
        double eyi = iDipole.getX(1);
        double ezi = iDipole.getX(2);
        double exj = jDipole.getX(0);
        double eyj = jDipole.getX(1);
        double ezj = jDipole.getX(2);


        a[0].E(0, ezj, -eyj);
        a[0].XE(iDipole);


        a[1].E(-ezj, 0, exj);
        a[1].XE(iDipole);

        a[2].E(eyj, -exj, 0);
        a[2].XE(iDipole);
        h.o1o2.E(a);

        h.o1o2.TE(-fac);//ij and ji are the same

        a[0].E(0, -ezi, eyi);
        a[0].XE(jDipole);

        a[1].E(ezi, 0, -exi);
        a[1].XE(jDipole);

        a[2].E(-eyi, exi, 0);
        a[2].XE(jDipole);
        h.o1o1.E(a);

        h.o1o1.TE(-fac);//ii

        a[0].E(0, -ezj, eyj);
        a[0].XE(iDipole);
        a[1].E(ezj, 0, -exj);
        a[1].XE(iDipole);
        a[2].E(-eyj, exj, 0);
        a[2].XE(iDipole);
        h.o2o2.E(a);


        h.o2o2.TE(-fac);//jj
        return h;
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
        if (r2 > cutoff * cutoff) return 0;

        iDipole.E(dipoleSource.getDipole(atom1));
        double idotj = iDipole.dot(dipoleSource.getDipole(atom2));
//        System.out.println(idotj+" "+(-fac*idotj));
        return -fac * idotj;
    }

    @Override
    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        double r2 = dr12.squared();
        if (r2 > cutoff * cutoff) {
            return 0;
        }

        iDipole.E(dipoleSource.getDipole(atom1));
        iDipole.XE(dipoleSource.getDipole(atom2));
        iDipole.TE(fac);
        t1.PE(iDipole);
        t2.PEa1Tv1(-1, iDipole);

        iDipole.E(dipoleSource.getDipole(atom1));
        double idotj = iDipole.dot(dipoleSource.getDipole(atom2));
//        System.out.println(idotj+" "+(-fac*idotj)+" data");
        return -fac * idotj;
    }

    @Override
    public double energy(IAtomList atoms) {
        return 0;
    }

    @Override
    public Vector[][] gradientAndTorque(IAtomList atoms) {
        return new Vector[0][];
    }

    @Override
    public double virial(IAtomList atoms) {
        return 0;
    }

    @Override
    public Vector[] gradient(IAtomList atoms) {
        return new Vector[0];
    }

    @Override
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return new Vector[0];
    }

    @Override
    public double u(double r2) {
        return 0;
    }
}
