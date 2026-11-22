/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.heatcapacity;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.Vector;

/**
 * Splits the kinetic energy of the box into center-of-mass (translation), rotation and vibration.
 * Rotation is taken with respect to the instantaneous inertia tensor of each molecule, and
 * vibration is whatever internal kinetic energy is not rotation.
 */
public class MeterModeEnergies {

    private final Box box;
    private final Vector[] dr;
    private final Vector drTmp;
    private final Vector vCOM, rCOM, dv;

    /** total kinetic energy in each mode, as of the last call to compute */
    public double keTrans, keRot, keVib;

    public MeterModeEnergies(Box box) {
        this.box = box;
        dr = new Vector[3];
        for (int i = 0; i < dr.length; i++) dr[i] = box.getSpace().makeVector();
        drTmp = box.getSpace().makeVector();
        vCOM = box.getSpace().makeVector();
        rCOM = box.getSpace().makeVector();
        dv = box.getSpace().makeVector();
    }

    public void compute() {
        keTrans = keRot = keVib = 0;
        Boundary boundary = box.getBoundary();
        for (IMolecule molecule : box.getMoleculeList()) {
            IAtomList atoms = molecule.getChildList();
            int n = atoms.size();
            double mTot = 0, keAtoms = 0;
            vCOM.E(0);
            for (int i = 0; i < n; i++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(i);
                double m = a.getType().getMass();
                mTot += m;
                vCOM.PEa1Tv1(m, a.getVelocity());
                keAtoms += 0.5 * m * a.getVelocity().squared();
            }
            vCOM.TE(1.0 / mTot);
            double keCOM = 0.5 * mTot * vCOM.squared();
            keTrans += keCOM;
            if (n == 1) continue;
            double keInternal = keAtoms - keCOM;

            // atoms might be in different images
            Vector r0 = atoms.get(0).getPosition();
            for (int i = 0; i < n; i++) {
                dr[i].Ev1Mv2(atoms.get(i).getPosition(), r0);
                boundary.nearestImage(dr[i]);
            }

            double rot = n == 2 ? rotationalDiatomic(atoms) : rotationalPolyatomic(atoms, mTot);
            keRot += rot;
            keVib += keInternal - rot;
        }
    }

    /** kinetic energy of the relative velocity perpendicular to the bond */
    private double rotationalDiatomic(IAtomList atoms) {
        IAtomKinetic a0 = (IAtomKinetic) atoms.get(0);
        IAtomKinetic a1 = (IAtomKinetic) atoms.get(1);
        double m0 = a0.getType().getMass(), m1 = a1.getType().getMass();
        double mu = m0 * m1 / (m0 + m1);
        dv.Ev1Mv2(a1.getVelocity(), a0.getVelocity());
        double vDotR = dv.dot(dr[1]);
        return 0.5 * mu * (dv.squared() - vDotR * vDotR / dr[1].squared());
    }

    /** (1/2) L . I^-1 . L, with L and I about the center of mass */
    private double rotationalPolyatomic(IAtomList atoms, double mTot) {
        int n = atoms.size();
        rCOM.E(0);
        for (int i = 0; i < n; i++) {
            rCOM.PEa1Tv1(atoms.get(i).getType().getMass(), dr[i]);
        }
        rCOM.TE(1.0 / mTot);
        double Lx = 0, Ly = 0, Lz = 0;
        double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
        for (int i = 0; i < n; i++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(i);
            double m = a.getType().getMass();
            drTmp.Ev1Mv2(dr[i], rCOM);
            dv.Ev1Mv2(a.getVelocity(), vCOM);
            double x = drTmp.getX(0), y = drTmp.getX(1), z = drTmp.getX(2);
            double vx = dv.getX(0), vy = dv.getX(1), vz = dv.getX(2);
            Lx += m * (y * vz - z * vy);
            Ly += m * (z * vx - x * vz);
            Lz += m * (x * vy - y * vx);
            Ixx += m * (y * y + z * z);
            Iyy += m * (x * x + z * z);
            Izz += m * (x * x + y * y);
            Ixy -= m * x * y;
            Ixz -= m * x * z;
            Iyz -= m * y * z;
        }
        // cofactors of the symmetric inertia tensor
        double cxx = Iyy * Izz - Iyz * Iyz;
        double cxy = Ixz * Iyz - Ixy * Izz;
        double cxz = Ixy * Iyz - Ixz * Iyy;
        double cyy = Ixx * Izz - Ixz * Ixz;
        double cyz = Ixy * Ixz - Ixx * Iyz;
        double czz = Ixx * Iyy - Ixy * Ixy;
        double det = Ixx * cxx + Ixy * cxy + Ixz * cxz;
        if (det <= 1e-12) return 0; // molecule is collinear; the bend potential makes this very rare
        double LIL = Lx * Lx * cxx + Ly * Ly * cyy + Lz * Lz * czz
                + 2 * (Lx * Ly * cxy + Lx * Lz * cxz + Ly * Lz * cyz);
        return 0.5 * LIL / det;
    }
}
