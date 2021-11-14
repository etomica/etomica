/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.util.random.RandomNumberGenerator;

/**
 * r^-12 potential plus a point dipole, force shifted, as used in
 * Kol, Laird and Leimkuhler, JCP 107 (1997) 2580-2588.
 *
 * @author Andrew Schultz
 */
public class P2LJDipoleAtomic implements Potential2Soft {

    private final Space space;

    public P2LJDipoleAtomic(Space space) {
        this(space, 1, 1, 1);
    }

    public P2LJDipoleAtomic(Space space, double sigma, double epsilon, double momentSquared) {
        this.space = space;
        setSigma(sigma);
        setEpsilon(epsilon);
        setDipoleMomentSquare(momentSquared);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        torque = new Vector[2];
        torque[0] = space.makeVector();
        torque[1] = space.makeVector();
        dr = space.makeVector();
        drunit = space.makeVector();
        work = space.makeVector();
        gradientAndTorque = new Vector[][]{gradient, torque};
    }

    public void setHardCoreDiamterSq(double val) {
        hsdiasq = val;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        cutoff2 = rCut * rCut;
        double r2 = rCut * rCut;
        double s2 = sigma2 / r2;
        double s6 = s2 * s2 * s2;
        double s1 = Math.sqrt(s2);
        fShift = epsilon48 * s6 * s6 * s1;
        uShift = -(epsilon4 * s6 * s6 + fShift / s1);
        dipoleFShift = 0.75 * s6 * s1;
        dipoleUShift = -1.75 * s2 * s1;
    }

    public double getRange() {
        return Math.sqrt(cutoff2);
    }

    public double energy(IAtomList pair) {
        IAtomOriented atom1 = (IAtomOriented) pair.get(0);
        IAtomOriented atom2 = (IAtomOriented) pair.get(1);

        // LJ contributation

        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / (r2);
        double s1 = Math.sqrt(s2);
        double s6 = s2 * s2 * s2;
        double ener = epsilon4 * s6 * s6 + fShift / s1 + uShift;

        if (momentSq != 0.0) {
            // v1 is the orientation of molecule 1
            Vector v1 = atom1.getOrientation().getDirection();

            // v2 is the orientation of molecule 2
            Vector v2 = atom2.getOrientation().getDirection();

            // we didn't normalize dr, so divide by r2 here
            double udd = v1.dot(v2) - 3.0 * v1.dot(dr) * v2.dot(dr) * s2;
            double fac = s2 * s1 + dipoleFShift / (s2 * s2) + dipoleUShift;
            ener += momentSq * fac * udd;
        }
        return ener;
    }

    public double getSigma() {
        return sigma;
    }

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s * s;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4 * epsilon;
        epsilon48 = 48 * epsilon;
        if (cutoff2 > 0) {
            setTruncationRadius(Math.sqrt(cutoff2));
        }
    }

    public void setDipoleMomentSquare(double newMomentSq) {
        momentSq = newMomentSq;
    }

    public double getDipoleMomentSquare() {
        return momentSq;
    }
    private double sigma, sigma2;
    private double epsilon, epsilon4, epsilon48;
    private double hsdiasq = 1.0 / Math.sqrt(2);
    private double momentSq;
    private Boundary boundary;
    private final Vector dr, drunit, work;
    private final Vector[] gradient;
    protected final Vector[] torque;
    protected final Vector[][] gradientAndTorque;
    protected double cutoff2;
    protected double dipoleFShift, dipoleUShift;
    protected double fShift, uShift;

    public static void main(String[] args) {
        RandomNumberGenerator random = new RandomNumberGenerator();
        Space3D space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesGeneral species = SpeciesSpheresRotating.create(space, new ElementSimple(sim), true, true);
        sim.addSpecies(species);
        Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
        sim.addBox(box);
        box.setNMolecules(species, 2);

        double rCut = 2.5;
        P2LJDipoleAtomic potential = new P2LJDipoleAtomic(sim.getSpace());
        potential.setTruncationRadius(rCut);
        potential.setEpsilon(1);
        potential.setDipoleMomentSquare(1);

        IAtomList leafAtoms = box.getLeafList();
//        IAtomOriented atom0 = (IAtomOriented)leafAtoms.getMolecule(0);
        IAtomOriented atom1 = (IAtomOriented) leafAtoms.get(1);

        Vector oldPosition = space.makeVector();
        Vector ran = space.makeVector();
//        atom0.getOrientation().randomRotation(sim.getRandom(), 1);
//        atom1.getOrientation().randomRotation(sim.getRandom(), 1);
//        atom1.getOrientation().randomRotation(sim.getRandom(), 1);

        atom1.getPosition().setX(0, 1);
//        for (int i=0; i<101; i++) {
//            atom1.getPosition().setX(0, 1+i*(rCut-1)*0.01);
//            double u = potential.energy(leafAtoms);
//            Vector[] gradient = potential.gradient(leafAtoms);
//            System.out.println(atom1.getPosition().get(0)+" "+u+" "+gradient[0]);
//        }
//        System.exit(1);
        for (int i = 0; i < 100; i++) {
            // calculate the gradient
            System.out.println("do gradient");
            oldPosition.E(atom1.getPosition());
            System.out.println("pos " + oldPosition + " " + Math.sqrt(oldPosition.squared()));
            // calculate the gradient numerically for 1 direction
            int d = random.nextInt(3);
            double h = 0.00001;
            System.out.println("d=" + d);
            atom1.getPosition().setX(d, oldPosition.getX(d) + h);
            double uplus = potential.energy(leafAtoms);
            System.out.println("U plus " + uplus);
            atom1.getPosition().setX(d, oldPosition.getX(d) - h);
            double uminus = potential.energy(leafAtoms);
            System.out.println("U minus " + uminus);
            double du = (uplus - uminus) / (2 * h);
            if (Double.isNaN(du)) {
                throw new RuntimeException("oops " + du + " " + uminus + " " + uplus);
            }
            System.out.println(du);
            System.out.println("success");
            do {
                // move the atom1 to an entirely random position within 5sigma of atom0
                atom1.getPosition().setRandomInSphere(random);
                atom1.getPosition().TE(rCut + 1);
                ran.setRandomSphere(random);
                atom1.getOrientation().setDirection(ran);
                System.out.println("direction = " + ran);
            }
            while (Double.isInfinite(potential.energy(leafAtoms)));
        }
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        // LJ contributation

        double r2 = dr12.squared();
        if (r2 < hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / r2;
        double s1 = Math.sqrt(s2);
        double s6 = s2 * s2 * s2;
        double ener = epsilon4 * s6 * s6 + fShift / s1 + uShift;

        if (momentSq != 0.0) {
            // v1 is the orientation of molecule 1
            Vector v1 = ((IAtomOriented) atom1).getOrientation().getDirection();

            // v2 is the orientation of molecule 2
            Vector v2 = ((IAtomOriented) atom2).getOrientation().getDirection();

            // dr12 is not normalized, so divide by r2 here
            double udd = v1.dot(v2) - 3.0 * v1.dot(dr12) * v2.dot(dr12) * s2;
            double fac = s2 * s1 + dipoleFShift / (s2 * s2) + dipoleUShift;
            ener += momentSq * fac * udd;
        }
        return ener;
    }

    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        // LJ contributation

        double r2 = dr12.squared();
        if (r2 < hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / r2;
        double s1 = Math.sqrt(s2);
        double s6 = s2 * s2 * s2;
        double ener = epsilon4 * s6 * s6 + fShift / s1 + uShift;

        // v1 is the orientation of molecule 1
        Vector v1 = ((IAtomOriented) atom1).getOrientation().getDirection();

        // v2 is the orientation of molecule 2
        Vector v2 = ((IAtomOriented) atom2).getOrientation().getDirection();

        // LJ contributation

        double rdudr = epsilon48 * s6 * s6 - fShift / s1;
        Vector ftmp1 = space.makeVector();
        ftmp1.Ea1Tv1(-rdudr / r2, dr12);

        if (momentSq != 0.0) {
            // normalize dr, the vector between the molecules
            drunit.E(dr12);
            double r12 = Math.sqrt(r2);
            drunit.TE(1 / r12);


            double fac = momentSq * (s2 * s1 + dipoleFShift / (s2 * s2) + dipoleUShift);
            double dfac = momentSq * (3 * s2 * s2 * s1 - 4 * dipoleFShift / (s2));
            double udd = v1.dot(v2) - 3.0 * v1.dot(drunit) * v2.dot(drunit);
            work.Ea1Tv1(v2.dot(drunit), v1);
            work.PEa1Tv1(v1.dot(drunit), v2);
            work.PEa1Tv1(-2 * v1.dot(drunit) * v2.dot(drunit) * s1, dr12);
            work.TE(3.0 * s1 * fac);
            ftmp1.ME(work);
            ftmp1.PEa1Tv1(-dfac * udd, dr12);

            work.E(v1);
            work.XE(v2);
            Vector ttmp = space.makeVector();
            ttmp.E(v1);
            ttmp.XE(drunit);
            ttmp.TE(3.0 * v2.dot(drunit));
            ttmp.ME(work);
            ttmp.TE(fac);
            t1.PE(ttmp);

            ttmp.E(v2);
            ttmp.XE(drunit);
            ttmp.TE(3.0 * v1.dot(drunit));
            ttmp.PE(work);
            ttmp.TE(fac);
            t2.PE(ttmp);

            ener += fac * udd;
        }

        // pairwise additive, so
        f1.PE(ftmp1);
        f2.ME(ftmp1);

        return ener;
    }
}
