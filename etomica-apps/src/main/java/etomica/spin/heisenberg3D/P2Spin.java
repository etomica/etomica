/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg3D;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.IPotentialTorque;
import etomica.potential.Potential2;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Magnetic spin potential, with an energy defined by
 * <p>
 * U = -J r1 dot r2
 * <p>
 * where J is a coupling parameter, and r1 and r2 are the vectors given by
 * atom.coord.position. It is expected (but not verified here) that these
 * vectors are normalized to unity, and that the simulation integrator's
 * algorithm enforces this constraint.
 *
 * @author weisong lin and David Kofke
 */

public class P2Spin extends Potential2 implements IPotentialTorque, IPotentialAtomicSecondDerivative {

    private static final long serialVersionUID = 1L;
    protected final Vector[] torque;
    protected final Tensor[] secondDerivative;
    protected final Vector[] a;
    private final Vector[][] gradientAndTorque;
    private final Vector[] gradient;
    protected Vector dr;
    private double coupling;

    public P2Spin(Space space) {
        this(space, 1.0);
    }

    public P2Spin(Space space, double coupling) {
        super(space);
        setCoupling(coupling);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        torque = new Vector[2];
        torque[0] = space.makeVector();
        torque[1] = space.makeVector();
        secondDerivative = new Tensor[3];
        this.secondDerivative[0] = space.makeTensor();
        this.secondDerivative[1] = space.makeTensor();
        this.secondDerivative[2] = space.makeTensor();
        gradientAndTorque = new Vector[][]{gradient, torque};
        a = new Vector[3];
        a[0] = space.makeVector();
        a[1] = space.makeVector();
        a[2] = space.makeVector();
        dr = space.makeVector();
    }

    /**
     * Returns the energy for the given pair of atoms.
     *
     * @param atoms
     * @throws ClassCastException if atoms is not an instance of AtomPair
     */

    public double energy(IAtomList atoms) {
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        return -coupling * atom1.getOrientation().getDirection().dot(atom2.getOrientation().getDirection());
    }

    /**
     * Returns 0, becuase potential operates on a lattice and range
     * should not be needed.  The PotentialMasterSite expects all Potentials
     * to have a range and uses the return value to determine whether or not
     * to use site iteration.
     */
    public double getRange() {
        return 0;
    }

    /**
     * * @return J the coupling parameter
     */
    public double getCoupling() {
        return coupling;
    }

    /**
     * set the coupling parameter J
     *
     * @param coupling
     */
    public void setCoupling(double coupling) {
        this.coupling = coupling;
    }

    /**
     * does nothing
     *
     * @param box
     */
    public void setBox(Box box) {

    }

    /**
     * no virial is use here
     *
     * @throws Exception when virial is used
     */
    public double virial(IAtomList atoms) {

        throw new RuntimeException("virial is not used in p2Spin");

    }

    /**
     * @param atoms
     * @return gradient and torque of given pair of atoms
     */

    public Vector[][] gradientAndTorque(IAtomList atoms) {//{0,0,torque}

        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);

        Vector ei = atom1.getOrientation().getDirection();
        Vector ej = atom2.getOrientation().getDirection();

        //torque_ij = j ei cross ej
        torque[0].E(ei);
        torque[0].XE(ej);
        torque[0].TE(coupling);
        torque[1].Ea1Tv1(-1, torque[0]);
//		System.out.println("torque1 = "  + torque[0]);
        return gradientAndTorque;
    }

    /**
     * do nothing
     */
    public Vector[][] gradientAndTorque(IAtomList atoms, Tensor pressureTensor) {
        return gradientAndTorque(atoms);
    }

    /**
     * compute the secondDerivative array of pair energy w.r.t theta1 or theta2
     * i.e d^2u/dtheta1_dtheta1 d^2u/dtheta1_dtheta2 and d^2u/dtheta2_dtheta2
     * theta1 is the angle between x axis and atom1's orientation etc.
     *
     * @param atoms given pair of atoms
     * @return secondDerivative array
     */
    public Tensor[] secondDerivative(IAtomList atoms) {
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);


        Vector ei = atom1.getOrientation().getDirection();
        Vector ej = atom2.getOrientation().getDirection();


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
        double[] dejdxjD = {0, -ezj, eyj};
        double[] dejdyjD = {ezj, 0, -exj};
        double[] dejdzjD = {-eyj, exj, 0};
        dejdxj.E(dejdxjD);
        dejdyj.E(dejdyjD);
        dejdzj.E(dejdzjD);
        secondDerivative[0].E(0);//ij
        secondDerivative[1].E(0);//ii
        secondDerivative[2].E(0);//jj

        //dTau/dxj= coupling* ei cross dejdxj with xj the angle rotate around x axis
        a[0].E(ei);
        a[0].XE(dejdxj);

        a[1].E(ei);
        a[1].XE(dejdyj);

        a[2].E(ei);
        a[2].XE(dejdzj);
        secondDerivative[0].E(a);
        secondDerivative[0].TE(coupling);//ij

        double[] deidxiD = {0, -ezi, eyi};
        double[] deidyiD = {ezi, 0, -exi};
        double[] deidziD = {-eyi, exi, 0};
        deidxi.E(deidxiD);
        deidyi.E(deidyiD);
        deidzi.E(deidziD);
        //dtau/dxi = J*deidxi cross ej
        a[0].E(deidxi);
        a[0].XE(ej);

        a[1].E(deidyi);
        a[1].XE(ej);

        a[2].E(deidzi);
        a[2].XE(ej);
        secondDerivative[1].E(a);
        secondDerivative[1].TE(coupling);//ii


        //dtaujidxj = J*dejdxj cross ei
        a[0].E(dejdxj);
        a[0].XE(ei);

        a[1].E(dejdyj);
        a[1].XE(ei);

        a[2].E(dejdzj);
        a[2].XE(ei);
        secondDerivative[2].E(a);
        secondDerivative[2].TE(coupling);//jj

//		System.out.println("out=====================debug in p2Spin3D==========================");
//		System.out.println("=====================compare with mathematica==========================");
//		System.out.println("ei = " + ei);
//		//TODO I notice that one component of ei or ej would be close to one!!!!
//		System.out.println("ej = " + ej);
//		System.out.println("Tij = "	+ secondDerivative[0]);
//		System.out.println("Tii = "	+ secondDerivative[1]);
//		System.out.println("Tjj = "	+ secondDerivative[2]);
//		System.exit(2);

        return secondDerivative;
    }

    /**
     * do nothing
     */
    public Vector[] gradient(IAtomList atoms) {
        throw new RuntimeException("don't need to use gradient");
    }

    /**
     * do nothing
     */

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        throw new RuntimeException("don't need to use gradient");
    }


}
