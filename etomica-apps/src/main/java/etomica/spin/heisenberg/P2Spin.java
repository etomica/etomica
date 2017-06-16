/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.IPotentialTorque;
import etomica.potential.Potential2;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space1d.Tensor1D;
import etomica.space1d.Vector1D;
//TODO

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
    private final Vector[][] gradientAndTorque;
    private final Vector[] gradient;
    protected Vector dr;
    protected Vector dr2;
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
        torque[0] = new Vector1D();
        torque[1] = new Vector1D();
        secondDerivative = new Tensor[3];
        this.secondDerivative[0] = new Tensor1D();
        this.secondDerivative[1] = new Tensor1D();
        this.secondDerivative[2] = new Tensor1D();
        gradientAndTorque = new Vector[][]{gradient, torque};
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

    public Vector[][] gradientAndTorque(IAtomList atoms) {

        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);

        double x1 = atom1.getOrientation().getDirection().getX(0);//cost1
        double y1 = atom1.getOrientation().getDirection().getX(1);//sint1
        double x2 = atom2.getOrientation().getDirection().getX(0);//cost2
        double y2 = atom2.getOrientation().getDirection().getX(1);//sint2

        //u=-J*cos(t1-t2) and  du/dt1 = J*sin(t1-t2) = J*(sint1*cost2- cost1*sint2 =y1*x2-x1*y2)
        double JSin = coupling * (y1 * x2 - x1 * y2);

        torque[0].E(-JSin);
        torque[1].E(JSin);
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
        double JCos = coupling * atom1.getOrientation().getDirection().dot(atom2.getOrientation().getDirection());


        secondDerivative[0].E(-JCos);//ij
        secondDerivative[1].E(JCos);//ii
        secondDerivative[2].E(JCos);//jj

//    	System.out.println(secondDerivative[0].component(0, 0));
//    	System.out.println(secondDerivative[1].component(0, 0));
//    	System.out.println(secondDerivative[2].component(0, 0));
//    	System.out.println("test for secondDerivative in p2Spin");
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
