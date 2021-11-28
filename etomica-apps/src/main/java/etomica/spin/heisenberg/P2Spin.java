/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotential2;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space1d.Tensor1D;

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

public class P2Spin implements IPotential2 {

    protected final Tensor[] secondDerivative;
    protected Vector dr;
    private double coupling;
    private final Hessian h;

    public P2Spin() {
        this(1.0);
    }

    public P2Spin(double coupling) {
        setCoupling(coupling);
        secondDerivative = new Tensor[3];
        this.secondDerivative[0] = new Tensor1D();
        this.secondDerivative[1] = new Tensor1D();
        this.secondDerivative[2] = new Tensor1D();
        h = new Hessian(new Tensor1D(), new Tensor1D(), new Tensor1D(), new Tensor1D(), new Tensor1D(), new Tensor1D());
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        IAtomOriented atom1o = (IAtomOriented) atom1;
        IAtomOriented atom2o = (IAtomOriented) atom2;
        return -coupling * atom1o.getOrientation().getDirection().dot(atom2o.getOrientation().getDirection());
    }

    @Override
    public double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
        IAtomOriented atom1o = (IAtomOriented) atom1;
        IAtomOriented atom2o = (IAtomOriented) atom2;
        return -coupling * atom1o.getOrientation().getDirection().dot(atom2o.getOrientation().getDirection());
    }

    @Override
    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        IAtomOriented atom1o = (IAtomOriented) atom1;
        IAtomOriented atom2o = (IAtomOriented) atom2;

        double x1 = atom1o.getOrientation().getDirection().getX(0);//cost1
        double y1 = atom1o.getOrientation().getDirection().getX(1);//sint1
        double x2 = atom2o.getOrientation().getDirection().getX(0);//cost2
        double y2 = atom2o.getOrientation().getDirection().getX(1);//sint2

        //u=-J*cos(t1-t2) and  du/dt1 = J*sin(t1-t2) = J*(sint1*cost2- cost1*sint2 =y1*x2-x1*y2)
        double JSin = coupling * (y1 * x2 - x1 * y2);

        t1.PE(-JSin);
        t2.PE(JSin);

        return -coupling * atom1o.getOrientation().getDirection().dot(atom2o.getOrientation().getDirection());
    }


    /**
     * Returns 0, because potential operates on a lattice and range
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
     * compute the secondDerivative array of pair energy w.r.t theta1 or theta2
     * i.e d^2u/dtheta1_dtheta1 d^2u/dtheta1_dtheta2 and d^2u/dtheta2_dtheta2
     * theta1 is the angle between x axis and atom1's orientation etc.
     *
     * @param atoms given pair of atoms
     * @return secondDerivative array
     */
    public Tensor[] secondDerivative(IAtomList atoms) {
        IAtomOriented atom1 = (IAtomOriented) atoms.get(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.get(1);
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

    @Override
    public Hessian d2u(Vector dr12, IAtom atom1, IAtom atom2) {
        double JCos = coupling * ((IAtomOriented)atom1).getOrientation().getDirection()
                            .dot(((IAtomOriented)atom2).getOrientation().getDirection());
        h.o1o2.E(-JCos);
        h.o1o1.E(JCos);
        h.o2o2.E(JCos);
        return h;
    }


}
