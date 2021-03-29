/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.molecule.DipoleSourceAtomic;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCallbackPhiSumFasterer implements PotentialCallback {
    protected Vector fieldE;
    protected Vector ei, ej;
    protected Vector Ai;
    protected Vector Aj;
    protected Vector dr;
    protected double secondDerivativeSum= 0;
    protected DipoleSourceAtomic dipoleSource;
    protected final Vector[] a;
    protected final Tensor iT;

    public PotentialCallbackPhiSumFasterer(Space space) {
        fieldE = space.makeVector();
        Ai = space.makeVector();
        Aj = space.makeVector();
        dr = space.makeVector();
        ei = space.makeVector();
        ej = space.makeVector();
        a = new Vector[3];
        a[0] = space.makeVector();
        a[1] = space.makeVector();
        a[2] = space.makeVector();
        iT = space.makeTensor();
        double [] xD = {1,0,0};
        double [] yD = {0,1,0};
        double [] zD = {0,0,1};
        a[0].E(xD);
        a[1].E(yD);
        a[2].E(zD);
        iT.E(a);
    }

    @Override
    public void pairComputeGeneral(Potential2Soft pij, IAtom atom1, IAtom atom2, Vector drij, Vector fij, Vector tij, Vector tji) {

        IPotentialPair.Hessian h = pij.d2u(drij, atom1, atom2);

        ei.E(dipoleSource.getDipole(atom1));
        ej.E(dipoleSource.getDipole(atom2));
        ei.normalize();
        ej.normalize();

        double traceij = h.o1o2.trace();
        double traceii = h.o1o1.trace();
        double tracejj = h.o2o2.trace();


        h.o1o2.transpose();
        h.o1o2.TE(-1);
        h.o1o1.transpose();
        h.o1o1.TE(-1);
        h.o2o2.transpose();
        h.o2o2.TE(-1);


        h.o1o2.PEa1Tt1(traceij, iT);
        h.o1o1.PEa1Tt1(traceii, iT);
        h.o2o2.PEa1Tt1(tracejj, iT);

        dr.E(ej);
        h.o1o2.transform(dr);
        secondDerivativeSum += 2*ei.dot(dr);//ij

        dr.E(ei);
        h.o1o1.transform(dr);

        secondDerivativeSum += ei.dot(dr);//ii

        dr.E(ej);
        h.o2o2.transform(dr);

        secondDerivativeSum += ej.dot(dr);//jj
    }

    public void setDipoleSource(DipoleSourceAtomic newDipoleSource) {
        dipoleSource = newDipoleSource;
    }

    public void zeroSum() {
        secondDerivativeSum = 0.0;
    }

    /**
     * Returns the current value of the energy sum.
     */
    public double getSum() {
        return secondDerivativeSum;
    }


}
