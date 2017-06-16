package etomica.spin.heisenberg3D;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCalculationPhiSum implements PotentialCalculationMolecular {
    protected final Vector ei, ej;
    protected final Vector[] a;
    protected final Tensor iT;
    protected Vector fieldE;
    protected Vector Ai;
    protected Vector Aj;
    protected Vector dr;
    protected double secondDerivativeSum = 0;
    protected DipoleSource dipoleSource;

    public PotentialCalculationPhiSum(Space space) {
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
        double[] xD = {1, 0, 0};
        double[] yD = {0, 1, 0};
        double[] zD = {0, 0, 1};
        a[0].E(xD);
        a[1].E(yD);
        a[2].E(zD);
        iT.E(a);
    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }
        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;

        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());


        double traceij = t[0].trace();
        double traceii = t[1].trace();
        double tracejj = t[2].trace();


        t[0].transpose();
        t[0].TE(-1);
        t[1].transpose();
        t[1].TE(-1);
        t[2].transpose();
        t[2].TE(-1);


        t[0].PEa1Tt1(traceij, iT);
        t[1].PEa1Tt1(traceii, iT);
        t[2].PEa1Tt1(tracejj, iT);


        dr.E(ej);
        t[0].transform(dr);
        secondDerivativeSum += 2 * ei.dot(dr);//ij

        dr.E(ei);
        t[1].transform(dr);
        secondDerivativeSum += ei.dot(dr);//ii

        dr.E(ej);
        t[2].transform(dr);
        secondDerivativeSum += ej.dot(dr);//jj

    }

    public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
        if (!(potential instanceof IPotentialMolecularSecondDerivative)) {
            return;
        }

        IPotentialMolecularSecondDerivative potentialSecondDerivative = (IPotentialMolecularSecondDerivative) potential;

        Tensor[] t = potentialSecondDerivative.secondDerivative(molecules);

        IMolecule molecule0 = molecules.getMolecule(0);
        IMolecule molecule1 = molecules.getMolecule(1);
        ei.E(dipoleSource.getDipole(molecule0));
        ej.E(dipoleSource.getDipole(molecule1));
        ei.normalize();
        ej.normalize();

        double traceij = t[0].trace();
        double traceii = t[1].trace();
        double tracejj = t[2].trace();


        t[0].transpose();
        t[0].TE(-1);
        t[1].transpose();
        t[1].TE(-1);
        t[2].transpose();
        t[2].TE(-1);


        t[0].PEa1Tt1(traceij, iT);
        t[1].PEa1Tt1(traceii, iT);
        t[2].PEa1Tt1(tracejj, iT);


        dr.E(ej);
        t[0].transform(dr);
        secondDerivativeSum += 2 * ei.dot(dr);//ij


        dr.E(ei);
        t[1].transform(dr);

        secondDerivativeSum += ei.dot(dr);//ii

        dr.E(ej);
        t[2].transform(dr);

        secondDerivativeSum += ej.dot(dr);//jj


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
