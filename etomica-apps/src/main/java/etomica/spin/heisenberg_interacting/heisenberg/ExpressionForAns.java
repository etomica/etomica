package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;
import etomica.space.Vector;

import static etomica.math.SpecialFunctions.besselI;

public class ExpressionForAns {
    protected Vector ei, ej;
    protected final double mu, J, bt, bJ, bmu; //TODO should I add final here
    protected double Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;

    public ExpressionForAns(Space space, double dipoleMagnitude, double interactionS, double beta) {
        ei = space.makeVector();
        ej = space.makeVector();
        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
    }

    public void doCalculation(IAtomList atoms, int nMax) {
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());
        double t1 = Math.acos(ei.getX(0));
        double t2 = Math.acos(ej.getX(0));
        double I0bJ = besselI(0, bJ);
        double I1bJ = besselI(1, bJ);
        double I2bJ = besselI(2, bJ);

        double test0 = 0;
        double test1 = 0;
        double test2 = 0;
        double test3 = 0;

        Axc0 = bmu * (I0bJ + I1bJ) * (-1 + Math.cos(t1));
        Axs0 = 0;
        dAxc0 = -bmu * (I0bJ + I1bJ) * Math.sin(t1);
        dAxs0 = 0;
        Axc1 = -0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);

        test0 = Axc1;

        if (nMax > 0) {
            for (int n = 1; n <= nMax; n++) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double InbJ = besselI(n, bJ);
                double Inm1bJ = besselI(n - 1, bJ);
                double Inm2bJ = besselI(n - 2, bJ);
                double Inp1bJ = besselI(n+1,bJ);
                Axc0 += 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                        + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1))
                        / (bJ + 4 * bJ * n4);

                Axs0 += (-2 * bmu * n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(t1) * Math.cos(n * t1)
                        - 2 * bmu * (bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.sin(t1) * Math.sin(n * t1)
                ) / (bJ + 4 * bJ * n4);

                dAxc0 += (2 * InbJ * (-Math.cos(n * t1) * Math.sin(t1) + (n - 2 * n3) * Math.cos(t1) * Math.sin(n * t1))
                        + Inm1bJ * (1 + n - 2 * n3) * Math.sin((n - 1) * t1)
                        + (-1+n-2*n3)*Inp1bJ*Math.sin((n+1)*t1)
                ) * bmu / (1 + 4 * n4);

                dAxs0 += ((1+2*n+2*n2)*Inm1bJ*Math.sin((n-1)*t1)
                        +2*InbJ*(-2*n*Math.cos(n*t1)*Math.sin(t1)+(1+2*n2)*Math.cos(t1)*Math.sin(n*t1))
                        +(1-2*n+2*n2)*Inp1bJ*Math.sin((n+1)*t1)
                        )*bmu/(1+4*n4);

                Axc1 += bmu * bmu * (
                        Inm2bJ * (Math.cos((n - 2) * t1) - Math.cosh(n * t1)) / (2 - 2 * n - n2)
                                +
                                (-4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                                        + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (n - 1) * Math.cosh(n * t1)))
                                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) / 4 / n2;

                test1 += Inm2bJ * (Math.cos((n - 2) * t1) - Math.cosh(n * t1)) * n2 / (2 - 2 * n - n2);

                test3 += Inm2bJ;
                System.out.println("n=" + n + "Inm2bJ= " + Inm2bJ);//TODO besselfunction not performing well when n=1 here!!!!
                test2 += bmu * bmu * (
                        (-4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                                + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                                + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (n - 1) * Math.cosh(n * t1)))
                        ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) / 4 / n2;

            }
        }
        System.out.println("~~~~~~~~~~~~~~~Debug only ~~~~~~~~~~~~~~~~~~~~~");
        System.out.println("nMax= " + nMax);
        System.out.println("bJ= " + bJ);
        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "}");
        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "}");
//        System.out.println("Axc00= " + bmu * (I0bJ + I1bJ) * (-1 + Math.cos(t1)));
//        System.out.println("Axc0= " + Axc0);
//        System.out.println("Axs0= " + Axs0);
//        System.out.println("dAxc0= " + (-bmu * (I0bJ + I1bJ) * Math.sin(t1)));
//        System.out.println("dAxc0= " + dAxc0);
//        System.out.println("dAxs0= " + dAxs0);
        System.out.println("Axc1= " + Axc1);
//        System.out.println("test0= " + test0);
        System.out.println("test1= " + test1);
//        System.out.println("test2= " + test2);
        System.out.println("test3= " + test3);

        //Axc0 passed tge test
//        int n = 2;
//        int n2 = n * n;
//        int n3 = n2 * n;
//        int n4 = n2 * n2;
//        double InbJ = besselI(n, bJ);
//        double Inm1bJ = besselI(n - 1, bJ);
//        double Axc0n = 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
//                + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1))
//                / (bJ + 4 * bJ * n4);
//        System.out.println("part1= " + (2 * bmu / (bJ + 4 * bJ * n4)));
//        System.out.println("part2= " + (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
//                + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1)));
//        System.out.println("part21= " + ((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1));
//        System.out.println("part22= " + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1));
//        System.out.println("Axc0n= " + Axc0n);

        System.exit(2);

    }

    public void zeroSum() {
        Axc0 = 0;
    }

    public double getSum() {
        return Axc0;
    }
}
