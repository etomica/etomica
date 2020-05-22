package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.PotentialCalculation;
import etomica.space.Tensor;
import etomica.space.Vector;

import static etomica.math.SpecialFunctions.besselI;

public class PotentialCalculationVSum implements PotentialCalculation {
    protected AtomLeafAgentManager<MeterMappedAveragingSum.MoleculeAgent> leafAgentManager;

    protected Vector ei, ej,vE,vEE,vDotGradV;
    protected  double mu, J, bt, bJ, bmu;
    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    protected double psix1, psix2, psix11, psix12, psix22, psi1x1, psi1x2, psi1x11, psi1x22, psix111, psix222, psix221, psix112;
    protected double psiy1, psiy2, psiy11, psiy12, psiy22, psi1y1, psi1y2, psi1y11, psi1y22, psiy111, psiy222, psiy221, psiy112;
    protected int nMax = 10;


    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }

        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.get(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.get(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());
//        System.out.println(ei);
        double t1 = Math.atan2(ei.getX(1), ei.getX(0));
        double t2 = Math.atan2(ej.getX(1), ej.getX(0));

        double cost1 = Math.cos(t1);
        double sint1 = Math.sin(t1);
        double sint1p2 = sint1 * sint1;
        double bJ2 = bJ * bJ;
        double bmu2 = bmu * bmu;
        double I0bJ = besselI(0, bJ);
        double I1bJ = besselI(1, bJ);
        double I2bJ = besselI(2, bJ);


        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                Axc0[n] = bmu * (I0bJ + I1bJ) * (-1 + Math.cos(t1));
                dAxc0[n] = -bmu * (I0bJ + I1bJ) * Math.sin(t1);
                Axs0[n] = 0;
                dAxs0[n] = 0;
                Axc1[n] = -0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);
                dAxc1[n] = -0.25 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * Math.sin(2 * t1);
                Axs1[n] = 0;
                dAxs1[n] = 0;

                d2Axc0[n] = -bmu * Math.cos(t1) * (I0bJ + I1bJ);
                d2Axs0[n] = 0;
                d3Axc0[n] = bmu * Math.sin(t1) * (I0bJ + I1bJ);
                d3Axs0[n] = 0;
                d2Axc1[n] = -0.5 * bmu * bmu * Math.cos(2 * t1) * (I0bJ + 2 * I1bJ + I2bJ);
                d2Axs1[n] = 0;
            }


            //x direction
            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double np1p3 = (n + 1) * (n + 1) * (n + 1);
                double np1p2 = (n + 1) * (n + 1);
                double nm2p2 = (n - 2) * (n - 2);
                double InbJ = besselI(n, bJ);
                double Inm1bJ = besselI(n - 1, bJ);//TODO
                double Inm2bJ = besselI(n - 2, bJ);
                double Inp1bJ = besselI(n + 1, bJ);
                double sinnt1 = Math.sin(n * t1);
                double cosnt1 = Math.cos(n * t2);
                double sinnm1t1 = Math.sin((n - 1) * t1);
                double sinnp1t1 = Math.sin((n + 1) * t1);
                double coshnt1 = Math.cosh(n * t1);

                Axc0[n] = 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                        + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1))
                        / (bJ + 4 * bJ * n4);

                dAxc0[n] = (2 * InbJ * (-Math.cos(n * t1) * Math.sin(t1) + (n - 2 * n3) * Math.cos(t1) * Math.sin(n * t1))
                        + Inm1bJ * (1 + n - 2 * n3) * Math.sin((n - 1) * t1)
                        + (-1 + n - 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1 + 4 * n4);

                Axs0[n] = ((1 + 2 * n + 2 * n2) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (-2 * n * Math.cos(n * t1) * Math.sin(t1) + (1 + 2 * n2) * Math.cos(t1) * Math.sin(n * t1))
                        + (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1 + 4 * n4);

                dAxs0[n] = (-2 * bmu * n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(t1) * Math.cos(n * t1)
                        - 2 * bmu * (bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.sin(t1) * Math.sin(n * t1)
                ) / (bJ + 4 * bJ * n4);

                Axc1[n] = (Inm2bJ * (Math.cos((n - 2) * t1) - Math.cosh(n * t1)) * n2 / (2 - 2 * n + n2)
                        + (-4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (n - 1) * Math.cosh(n * t1)))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu * bmu / 4.0 / n2;

                dAxc1[n] = 0.25 * bmu * bmu * (
                        -(n - 2) * Inm2bJ * Math.sin((n - 2) * t1) / (2 - 2 * n + n2)
                                +
                                (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                                        - 2 * n * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (bJ - 1 - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)

                                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Axs1[n] = ((n2 * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n2 * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) - (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * Math.sin((n + 2) * t1) * ((-bJ + bJ * bJ - n * bJ) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n2 + 2 * n) * InbJ)
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu * bmu / 4.0 / n2;

                dAxs1[n] = 0.25 * bmu * bmu * (((n - 2) * Inm2bJ * Math.cos((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                        - 4 * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Axc0[n] = (-(1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.cos((n + 1) * t1)
                        + (-1 + n) * Inm1bJ * ((1 + n) * Math.cos((n - 1) * t1) - 2 * n3 * Math.cos((n - 1) * t1))
                        + InbJ * (-2 * (1 - n2 + 2 * n4) * Math.cos(t1) * Math.cos(n * t1) + 4 * n3 * Math.sin(t1) * Math.sin(n * t1))
                ) * bmu / (1 + 4 * n4);

                d2Axs0[n] = (-2 * bJ * bmu * Inm1bJ * (2 * n3 * Math.cos(n * t1) * Math.sin(t1) + (1 - n2 + 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                        + 2 * bmu * InbJ * (n * (1 - n2 - 2 * bJ * n2 + 2 * n3 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + (n * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 2 * n4)) * Math.cos(t1) * Math.sin(n * t1))
                ) / (bJ + 4 * bJ * n4);

                d3Axc0[n] = (2 * InbJ * ((1 - n2 + 4 * n4) * Math.cos(n * t1) * Math.sin(t1) + n * (1 + n2 + 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                        - (-1 + n) * (-1 + n) * Inm1bJ * ((1 + n) * Math.sin((n - 1) * t1) + 2 * n3 * Math.sin((1 - n) * t1))
                ) * bmu / (1.0 + 4 * n4);

                d3Axs0[n] = -(bJ * Inm1bJ * (n * (1 + n2 + 2 * n4) * Math.cos(t1) * Math.cos(n * t1) + (-1 + n2 - 4 * n4) * Math.sin(t1) * Math.sin(n * t1))
                        + InbJ * (-n * Math.cos(t1) * Math.cos(n * t1) * ((1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) - bJ * (1 + n2 + 2 * n4)) + Math.sin(t1) * Math.sin(n * t1) * (n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 4 * n4)))
                ) * 2 * bmu / (bJ + 4 * bJ * n4);

                d2Axc1[n] = 0.25 * bmu * bmu * (-(n - 2) * (n - 2) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                d2Axs1[n] = 0.25 * bmu * bmu * (-((n - 2) * (n - 2) * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (n + 2) * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

            }
        }


        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                Ayc0[n] = bmu * (I0bJ + I1bJ) * Math.sin(t1);
                dAyc0[n] = bmu * (I0bJ + I1bJ) * Math.cos(t1);
                Ays0[n] = 0;
                dAys0[n] = 0;
                Ayc1[n] = 0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);
                dAyc1[n] = 0.25 * bmu * bmu * Math.sin(2 * t1) * (I0bJ + 2 * I1bJ + I2bJ);
                Ays1[n] = 0;
                dAys1[n] = 0;

                d2Ayc0[n] = -bmu * Math.sin(t1) * (I0bJ + I1bJ);
                d2Ays0[n] = 0;
                d3Ayc0[n] = -bmu * Math.cos(t1) * (I0bJ + I1bJ);
                d3Ays0[n] = 0;
                d2Ayc1[n] = 0.5 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * Math.cos(2 * t1);
                d2Ays1[n] = 0;
            }


            //y direction
            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double InbJ = besselI(n, bJ);
                double Inm1bJ = besselI(n - 1, bJ);
                double Inm2bJ = besselI(n - 2, bJ);
                double Inp1bJ = besselI(n + 1, bJ);
                Ayc0[n] = (2 * bmu * Math.cos(n * t1) * Math.sin(t1) * ((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * n2 + 2 * bJ * n2 - 2 * n3) * InbJ)
                        - 2 * bmu * n * Math.cos(t1) * Math.sin(n * t1) * (2 * bJ * Inm1bJ + (1 + 2 * bJ + 2 * (-1 + n) * n) * InbJ)
                ) / (bJ + 4 * bJ * n4);

                dAyc0[n] = 2.0 * bmu * ((bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                        + n * Math.sin(t1) * Math.sin(n * t1) * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (n - bJ) * InbJ))
                ) / (bJ + 4 * bJ * n4);

                Ays0[n] = (bJ * (1 - 2 * n + 2 * n2) * Inp1bJ * (-Math.cos((n + 1) * t1) + Math.cosh(n * t1))
                        + bJ * Inm1bJ * ((1 + 2 * n + 2 * n2) * Math.cos((n - 1) * t1) + (-1 + 2 * n - 2 * n2) * Math.cosh(n * t1))
                        + 2 * InbJ * (2 * bJ * n * Math.cos(t1) * Math.cos(n * t1) + n * (1 - 2 * n + 2 * n2) * Math.cosh(n * t1) + bJ * (1 + 2 * n2) * Math.sin(t1) * Math.sin(n * t1))
                ) * bmu / (bJ + 4 * bJ * n4);

                dAys0[n] = ((1 + n - 2 * n3) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (n * (-1 + 2 * n2) * Math.cos(n * t1) * Math.sin(t1) + Math.cos(t1) * Math.sin(n * t1))
                        + (1 - n + 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1.0 + 4 * n4);


                Ayc1[n] = (n2 * Inm2bJ * (-Math.cos((n - 2) * t1) + Math.cosh(n * t1)) / (2 - 2 * n + n2)
                        - (4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (-1 + n) * Math.cosh(n * t1)))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu * bmu / 4 / n2;

                dAyc1[n] = 0.25 * bmu * bmu * (((-2 + n) * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Ays1[n] = (-bJ * n2 * (2 + n * (2 + n)) * I0bJ * Math.sin((n - 2) * t1) * ((-1 + bJ + n) * Inm1bJ + bJ * InbJ)
                        - 2 * bJ * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        - n2 * (2 - 2 * n + n2) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n * (1 + n)) * InbJ)
                ) * bmu * bmu / (2 * bJ * bJ * n2 * (4 + n4) * I0bJ);

                dAys1[n] = 0.25 * bmu * bmu * (-(-2 + n) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ) - 4 * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                d2Ayc0[n] = 2 * bmu * (-(bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * Math.cos(n * t1) * Math.sin(t1) + n2 * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(n * t1) * Math.sin(t1) - n * (bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * Math.cos(t1) * Math.sin(n * t1) + n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(t1) * Math.sin(n * t1)) / (bJ + 4 * bJ * n4);

                d2Ays0[n] = (bmu * ((-1 + n) * (1 + n - 2 * n3) * Inm1bJ * Math.cos((-1 + n) * t1) + (1 + n) * (1 - n + 2 * n3) * Inp1bJ * Math.cos(t1 + n * t1) + 2 * InbJ * (n * Math.cos(t1) * Math.cos(n * t1) + n * (-1 + 2 * n2) * Math.cos(t1) * Math.cos(n * t1) - Math.sin(t1) * Math.sin(n * t1) - n2 * (-1 + 2 * n2) * Math.sin(t1) * Math.sin(n * t1)))) / (1 + 4 * n4);


                d3Ayc0[n] = 2 * bmu * (bJ * Inm1bJ * (-(1 - n2 + 4 * n4) * Math.cos(t1) * Math.cos(n * t1) + n * (1 + n2 + 2 * n4) * Math.sin(t1) * Math.sin(n * t1))
                        + InbJ * (-(-n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 - n2 + 4 * n4)) * Math.cos(t1) * Math.cos(n * t1) + n * (-(1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 + n2 + 2 * n4)) * Math.sin(t1) * Math.sin(n * t1))
                ) / (bJ + 4 * bJ * n4);


                d3Ays0[n] = -bmu * (-(-1 + n) * (-1 + n) * (-1 + n) * (1 + 2 * n + 2 * n2) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (n * (1 + n2 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + (1 - n2 + 4 * n4) * Math.cos(t1) * Math.sin(n * t1))
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                ) / (1 + 4 * n4);

                d2Ayc1[n] = 0.25 * bmu * bmu * ((-2 + n) * (-2 + n) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * (-2 + n) * (-2 + n) * n * (2 + 2 * n + n2) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Ays1[n] = 0.25 * bmu * bmu * ((n - 2) * (n - 2) * Inm2bJ * Math.sin((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );
            }
        }


        double p0 = Math.exp(bJ * Math.cos(t1 - t2));
        double pM1 = 1 / p0;
        double pM2 = pM1 / p0;
        double lnp1 = -bJ * Math.sin(t1 - t2);
        double lnp11 = -bJ * Math.cos(t1 - t2);
        double lnp2 = -lnp1;
        double px1 = bmu * p0 * (Math.cos(t1) + Math.cos(t2));
        double py1 = bmu * p0 * (Math.sin(t1) + Math.sin(t2));
        psix1 = 0;
        psix2 = 0;
        psix11 = 0;
        psix12 = 0;
        psix22 = 0;
        psi1x1 = 0;
        psi1x2 = 0;
        psiy1 = 0;
        psiy2 = 0;
        psiy11 = 0;
        psiy12 = 0;
        psiy22 = 0;
        psi1y1 = 0;
        psi1y2 = 0;
        psi1x11 = 0;
        psi1y11 = 0;
        psi1x22 = 0;
        psi1y22 = 0;
        psix111 = 0;
        psiy111 = 0;
        psix112 = 0;
        psiy112 = 0;
        psix221 = 0;
        psiy221 = 0;
        psix222 = 0;
        psiy222 = 0;
        psi1x11 = 0;
        psi1y11 = 0;
        psi1x22 = 0;
        psi1y22 = 0;
        for (int n = 0; n <= nMax; n++) {
            psix1 += dAxs0[n] * Math.sin(n * t2) + dAxc0[n] * Math.cos(n * t2);
            psiy1 += dAys0[n] * Math.sin(n * t2) + dAyc0[n] * Math.cos(n * t2);

            psix2 += n * Axs0[n] * Math.cos(n * t2) - n * Axc0[n] * Math.sin(n * t2);
            psiy2 += n * Ays0[n] * Math.cos(n * t2) - n * Ayc0[n] * Math.sin(n * t2);

            psi1x1 += dAxs1[n] * Math.sin(n * t2) + dAxc1[n] * Math.cos(n * t2);
            psi1y1 += dAys1[n] * Math.sin(n * t2) + dAyc1[n] * Math.cos(n * t2);


            psi1x2 += n * Axs1[n] * Math.cos(n * t2) - n * Axc1[n] * Math.sin(n * t2);
            psi1y2 += n * Ays1[n] * Math.cos(n * t2) - n * Ayc1[n] * Math.sin(n * t2);

            psix11 += d2Axs0[n] * Math.sin(n * t2) + d2Axc0[n] * Math.cos(n * t2);
            psiy11 += d2Ays0[n] * Math.sin(n * t2) + d2Ayc0[n] * Math.cos(n * t2);

            psix112 += (n * d2Axs0[n] * Math.cos(n * t2)) - (n * d2Axc0[n] * Math.sin(n * t2));
            psiy112 += (n * d2Ays0[n] * Math.cos(n * t2)) - (n * d2Ayc0[n] * Math.sin(n * t2));

            psix111 += d3Axs0[n] * Math.sin(n * t2) + d3Axc0[n] * Math.cos(n * t2);
            psiy111 += d3Ays0[n] * Math.sin(n * t2) + d3Ayc0[n] * Math.cos(n * t2);

            psi1x11 += d2Axs1[n] * Math.sin(n * t2) + d2Axc1[n] * Math.cos(n * t2);
            psi1y11 += d2Ays1[n] * Math.sin(n * t2) + d2Ayc1[n] * Math.cos(n * t2);

            psix12 += n * dAxs0[n] * Math.cos(n * t2) - n * dAxc0[n] * Math.sin(n * t2);
            psiy12 += n * dAys0[n] * Math.cos(n * t2) - n * dAyc0[n] * Math.sin(n * t2);

            psix22 += -n * n * Axs0[n] * Math.sin(n * t2) - n * n * Axc0[n] * Math.cos(n * t2);
            psiy22 += -n * n * Ays0[n] * Math.sin(n * t2) - n * n * Ayc0[n] * Math.cos(n * t2);

            psix221 += -n * n * dAxs0[n] * Math.sin(n * t2) - n * n * dAxc0[n] * Math.cos(n * t2);
            psiy221 += -n * n * dAys0[n] * Math.sin(n * t2) - n * n * dAyc0[n] * Math.cos(n * t2);

            psix222 += -n * n * n * Axs0[n] * Math.cos(n * t2) + n * n * n * Axc0[n] * Math.sin(n * t2);
            psiy222 += -n * n * n * Ays0[n] * Math.cos(n * t2) + n * n * n * Ayc0[n] * Math.sin(n * t2);

            psi1x22 += -n * n * Axs1[n] * Math.sin(n * t2) - n * n * Axc1[n] * Math.cos(n * t2);
            psi1y22 += -n * n * Ays1[n] * Math.sin(n * t2) - n * n * Ayc1[n] * Math.cos(n * t2);
        }

        double vEx1 = psix1 / p0;
        double vEx2 = psix2 / p0;
        double vEy1 = psiy1 / p0;
        double vEy2 = psiy2 / p0;

        double vEEx1 = psi1x1 / p0 - px1 * vEx1 / p0;
        double vEEx2 = psi1x2 / p0 - px1 * vEx2 / p0;

        double vEEy1 = psi1y1 / p0 - py1 * vEy1 / p0;
        double vEEy2 = psi1y2 / p0 - py1 * vEy2 / p0;

//        double vDotGradvx1 = pM2 * (-lnp1 * psix1 * (psix1 - psix2) + psix1 * psix11 + psix2 * psix12);
//        double vDotGradvy1 = pM2 * (-lnp1 * psiy1 * (psiy1 - psiy2) + psiy1 * psiy11 + psiy2 * psiy12);
//        double vDotGradvx2 = pM2 * (-lnp1 * psix2 * (psix1 - psix2) + psix1 * psix12 + psix2 * psix22);
//        double vDotGradvy2 = pM2 * (-lnp1 * psiy2 * (psiy1 - psiy2) + psiy1 * psiy12 + psiy2 * psiy22);

        vE.setX(0,vEx1);
        vE.setX(1,vEy1);
        vEE.setX(0,vEEx1);
        vEE.setX(1,vEEy1);
//        vDotGradV.setX(0,vDotGradvx1);
//        vDotGradV.setX(1,vDotGradvy1);
        leafAgentManager.getAgent(atom1).vE().PE(vE);
        leafAgentManager.getAgent(atom1).vEE().PE(vEE);
//        leafAgentManager.getAgent(atom1).vDotGradV().PE(vDotGradV);

        vE.setX(0,vEx2);
        vE.setX(1,vEy2);
        vEE.setX(0,vEEx2);
        vEE.setX(1,vEEy2);
//        vDotGradV.setX(0,vDotGradvx2);
//        vDotGradV.setX(1,vDotGradvy2);
        leafAgentManager.getAgent(atom2).vE().PE(vE);//TODO it should be TODO right
        leafAgentManager.getAgent(atom2).vEE().PE(vEE);
//        leafAgentManager.getAgent(atom2).vDotGradV().PE(vDotGradV);





    }

    public void setAgentManager(AtomLeafAgentManager agentManager) {
        leafAgentManager = agentManager;
    }

    public void reset() {
        leafAgentManager.getAgents().values().forEach((agent) -> agent.phi().E(0));
    }


}
