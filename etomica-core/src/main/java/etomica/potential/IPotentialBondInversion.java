package etomica.potential;
/**
 * Interface for a potential capable of return bond torsion energies and derivatives.
 * Because it only takes cos(theta), the u(-theta) = u(theta) and du/dtheta = 0 at
 * theta=0 and theta=pi.
 */
public interface IPotentialBondInversion {

        /**
         * Returns the energy for the given value of cos(theta).
         */
        double u(double costheta);

        /**
         * Returns the energy (u) and du/d(cos(theta)) (du) as out parameters for the
         * given value of cos(theta).
         */
        void udu(double costheta, double[] u, double[] du);

}
