/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.space.Space;

/**
 * This class represents a Mayer function for any pair of species from a
 * hard sphere mixture.  The function is 0 when the spheres do not overlap,
 * and a constant value when they do.  The class requires the # of types
 * (species) and the diameter of each.  It takes an optional parameter to set
 * non-additivity.  There are convenience methods to create an instance of this
 * class for a reference or target system.
 * <p>
 * For the reference system, the function is positive and is divided by the
 * volume of the integration for the pair, such that the integral of the
 * function will be 1.
 * <p>
 * For the target system, the function is -1.
 *
 * @author Pavan and Andrew
 */
public class MayerHSMixture implements MayerFunction {

    protected final double[][] sigma2;
    protected final double[][] allVal;

    public static MayerHSMixture makeReferenceF(Space space, int nPoints, double[][] pairSigma) {
        return new MayerHSMixture(nPoints, pairSigma, 1, true, space);
    }

    public static MayerHSMixture makeTargetF(int nPoints, double[][] pairSigma) {
        return new MayerHSMixture(nPoints, pairSigma, -1, false, null);
    }

    protected MayerHSMixture(int nPoints, double[][] pairSigma, double val, boolean isRef, Space space) {
        sigma2 = new double[nPoints][nPoints];
        allVal = new double[nPoints][nPoints];

        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j < nPoints; j++) {
                sigma2[i][j] = sigma2[j][i] = pairSigma[i][j] * pairSigma[i][j];
                double vij = isRef ? space.sphereVolume(pairSigma[i][j]) : 1;
                allVal[i][j] = allVal[j][i] = val / vij;
            }
        }
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        int i = pair.getMolecule(0).getIndex();
        int j = pair.getMolecule(1).getIndex();
        return (r2 < sigma2[i][j]) ? allVal[i][j] : 0.0;
    }

    public IPotential getPotential() {
        return null;
    }

    public void setBox(Box newBox) {
    }
}
