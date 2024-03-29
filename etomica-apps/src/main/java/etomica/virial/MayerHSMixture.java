/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

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
    protected final Boundary boundary;
    protected final Vector dr;

    public static MayerHSMixture makeReferenceF(Space space, int nPoints, double[][] pairSigma, Boundary b) {
        return new MayerHSMixture(nPoints, pairSigma, 1, true, space, b);
    }

    public static MayerHSMixture makeTargetF(Space space, int nPoints, double[][] pairSigma, Boundary b) {
        return new MayerHSMixture(nPoints, pairSigma, -1, false, space, b);
    }

    protected MayerHSMixture(int nPoints, double[][] pairSigma, double val, boolean isRef, Space space, Boundary b) {
        sigma2 = new double[nPoints][nPoints];
        allVal = new double[nPoints][nPoints];
        boundary = b;
        dr = space != null ? space.makeVector() : null;

        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j < nPoints; j++) {
                sigma2[i][j] = sigma2[j][i] = pairSigma[i][j] * pairSigma[i][j];
                //double vij = isRef ? space.sphereVolume(pairSigma[i][j]) : 1;
                double vij = 0;
                if (isRef) {
                    if (b==null || b.getBoxSize().getX(0) >= 2*pairSigma[i][j]) {
                        vij = space.sphereVolume(pairSigma[i][j]);
                    }
                    else if (Math.sqrt(2)*b.getBoxSize().getX(0) >= 2*pairSigma[i][j]){
                        double R = pairSigma[i][j];
                        double h = R - b.getBoxSize().getX(0)/2;
                        double vijExcluded = Math.PI*h*h*(3*R - h)/3;
                        vij = space.sphereVolume(pairSigma[i][j]) -  6*vijExcluded;
                    }
                    else if (Math.sqrt(3)*b.getBoxSize().getX(0) > 2*pairSigma[i][j]) {//Math.sqrt(2)*b.getBoxSize().getX(0) < 2*pairSigma[i][j]) {
                        throw new RuntimeException("Box length not handled. Try increasing or decreasing box length.");
                    }
                    else {//if (Math.sqrt(3)*b.getBoxSize().getX(0) <= 2*pairSigma[i][j]) {
                        vij = b.volume();
                    }
                }
                else {
                    vij = 1;
                }
                allVal[i][j] = allVal[j][i] = val / vij;
            }
        }
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        int i = pair.get(0).getIndex();
        int j = pair.get(1).getIndex();
        if (boundary != null) {
            dr.Ev1Mv2(pair.get(0).getChildList().get(0).getPosition(),
                    pair.get(1).getChildList().get(0).getPosition());
            boundary.nearestImage(dr);
            r2 = dr.squared();
        }
        return (r2 < sigma2[i][j]) ? allVal[i][j] : 0.0;
    }

    public void setBox(Box newBox) {
    }
}
