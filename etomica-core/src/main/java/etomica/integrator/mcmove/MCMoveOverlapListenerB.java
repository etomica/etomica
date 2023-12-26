/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.box.Box;

import java.util.Arrays;

public class MCMoveOverlapListenerB extends MCMoveOverlapListener {

    public MCMoveOverlapListenerB(MCMoveInsertDeleteBiased mcMove, int numAlpha, double daDef, double daSpan) {
        super(mcMove, numAlpha, daDef, 0, daSpan);
    }

    public void actionPerformed(MCMoveEvent event) {
        if (event.getMCMove() != mcMove) return;
        Box box = mcMove.getBox();
        int numAtoms = box.getNMolecules(mcMove.getSpecies());
        if (event instanceof MCMoveTrialFailedEvent) {
            // trial failed, but we were still here.  we need to increment our sums here
            // for the histogram.
            if (sumInsert.length < numAtoms+1) {
                sumInsert = Arrays.copyOf(sumInsert, numAtoms+1);
                numInsert = Arrays.copyOf(numInsert, numAtoms+1);
                sumDelete = Arrays.copyOf(sumDelete, numAtoms+1);
                numDelete = Arrays.copyOf(numDelete, numAtoms+1);
            }
//            if (mcMove.lastMoveInsert()) {
//                numInsert[numAtoms]++;
//            }
//            else {
//                numDelete[numAtoms]++;
//            }
        }
        else if (event instanceof MCMoveTrialInitiatedEvent) {
            // x = V/N*Math.exp(-beta*deltaU)
            double x = mcMove.getChi(temperature) * Math.exp(-mcMove.getLnBiasDiff());
            if (mcMove.lastMoveInsert()) {
                numAtoms--;
            }
            if (minNumAtoms > numAtoms) minNumAtoms = numAtoms;
            if (sumInsert.length < numAtoms+1) {
                sumInsert = Arrays.copyOf(sumInsert, numAtoms+1);
                numInsert = Arrays.copyOf(numInsert, numAtoms+1);
                sumDelete = Arrays.copyOf(sumDelete, numAtoms+1);
                numDelete = Arrays.copyOf(numDelete, numAtoms+1);
            }
            if (sumInsert[numAtoms] == null) {
                sumInsert[numAtoms] = new double[numAlpha];
                sumDelete[numAtoms] = new double[numAlpha];
            }
            if (mcMove.lastMoveInsert()) {
                double da = daDef; // + Math.log(((double)(numAtomsLattice-numAtoms))/(numAtoms+1));
                for (int i=0; i<numAlpha; i++) {
                    double iLnAlpha = da + daSpan*(i-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0);
                    sumInsert[numAtoms][i] += 1.0/(1 + Math.exp(iLnAlpha)/x);
                }
                numInsert[numAtoms]++;
            }
            else {
                double da = daDef; // + Math.log(((double)(numAtomsLattice-numAtoms+1))/numAtoms);
                for (int i=0; i<numAlpha; i++) {
                    double iLnAlpha = da + daSpan*(i-(numAlpha-1)*0.5)*2.0/(numAlpha-1.0);
                    sumDelete[numAtoms][i] += 1.0/(Math.exp(iLnAlpha) + 1.0/x);
                }
                numDelete[numAtoms]++;
            }
        }
    }
}
