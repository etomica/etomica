/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.box.Box;
import etomica.lattice.crystal.Primitive;
import etomica.molecule.MoleculePositionCOMPBC;
import etomica.potential.PotentialCallbackMoleculeHessian;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

public class LatticeDynamics implements PotentialCallbackMoleculeHessian.HessianConsumer {

    protected final int molD;
    protected final WaveVectorFactory kFactory;
    protected final Vector[] waveVectors;

    protected final Box box;
    protected final Tensor[][][][] matrix;
    protected final int numBasis;
    protected final Vector[] cellPos;

    public LatticeDynamics(SpeciesManager sm, Box box, Primitive primitive, int numBasis) {
        this.box = box;
        this.numBasis = numBasis;
        kFactory = new WaveVectorFactorySimple(primitive, box.getSpace());
        kFactory.makeWaveVectors(box);
        double[] kCoefficients = kFactory.getCoefficients(); //kCoefficients=0.5 non-deg.; = 1 degenerate twice!

        molD = sm.getSpecies(0).getLeafAtomCount() > 1 ? 2 : 1;
        waveVectors = kFactory.getWaveVectors();
        int numWaveVectors = waveVectors.length;
        matrix = new Tensor[numWaveVectors][numBasis*molD][numBasis*molD][2];
        for (int i=0; i<numWaveVectors; i++) {
            for (int j=0; j<numBasis*molD; j++) {
                for (int k=0; k<numBasis*molD; k++) {
                    for (int l=0; l<2; l++) {
                        matrix[i][j][k][l] = box.getSpace().makeTensor();
                    }
                }
            }
        }
        int numCells = box.getMoleculeList().size() / numBasis;
        cellPos = new Vector[numCells];
        for (int i=0; i<numCells; i++) {
            cellPos[i] = box.getSpace().makeVector();
            // complete failure for mixtures
            cellPos[i].E(MoleculePositionCOMPBC.com(box.getBoundary(), box.getMoleculeList().get(i*numBasis)));
        }
    }

    @Override
    public void takeHessian(int i, int j, Tensor tt, Tensor tr, Tensor rt, Tensor rr) {
        // we take i to be the molecule within the first unit cell and j to be any molecule in the box
        if (i >= numBasis) return;
        for (int iwv=0; iwv<waveVectors.length; iwv++) {
            double kdotr = waveVectors[iwv].dot(cellPos[j/numBasis]);
            double c = Math.cos(kdotr);
            double s = -Math.sin(kdotr);
            matrix[iwv][molD*i][molD*(j%numBasis)][0].PEa1Tt1(c, tt);
            matrix[iwv][molD*i][molD*(j%numBasis)][1].PEa1Tt1(s, tt);

            if (molD > 1) {
                matrix[iwv][molD * i][molD * (j % numBasis) + 1][0].PEa1Tt1(c, tr);
                matrix[iwv][molD * i][molD * (j % numBasis) + 1][1].PEa1Tt1(s, tr);

                matrix[iwv][molD * i + 1][molD * (j % numBasis)][0].PEa1Tt1(c, rt);
                matrix[iwv][molD * i + 1][molD * (j % numBasis)][1].PEa1Tt1(s, rt);

                matrix[iwv][molD * i + 1][molD * (j % numBasis) + 1][0].PEa1Tt1(c, rr);
                matrix[iwv][molD * i + 1][molD * (j % numBasis) + 1][1].PEa1Tt1(s, rr);
            }
        }
    }

    public Tensor[][][][] getMatrix() {
        return matrix;
    }
}
