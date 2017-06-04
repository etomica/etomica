/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.action.BoxInflate;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

public class ConfigurationReplicate {

    /**
     * Replicates the configuration from config within boxBig.  There are reps
     * replicates in each direction.
     */
    public static void replicate(Configuration config, Box boxBig, int[] reps, Space space) {
        Simulation sim = new Simulation(space);
        Box box0 = new Box(space);
        sim.addBox(box0);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);

        int numAtoms1 = boxBig.getLeafList().getAtomCount();
        int numAtoms0 = numAtoms1;
        for (int i=0; i<reps.length; i++) {
            numAtoms0 /= reps[i];
        }
        box0.setNMolecules(species, numAtoms0);
        BoxInflate inflater = new BoxInflate(box0, space);
        inflater.setTargetDensity(numAtoms1/boxBig.getBoundary().volume());
        inflater.actionPerformed();
        config.initializeCoordinates(box0);
        IAtomList leafList0 = box0.getLeafList();
        IAtomList leafList1 = boxBig.getLeafList();
        double[] xyzShift = new double[3];
        for (int i=0; i<reps[0]; i++) {
            xyzShift[0] = box0.getBoundary().getBoxSize().getX(0)*(-0.5*(reps[0]-1) + i);
            for (int j=0; j<reps[1]; j++) {
                xyzShift[1] = box0.getBoundary().getBoxSize().getX(1)*(-0.5*(reps[1]-1) + j);
                for (int k=0; k<reps[2]; k++) {
                    xyzShift[2] = box0.getBoundary().getBoxSize().getX(2)*(-0.5*(reps[2]-1) + k);
                    int start1 = numAtoms0*(i*reps[2]*reps[1] + j*reps[2] + k);
                    for (int iAtom = 0; iAtom<numAtoms0; iAtom++) {
                        Vector p0 = leafList0.getAtom(iAtom).getPosition();
                        Vector p1 = leafList1.getAtom(start1+iAtom).getPosition();
                        for (int l=0; l<3; l++) {
                            p1.setX(l, p0.getX(l) + xyzShift[l]);
                        }
                    }
                }
            }
        }
    }

}
