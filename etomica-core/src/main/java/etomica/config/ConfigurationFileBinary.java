/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.action.BoxInflate;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

/**
 * reads configuration coordinates from a binary file and assigns them to the
 * leaf atoms in a box.  The file format can be written by
 * WriteConfigurationBinary.
 */
public class ConfigurationFileBinary implements Configuration {

    public ConfigurationFileBinary(String aConfName) {
        confName = aConfName;
    }
    
    public void initializeCoordinates(Box box) {
        String fileName = confName+".pos";

        double[][] x;
        FileInputStream fis;
        try {
            fis = new FileInputStream(fileName);
        }
        catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            ObjectInputStream in = new ObjectInputStream(fis);
            x = (double[][])in.readObject();
            in.close();
            fis.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
        finally {
            try {
                fis.close();
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.get(iLeaf);
            Vector p = a.getPosition();
            for (int i=0; i<x[iLeaf].length; i++) {
                p.setX(i,x[iLeaf][i]);
            }
        }
    }

    /**
     * Reads the configuration from the file and then replicates the
     * configuration reps times in each direction.  The new configuration is
     * then written out to outFilename (also as binary).  The density and
     * original number of atoms are required as input.
     */
    public static void replicate(Configuration config, Box box1, int[] reps, Space space) {
        Simulation sim = new Simulation(space);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        Box box0 = new Box(space);
        sim.addBox(box0);

        int numAtoms1 = box1.getLeafList().size();
        int numAtoms0 = numAtoms1;
        for (int i=0; i<reps.length; i++) {
            numAtoms0 /= reps[i];
        }
        box0.setNMolecules(species, numAtoms0);
        
        Vector boundaryBox0 = space.makeVector();
        boundaryBox0.setX(0, box1.getBoundary().getBoxSize().getX(0)/reps[0]);
        boundaryBox0.setX(1, box1.getBoundary().getBoxSize().getX(1)/reps[1]);
        boundaryBox0.setX(2, box1.getBoundary().getBoxSize().getX(2)/reps[2]);
        box0.getBoundary().setBoxSize(boundaryBox0);
        config.initializeCoordinates(box0);
        IAtomList leafList0 = box0.getLeafList();//184 = 4*46
        IAtomList leafList1 = box1.getLeafList();//1472 = 4*46*8
        double[] xyzShift = new double[3];
        for (int i=0; i<reps[0]; i++) {
            xyzShift[0] = box0.getBoundary().getBoxSize().getX(0)*(-0.5*(reps[0]-1) + i); 
            for (int j=0; j<reps[1]; j++) {
                xyzShift[1] = box0.getBoundary().getBoxSize().getX(1)*(-0.5*(reps[1]-1) + j); 
                for (int k=0; k<reps[2]; k++) {
                    xyzShift[2] = box0.getBoundary().getBoxSize().getX(2)*(-0.5*(reps[2]-1) + k);
                    int start1 = numAtoms0*(i*reps[2]*reps[1] + j*reps[2] + k);
                    for (int iAtom = 0; iAtom<numAtoms0; iAtom++) {
                        Vector p0 = leafList0.get(iAtom).getPosition();
                        Vector p1 = leafList1.get(start1+iAtom).getPosition();
                        for (int l=0; l<3; l++) {
                            p1.setX(l, p0.getX(l) + xyzShift[l]);
                        }
                    }
                }
            }
        }
    }

    public static void rescale(Configuration config, Box box1, double density1, Space space) {
        Simulation sim = new Simulation(space);
        Box box0 = new Box(space);
        sim.addBox(box0);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        box0.setNMolecules(species, box1.getLeafList().size());
        BoxInflate inflater = new BoxInflate(box0, space);
        config.initializeCoordinates(box0);
        inflater.setTargetDensity(density1);
        inflater.actionPerformed();
        
        IAtomList atoms0 = box0.getLeafList();
        IAtomList atoms1 = box1.getLeafList();
        for (int i = 0; i<atoms0.size(); i++) {
            atoms1.get(i).getPosition().E(atoms0.get(i).getPosition());
        }
    }

    protected final String confName;
}
