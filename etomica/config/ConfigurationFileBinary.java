package etomica.config;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

import etomica.action.BoxInflate;
import etomica.action.WriteConfigurationBinary;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresMono;

/**
 * reads configuration coordinates from a binary file and assigns them to the
 * leaf atoms in a box.  The file format can be written by
 * WriteConfigurationBinary.
 */
public class ConfigurationFileBinary implements Configuration {

    public ConfigurationFileBinary(String aConfName) {
        confName = aConfName;
    }
    
    public void initializeCoordinates(IBox box) {
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
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            IVectorMutable p = a.getPosition();
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
    public void replicate(int numAtoms0, double density, int[] reps, String outConfigname, ISpace space) {
        Simulation sim = new Simulation(space);
        IBox box0 = new Box(space);
        sim.addBox(box0);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        box0.setNMolecules(species, numAtoms0);
        BoxInflate inflater = new BoxInflate(box0, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        initializeCoordinates(box0);

        int numAtoms1 = numAtoms0;
        for (int i=0; i<reps.length; i++) {
            numAtoms1 *= reps[i];
        }
        IBox box1 = new Box(space);
        sim.addBox(box1);
        box1.setNMolecules(species, numAtoms1);
        inflater.setBox(box1);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        IAtomList leafList0 = box0.getLeafList();
        IAtomList leafList1 = box1.getLeafList();
        double[] xyzShift = new double[3];
        for (int i=0; i<reps[0]; i++) {
            xyzShift[0] = box0.getBoundary().getBoxSize().getX(0)*(-0.5*(reps[0]-1) + i); 
            for (int j=0; j<reps[1]; j++) {
                xyzShift[1] = box0.getBoundary().getBoxSize().getX(1)*(-0.5*(reps[1]-1) + j); 
                for (int k=0; k<reps[2]; k++) {
                    xyzShift[2] = box0.getBoundary().getBoxSize().getX(2)*(-0.5*(reps[2]-1) + k);
                    int start1 = numAtoms0*(i*reps[2]*reps[1] + j*reps[2] + k);
                    for (int iAtom = 0; iAtom<numAtoms0; iAtom++) {
                        IVectorMutable p0 = leafList0.getAtom(iAtom).getPosition();
                        IVectorMutable p1 = leafList1.getAtom(start1+iAtom).getPosition();
                        for (int l=0; l<3; l++) {
                            p1.setX(l, p0.getX(l) + xyzShift[l]);
                        }
                    }
                }
            }
        }
        WriteConfigurationBinary writeConfig = new WriteConfigurationBinary(space);
        writeConfig.setFileName(outConfigname+".pos");
        writeConfig.setBox(box1);
        writeConfig.actionPerformed();
    }

    /**
     * Reads the configuration from the file and then rescales the
     * configuration to the new density (density1).  The new configuration is
     * then written out to outFilename (also as binary).  The number of atoms
     * and original density are required as input.
     */
    public void rescale(int numAtoms, double density0, double density1, String outConfigname, ISpace space) {
        Simulation sim = new Simulation(space);
        IBox box = new Box(space);
        sim.addBox(box);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density0);
        inflater.actionPerformed();
        initializeCoordinates(box);

        inflater.setTargetDensity(density1);
        inflater.actionPerformed();
        WriteConfigurationBinary writeConfig = new WriteConfigurationBinary(space);
        writeConfig.setFileName(outConfigname+".pos");
        writeConfig.setBox(box);
        writeConfig.actionPerformed();
    }

    protected final String confName;
}
