package etomica.config;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * reads configuration coordinates from a file and assigns them to the leaf atoms in a phase
 */
public class ConfigurationFile extends Configuration {

    public ConfigurationFile(String aConfName) {
        confName = aConfName;
    }
    
    public void initializeCoordinates(Phase phase) {
        AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms(phase);
        String fileName = confName+".pos";
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            atomIterator.reset();
            while (atomIterator.hasNext()) {
                Atom atom = atomIterator.nextAtom();
                setPosition((AtomLeaf)atom,bufReader.readLine());
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
        }
    }
        
    private void setPosition(AtomLeaf atom, String string) {
        String[] coordStr = string.split(" +");
        IVector pos = atom.getPosition();
        for (int i=0; i<pos.getD(); i++) {
            pos.setX(i, Double.valueOf(coordStr[i]).doubleValue());
        }
    }
    
    private static final long serialVersionUID = 2L;
    private String confName;
}
