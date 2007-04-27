package etomica.config;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
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
        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        String fileName = confName+".pos";
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            int nLeaf = leafList.size();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                AtomLeaf a = (AtomLeaf)leafList.get(iLeaf);
                setPosition(a,bufReader.readLine());
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
