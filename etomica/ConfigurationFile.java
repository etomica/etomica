package etomica;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * reads configuration coordinates from a file and assigns them to the leaf atoms in a phase
 */
public class ConfigurationFile extends Configuration {

    public ConfigurationFile(Space space, String aConfName) {
        super(space);
        confName = aConfName;
    }
    
    public void initializePositions(AtomIterator[] iterators){
        String fileName = confName+".pos";
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            AtomIterator iterator = iterators[0];
            int iIterator = 0;
            iterator.reset();
            while (!iterator.hasNext() && iIterator < iterators.length-1) {
                iterator = iterators[++iIterator];
                iterator.reset();
            }
            while (iterator.hasNext()) {
                String string = bufReader.readLine();
                Atom atom = iterator.nextAtom();
                Space.Vector newPos = (Space.Vector)atom.coord.position().clone();
                String[] coordStr = string.split(" +");
                double[] coord = new double[coordStr.length];
                for (int i=0; i<coord.length; i++) {
                    coord[i] = Double.valueOf(coordStr[i]).doubleValue();
                }
                newPos.E(coord);
                atom.coord.displaceTo(newPos);
                while (!iterator.hasNext() && iIterator < iterators.length-1) {
                    iterator = iterators[++iIterator];
                    iterator.reset();
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
        }
    }
        
    private String confName;
}
