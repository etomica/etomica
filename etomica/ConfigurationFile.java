package etomica;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.iterator.AtomIteratorTree;
import etomica.space.Vector;

/**
 * reads configuration coordinates from a file and assigns them to the leaf atoms in a phase
 */
public class ConfigurationFile extends Configuration {

    public ConfigurationFile(Space space, String aConfName) {
        super(space);
        confName = aConfName;
        newPos = space.makeVector();
    }
    
    public void initializePositions(AtomIterator[] iterators){
        if (iterators.length == 0) return;
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
            AtomIteratorTree leafIterator = new AtomIteratorTree();
            int iIterator = 0;
            iterator.reset();
            while (!iterator.hasNext() && iIterator < iterators.length-1) {
                iterator = iterators[++iIterator];
                iterator.reset();
            }
            while (iterator.hasNext()) {
                Atom molecule = iterator.nextAtom();
                if (molecule.node.isLeaf()) {
                    setPosition(molecule,bufReader.readLine());
                }
                else {
                    leafIterator.setRoot(molecule);
                    leafIterator.reset();
                    while (leafIterator.hasNext()) {
                        setPosition(leafIterator.nextAtom(),bufReader.readLine());
                    }
                }
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
    
    private void setPosition(Atom atom, String string) {
        String[] coordStr = string.split(" +");
        double[] coord = new double[coordStr.length];
        for (int i=0; i<coord.length; i++) {
            coord[i] = Double.valueOf(coordStr[i]).doubleValue();
        }
        newPos.E(coord);
        atom.coord.position().E(newPos);
    }
    
    private String confName;
    private Vector newPos;
}
