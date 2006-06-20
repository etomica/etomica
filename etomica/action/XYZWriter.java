package etomica.action;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;

/**
 * Action that dumps a phase's configuration to an XYZ file.  Arbitrary but 
 * unique elements are assigned to each atom type.  After writing the PDB file,
 * writeRasmolScript can be called to write a script that will properly 
 * initialize the atomic radii.
 */
public class XYZWriter implements Action, Serializable {

    public XYZWriter(Phase aPhase) {
        iterator = new AtomIteratorLeafAtoms(aPhase);
    }

    public String getLabel() {
        return "XYZ Writer";
    }
    
    /**
     * Sets the file to write to.  This method (or setFileName) must be called
     * before calling actionPerformed and again before calling 
     * writeRasmolScript.
     */
    public void setFile(File newFile) {
        file = newFile;
    }
    
    /**
     * Sets the file name to write to.  This method (or setFile) must be called
     * before calling actionPerformed and again before calling 
     * writeRasmolScript.
     */
    public void setFileName(String fileName) {
        file = new File(fileName);
    }

    public void actionPerformed() {
        if (file == null) {
            throw new IllegalStateException("must call setFile or setFileName before actionPerformed");
        }
        FileWriter fileWriter;
        try { 
            fileWriter = new FileWriter(file);
        }catch(IOException e) {
            System.err.println("Cannot open "+file.getPath()+", caught IOException: " + e.getMessage());
            return;
        }
        try {
            iterator.reset();
            elementAtomType.clear();
            fileWriter.write(Integer.toString(iterator.size())+"\n");
            fileWriter.write("#\n");
            while (iterator.hasNext()) {
                AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
                Iterator elementIterator = elementAtomType.iterator();
                int elementIndex = -1;
                while (elementIterator.hasNext()) {
                    ElementLinker thisElement = (ElementLinker)elementIterator.next();
                    if (thisElement.type == atom.type) {
                        elementIndex = thisElement.elementIndex;
                    }
                }
                if (elementIndex == -1) {
                    ElementLinker thisElement = new ElementLinker(elementCount,atom.type);
                    elementIndex = thisElement.elementIndex;
                    elementCount++;
                    elementAtomType.add(thisElement);
                }
                fileWriter.write(elements[elementIndex]+" "+atom.coord.position().x(0)+" "+atom.coord.position().x(1)+" "+atom.coord.position().x(2)+"\n");
            }
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+file.getPath()+", caught IOException: " + e.getMessage());
        }
    }

    /**
     * Writes a script for rasmol that initializes the radii of each atom type.
     */
    public void writeRasmolScript() {
        if (file == null) {
            throw new IllegalStateException("must call setFile or setFileName before actionPerformed");
        }
        if (file.getAbsolutePath().matches("\\.xyz")) {
            throw new IllegalStateException("must call setFile or setFileName before writeRasmolScript");
        }
        FileWriter fileWriter;
        try { 
            fileWriter = new FileWriter(file);
        }catch(IOException e) {
            System.err.println("Cannot open "+file.getPath()+", caught IOException: " + e.getMessage());
            return;
        }
        try {
            Iterator elementIterator = elementAtomType.iterator();
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = (ElementLinker)elementIterator.next();
                fileWriter.write("select elemno="+elementNum[thisElement.elementIndex]+"\n");
                fileWriter.write("spacefill "+((AtomTypeSphere)thisElement.type).diameter(null)*0.5);
            }
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+file.getPath()+", caught IOException: " + e.getMessage());
        }
    }
    
    private File file;
    private static final char[] elements = new char[] {'H', 'O', 'F', 'N', 'C', 'P', 'S'};
    private static final int[] elementNum = new int[] {1, 8, 9, 7, 6, 15, 16};
    private int elementCount = 0;
    private final LinkedList elementAtomType = new LinkedList();
    private final AtomIteratorLeafAtoms iterator;
    
    private static final class ElementLinker implements Serializable {
        public final int elementIndex;
        public final AtomType type;
        public ElementLinker(int aElementIndex, AtomType aType) {
            elementIndex = aElementIndex;
            type = aType;
        }
    }
}
