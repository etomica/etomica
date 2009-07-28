package etomica.action;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.Formatter;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;

/**
 * Action that dumps a box's configuration to an PDB file.  Arbitrary but 
 * unique elements are assigned to each atom type.  After writing the PDB file,
 * writeRasmolScript can be called to write a script that will properly 
 * initialize the atomic radii.
 */
public class PDBWriter implements IAction, Serializable {

    public PDBWriter() {
        try {
            Class.forName("java.util.Formatter");
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException("JRE 1.5 is required for PDBWriter.  Try XYZWriter instead");
        }
    }        
    
    public PDBWriter(IBox aBox) {
        this();
        setBox(aBox);
    }

    public void setBox(IBox newBox) {
        leafList = newBox.getLeafList();
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
        Formatter formatter;
        try {
            formatter = new Formatter(file);
        }
        catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
        int atomCount = 0;
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned atom = (IAtomPositioned)leafList.getAtom(iLeaf);
            Iterator<ElementLinker> elementIterator = elementAtomType.iterator();
            int elementIndex = -1;
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = elementIterator.next();
                if (thisElement.type == ((IAtom)atom).getType()) {
                    elementIndex = thisElement.elementIndex;
                }
            }
            if (elementIndex == -1) {
                ElementLinker thisElement = new ElementLinker(elementCount,((IAtom)atom).getType());
                elementIndex = thisElement.elementIndex;
                elementCount++;
                elementAtomType.add(thisElement);
            }
            formatter.format("ATOM%7d%3s                %8.3f%8.3f%8.3f\n", new Object[]{new Integer(atomCount), new Character(elements[elementIndex]), 
                    new Double(atom.getPosition().getX(0)), new Double(atom.getPosition().getX(1)), new Double(atom.getPosition().getX(2))});
            atomCount++;
        }
        formatter.close();
    }

    /**
     * Writes a script for rasmol that initializes the radii of each atom type.
     */
    public void writeRasmolScript() {
        if (file == null) {
            throw new IllegalStateException("must call setFile or setFileName before actionPerformed");
        }
        if (file.getAbsolutePath().matches("\\.pdb")) {
            throw new IllegalStateException("must call setFile or setFileName before writeRasmolScript");
        }
        FileWriter fileWriter;
        try { 
            fileWriter = new FileWriter(file);
        }catch(IOException e) {
            throw new RuntimeException(e);
        }
        try {
            Iterator<ElementLinker> elementIterator = elementAtomType.iterator();
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = elementIterator.next();
                fileWriter.write("select elemno="+elementNum[thisElement.elementIndex]+"\n");
                fileWriter.write("spacefill "+((IAtomTypeSphere)thisElement.type).getDiameter()*0.5);
            }
            fileWriter.close();
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }
    
    private static final long serialVersionUID = 1L;
    private File file;
    private static final char[] elements = new char[] {'H', 'O', 'F', 'N', 'C', 'P', 'S'};
    private static final int[] elementNum = new int[] {1, 8, 9, 7, 6, 15, 16};
    private int elementCount = 0;
    private final LinkedList<ElementLinker> elementAtomType = new LinkedList<ElementLinker>();
    private IAtomList leafList;
    
    private static final class ElementLinker implements Serializable {
        private static final long serialVersionUID = 1L;
        public final int elementIndex;
        public final IAtomType type;
        public ElementLinker(int aElementIndex, IAtomType aType) {
            elementIndex = aElementIndex;
            type = aType;
        }
    }
}
