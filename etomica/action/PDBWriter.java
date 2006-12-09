package etomica.action;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;

/**
 * Action that dumps a phase's configuration to an PDB file.  Arbitrary but 
 * unique elements are assigned to each atom type.  After writing the PDB file,
 * writeRasmolScript can be called to write a script that will properly 
 * initialize the atomic radii.
 */
public class PDBWriter implements Action, Serializable {

    public PDBWriter(Phase aPhase) {
        iterator = new AtomIteratorLeafAtoms(aPhase);
        try {
            Class.forName("java.util.Formatter");
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException("JRE 1.5 is required for PDBWriter.  Try XYZWriter instead");
        }
    }

    public String getLabel() {
        return "PDB Writer";
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
        Class formatterClass;
        try {
            formatterClass = Class.forName("java.util.Formatter");
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException("JRE 1.5 is required for PDBWriter.  Try XYZWriter instead");
        }
        Object formatter;
        Method formatMethod;
        try {
            formatMethod = formatterClass.getDeclaredMethod("format",new Class[]{String.class,new Object[0].getClass()});
            formatter = formatterClass.getDeclaredConstructor(new Class[]{File.class}).newInstance(new Object[]{file});
        }
        catch (NoSuchMethodException e) {
            throw new RuntimeException(e);
        }
        catch (InvocationTargetException e) {
            throw new RuntimeException(e);
        }
        catch (InstantiationException e) {
            throw new RuntimeException(e);
        }
        catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
        iterator.reset();
        elementAtomType.clear();
        int atomCount = 0;
        while (iterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
            Iterator elementIterator = elementAtomType.iterator();
            int elementIndex = -1;
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = (ElementLinker)elementIterator.next();
                if (thisElement.type == atom.getType()) {
                    elementIndex = thisElement.elementIndex;
                }
            }
            if (elementIndex == -1) {
                ElementLinker thisElement = new ElementLinker(elementCount,atom.getType());
                elementIndex = thisElement.elementIndex;
                elementCount++;
                elementAtomType.add(thisElement);
            }
            try {
                formatMethod.invoke(formatter,new Object[]{"ATOM%7d%3s                %8.3f%8.3f%8.3f\n", new Object[]{new Integer(atomCount), new Character(elements[elementIndex]), 
                        new Double(atom.getCoord().position().x(0)), new Double(atom.getCoord().position().x(1)), new Double(atom.getCoord().position().x(2))}});
            }
            catch (InvocationTargetException e) {
                throw new RuntimeException(e);
            }
            catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }
            atomCount++;
        }
        try {
            formatterClass.getDeclaredMethod("close",new Class[]{}).invoke(formatter,null);
        }
        catch (NoSuchMethodException e) {
            throw new RuntimeException(e);
        }
        catch (InvocationTargetException e) {
            throw new RuntimeException(e);
        }
        catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
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
            Iterator elementIterator = elementAtomType.iterator();
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = (ElementLinker)elementIterator.next();
                fileWriter.write("select elemno="+elementNum[thisElement.elementIndex]+"\n");
                fileWriter.write("spacefill "+((AtomTypeSphere)thisElement.type).getDiameter()*0.5);
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
    private final LinkedList elementAtomType = new LinkedList();
    private final AtomIteratorLeafAtoms iterator;
    
    private static final class ElementLinker implements Serializable {
        private static final long serialVersionUID = 1L;
        public final int elementIndex;
        public final AtomType type;
        public ElementLinker(int aElementIndex, AtomType aType) {
            elementIndex = aElementIndex;
            type = aType;
        }
    }
}
