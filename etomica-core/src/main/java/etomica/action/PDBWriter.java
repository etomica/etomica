/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;

import java.io.*;
import java.util.Formatter;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Action that dumps a box's configuration to an PDB file.  Arbitrary but 
 * unique elements are assigned to each atom type.  After writing the PDB file,
 * writeRasmolScript can be called to write a script that will properly 
 * initialize the atomic radii.
 */
public class PDBWriter implements IAction {

    private static final char[] elements = new char[]{'H', 'O', 'F', 'N', 'C', 'P', 'S'};
    private static final int[] elementNum = new int[]{1, 8, 9, 7, 6, 15, 16};
    private final LinkedList<ElementLinker> elementAtomType = new LinkedList<>();
    private File file;
    private int elementCount = 0;
    private IAtomList leafList;
    
    public PDBWriter() {
    }
    public PDBWriter(Box aBox) {
        this();
        setBox(aBox);
    }

    public void setBox(Box newBox) {
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
        for (IAtom atom : leafList) {
            Iterator<ElementLinker> elementIterator = elementAtomType.iterator();
            int elementIndex = -1;
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = elementIterator.next();
                if (thisElement.type == atom.getType()) {
                    elementIndex = thisElement.elementIndex;
                }
            }
            if (elementIndex == -1) {
                ElementLinker thisElement = new ElementLinker(elementCount, atom.getType());
                elementIndex = thisElement.elementIndex;
                elementCount++;
                elementAtomType.add(thisElement);
            }
            formatter.format("ATOM%7d%3s                %8.3f%8.3f%8.3f\n", atomCount, elements[elementIndex],
                    atom.getPosition().getX(0), atom.getPosition().getX(1), atom.getPosition().getX(2));
            atomCount++;
        }
        formatter.close();
    }

    /**
     * Writes a script for rasmol that initializes the radii of each atom type.
     * The DiameterHashByType is used to provide radii.
     */
    public void writeRasmolScript(DiameterHashByType diameterHash) {
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
            for (ElementLinker thisElement : elementAtomType) {
                fileWriter.write("select elemno=" + elementNum[thisElement.elementIndex] + "\n");
                fileWriter.write("spacefill " + diameterHash.getDiameter(thisElement.type) * 0.5);
            }
            fileWriter.close();
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }
    
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
