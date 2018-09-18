/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.ElementChemical;
import etomica.chem.elements.IElement;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Action that dumps a box's configuration to an XYZ file.  Arbitrary but 
 * unique elements are assigned to each atom type.  After writing the PDB file,
 * writeRasmolScript can be called to write a script that will properly 
 * initialize the atomic radii.
 */
public class XYZWriter implements IAction, Serializable {

    private static final long serialVersionUID = 1L;
    private static final String[] elements = new String[]{"H", "O", "F", "N", "C", "P", "S"};
    private static final int[] elementNum = new int[]{1, 8, 9, 7, 6, 15, 16};
    private final LinkedList<ElementLinker> elementAtomType;
    private final IAtomList leafList;
    private File file;
    private int elementCount = 0;
    private boolean doAppend;
    public XYZWriter(Box aBox) {
        leafList = aBox.getLeafList();
        elementAtomType = new LinkedList<ElementLinker>();
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

    public void setIsAppend(boolean newDoAppend) {
        doAppend = newDoAppend;
    }

    public boolean isAppend() {
        return doAppend;
    }

    public void actionPerformed() {
        if (file == null) {
            throw new IllegalStateException("must call setFile or setFileName before actionPerformed");
        }
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(file,doAppend);
        }catch(IOException e) {
            System.err.println("Cannot open "+file.getPath()+", caught IOException: " + e.getMessage());
            return;
        }
        try {
            fileWriter.write(Integer.toString(leafList.size())+"\n");
            fileWriter.write("#\n");
            int nLeaf = leafList.size();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom atom = leafList.get(iLeaf);
                IElement element = atom.getType().getElement();
                String symbol = element.getSymbol();
                if (!(element instanceof ElementChemical)) {
                    Iterator<ElementLinker> elementIterator = elementAtomType.iterator();
                    int elementIndex = -1;
                    while (elementIterator.hasNext()) {
                        ElementLinker thisElement = elementIterator.next();
                        if (thisElement.type == atom.getType()) {
                            elementIndex = thisElement.elementIndex;
                            break;
                        }
                    }
                    if (elementIndex == -1) {
                        ElementLinker thisElement = new ElementLinker(elementCount, atom.getType());
                        elementIndex = thisElement.elementIndex;
                        elementCount++;
                        elementAtomType.add(thisElement);
                    }
                    symbol = elements[elementIndex];
                }
                fileWriter.write(symbol+" "+atom.getPosition().getX(0)+" "+atom.getPosition().getX(1)+" "+atom.getPosition().getX(2)+"\n");
            }
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+file.getPath()+", caught IOException: " + e.getMessage());
        }
    }

    /**
     * Writes a script for rasmol that initializes the radii of each atom type.
     * The DiameterHashByType is used to provide radii.
     */
    public void writeRasmolScript(DiameterHashByType diameterHash) {
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
            Iterator<ElementLinker> elementIterator = elementAtomType.iterator();
            while (elementIterator.hasNext()) {
                ElementLinker thisElement = elementIterator.next();
                fileWriter.write("select elemno="+elementNum[thisElement.elementIndex]+"\n");
                fileWriter.write("spacefill "+diameterHash.getDiameter(thisElement.type)*0.5);
            }
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+file.getPath()+", caught IOException: " + e.getMessage());
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
