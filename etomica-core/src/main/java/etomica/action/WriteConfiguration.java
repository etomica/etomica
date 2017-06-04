/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * Dumps a box's configuration to a file.  The coordinates are written in a 
 * format that can be read in by ConfigurationFile.  The output file has a 
 * "pos_new" extension, which should be renamed to "pos" for use with
 * ConfigurationFile.
 */
public class WriteConfiguration implements IAction {

	public WriteConfiguration(Space space) {
        writePosition = space.makeVector();
        setDoApplyPBC(true);
	}

    /**
     * Sets the configuration name.  The file written to is newConfName.pos_new
     */
    public void setConfName(String newConfName) {
        confName = newConfName;
        fileName = newConfName+".pos_new";
    }

    /**
     * Returns the configuration name.  The file written to is confName.pos_new
     */
    public String getConfName() {
        return confName;
    }

    public void setFileName(String newFileName) {
        fileName = newFileName;
    }

    /**
     * Sets the box whose atom coordinates get written to the file.
     */
    public void setBox(Box newBox) {
        box = newBox;
    }

    /**
     * Returns the box whose atom coordinates get written to the file.
     */
    public Box getBox() {
        return box;
    }

    /**
     * Directs the writer to apply periodic boundary conditions or not (true 
     * by default).
     */
    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }

    /**
     * Returns true if PBC are applied to coordinates written to the file.
     */
    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }

    /**
     * Writes the leaf Atom coordinates to the file confName.pos_new.  If the
     * file exists, it is overwritten.
     */
    public void actionPerformed() {
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(fileName);
        }catch(IOException e) {
            System.err.println("Cannot open "+fileName+", caught IOException: " + e.getMessage());
            return;
        }
        try {
            IAtomList leafList = box.getLeafList();
            int nLeaf = leafList.getAtomCount();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.getAtom(iLeaf);
                writeAtom(fileWriter, a);
            }
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
        }
    }

    protected void writeAtom(FileWriter fileWriter, IAtom a) throws IOException {
        writePosition.E(a.getPosition());
        if (doApplyPBC) {
            Vector shift = box.getBoundary().centralImage(writePosition);
            if (!shift.isZero()) {
                writePosition.PE(shift);
            }
        }
        
        fileWriter.write(writePosition.getX(0)+"");
        for (int i=1; i<writePosition.getD(); i++) {
            fileWriter.write(" "+writePosition.getX(i));
        }
        fileWriter.write("\n");
    }

    protected String confName, fileName;
    protected Box box;
    protected boolean doApplyPBC;
    protected final Vector writePosition;

}
