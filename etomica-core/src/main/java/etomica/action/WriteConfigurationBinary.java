/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * Dumps a box's configuration to a file.  The coordinates are serialized to a
 * file as a 2D array of doubles.  The output file can be read by
 * ConfigurationFileBinary and has a "pos_new" extension by default, which
 * should be renamed to "pos" for use with ConfigurationFileBinary.
 */
public class WriteConfigurationBinary implements IAction {

	public WriteConfigurationBinary(Space space) {
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
        IAtomList leafList = box.getLeafList();
        int dim = 0;
        int nLeaf = leafList.getAtomCount();
        if (nLeaf > 0) dim = leafList.getAtom(0).getPosition().getD();
        double[][] x = new double[nLeaf][dim];
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            writePosition.E(a.getPosition());
            if (doApplyPBC) {
                Vector shift = box.getBoundary().centralImage(writePosition);
                if (!shift.isZero()) {
                    writePosition.PE(shift);
                }
            }
            
            writePosition.assignTo(x[iLeaf]);
        }
        try {
            FileOutputStream fos = null;
            ObjectOutputStream out = null;
            fos = new FileOutputStream(fileName);
            out = new ObjectOutputStream(fos); 
            out.writeObject(x);
            out.close();
            fos.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private String confName, fileName;
    private Box box;
    private boolean doApplyPBC;
    protected final Vector writePosition;

}
