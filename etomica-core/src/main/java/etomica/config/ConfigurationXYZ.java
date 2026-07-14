/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * reads configuration coordinates from an xyz file and assigns them to the leaf atoms in a box
 */
public class ConfigurationXYZ implements Configuration {

    public ConfigurationXYZ(String aFilename) {
        filename = aFilename;
    }

    public void initializeCoordinates(Box box) {
        IAtomList leafList = box.getLeafList();
        FileReader fileReader;
        try {
            fileReader = new FileReader(filename);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+filename+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            int nLeaf = leafList.size();
            String line1 = bufReader.readLine();
            int nLeafCheck = Integer.parseInt(line1);
            if (nLeafCheck != nLeaf) {
                throw new RuntimeException("I have "+nLeaf+" atoms, but the XYZ file has "+nLeafCheck);
            }
            bufReader.readLine();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.get(iLeaf);
                setPosition(a,bufReader.readLine());
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+filename+", caught IOException: " + e.getMessage());
        }
    }

    protected void setPosition(IAtom atom, String string) {
        String[] coordStr = string.split("[ \t]+");
        Vector pos = atom.getPosition();
        for (int i=1; i<=pos.getD(); i++) {
            pos.setX(i-1, Double.valueOf(coordStr[i]).doubleValue());
        }
    }

    protected String filename;
}
