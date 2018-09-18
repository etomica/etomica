/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;

/**
 * reads configuration coordinates from a file and assigns them to the leaf atoms in a box
 */
public class ConfigurationFile implements Configuration, java.io.Serializable {

    public ConfigurationFile(String aConfName) {
        confName = aConfName;
    }
    
    public void initializeCoordinates(Box box) {
        IAtomList leafList = box.getLeafList();
        String fileName = confName+".pos";
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            int nLeaf = leafList.size();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.get(iLeaf);
                setPosition(a,bufReader.readLine());
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
    }
        
    protected void setPosition(IAtom atom, String string) {
        String[] coordStr = string.split("[ \t]+");
        Vector pos = atom.getPosition();
        for (int i=0; i<pos.getD(); i++) {
            pos.setX(i, Double.valueOf(coordStr[i]).doubleValue());
        }
    }
    
    private static final long serialVersionUID = 2L;
    protected String confName;
}
