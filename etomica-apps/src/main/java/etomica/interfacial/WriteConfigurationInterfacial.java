/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.WriteConfiguration;
import etomica.api.IAtom;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.space.ISpace;

/**
 * Dumps a box's configuration to a file.  The coordinates are written in a 
 * format that can be read in by ConfigurationFile.  The output file has a 
 * "pos_new" extension, which should be renamed to "pos" for use with
 * ConfigurationFile.
 */
public class WriteConfigurationInterfacial extends WriteConfiguration {

    protected ISpecies species;
    protected IVector shift;
    
	public WriteConfigurationInterfacial(ISpace space) {
	    super(space);
	}
	
	public void setSpecies(ISpecies species) {
	    this.species = species;
	}

	public void setShift(IVector shift) {
	    this.shift = shift;
	}

    protected void writeAtom(FileWriter fileWriter, IAtom a) throws IOException {
        if (a.getType().getSpecies() != species) return;
        writePosition.Ev1Pv2(a.getPosition(), shift);
        
        fileWriter.write((a.getParentGroup().getIndex()+1)+"");
        for (int i=0; i<writePosition.getD(); i++) {
            fileWriter.write(" "+writePosition.getX(i));
        }
        fileWriter.write("\n");
    }

}
