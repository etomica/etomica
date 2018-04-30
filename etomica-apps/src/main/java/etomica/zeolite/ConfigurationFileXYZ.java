/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.space.Vector;
import etomica.box.Box;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.Configuration;
import etomica.space.Space;
import etomica.util.Resources;

public class ConfigurationFileXYZ implements Configuration, java.io.Serializable {

		public ConfigurationFileXYZ(String aConfName, Space _space){
			super();
			this.space = _space;
			confName = aConfName;
		}
		
		public void initializeCoordinates(Box box) {
            min = space.makeVector();
            max = space.makeVector();
            dim = space.makeVector();
	        String fileName = confName;
            AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms(box);
			try (BufferedReader bufReader = Resources.openResourceFile(confName, this.getClass())) {
				atomIterator.reset();
				//Skips the first line, which contains numAtoms
				bufReader.readLine();
				for (IAtom atom = atomIterator.nextAtom();
					 atom != null; atom = atomIterator.nextAtom()) {
					setPosition(atom, bufReader.readLine());
				}
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			atomIterator.reset();
	        
            for (IAtom atom = atomIterator.nextAtom();
                 atom != null; atom = atomIterator.nextAtom()) {
	        	translatePosition(atom);
	        }
	        
		}
		
		private void setPosition(IAtom atom, String string) {
	        String[] coordStr = string.split(" +");
            Vector pos = atom.getPosition();
            for (int i=0; i<pos.getD(); i++) {
	            double coord = Double.valueOf(coordStr[i]).doubleValue();
	            if(coord<min.getX(i)) min.setX(i,coord);
	            if(coord>max.getX(i)) max.setX(i,coord);
                pos.setX(i, coord);
	        }
            dim.E(max);
            dim.ME(min);
	    }

        private void translatePosition(IAtom atom){
            atom.getPosition().ME(min);
            atom.getPosition().PEa1Tv1(-0.5,dim);
		}
        
		public int[] getNumAtoms(){
			String numLine;
			try (BufferedReader bufReader = Resources.openResourceFile(confName, this.getClass())) {
				numLine = bufReader.readLine();
			} catch (IOException e) {
			    throw new RuntimeException(e);
			}
			String[] coordStr = numLine.split(" +");
	        double[] coord = new double[coordStr.length];
	        nAtomsList = new int[coordStr.length];
	        for (int i=0; i<coord.length; i++) {
	            nAtomsList[i] = Integer.parseInt(coordStr[i]);
	        }
	        return nAtomsList;
		}
		
		public Vector getUpdatedDimensions(){
			updatedDimensions = space.makeVector();
			updatedDimensions.E(dim);
			updatedDimensions.TE(1.01);
			return updatedDimensions;
		}
		
        private static final long serialVersionUID = 2L;
		private Vector updatedDimensions;
		private Vector min;
		private Vector max;
		private Vector dim;
		private int[] nAtomsList;
		private String confName;
		private final Space space;

}
