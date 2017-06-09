/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.action.WriteConfigurationDLPOLY;
import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.space.Space;

public class PotentialDLPOLY extends PotentialMolecular {
	
	public PotentialDLPOLY(Space space){
		super(0, space);
	}
	
	public double energy(IMoleculeList atoms) {

		configDLPOLY.actionPerformed();
		
		try{
			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec("DLMULTI");
			int exitVal = proc.waitFor();
			FileReader fileReader = new FileReader("ConfigEnergy");			
			BufferedReader bufReader = new BufferedReader(fileReader);
			
			String line = bufReader.readLine();
			return Double.parseDouble(line);
			
			
		}catch (IOException e){
			throw new RuntimeException(e);
		}catch (InterruptedException err){
			throw new RuntimeException(err);
		}
	}

	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}


	public void setBox(Box box) {
		configDLPOLY.setBox(box);
	}
	
	public WriteConfigurationDLPOLY getConfigDLPOLY() {
		return configDLPOLY;
	}

	public void setConfigDLPOLY(WriteConfigurationDLPOLY configDLPOLY) {
		this.configDLPOLY = configDLPOLY;
	}
	
	private WriteConfigurationDLPOLY configDLPOLY;

}
