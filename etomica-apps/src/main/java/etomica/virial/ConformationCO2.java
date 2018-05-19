/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.units.Electron;

 /**
  *  Conformation for CO2
  *  Reference paper: Partial molar ~~, Stubbs Darke-Wilhelm, Siepmann, JCP
  * 
 * @author Shu Yang
 *
 */
public class ConformationCO2 implements IConformation, java.io.Serializable{
	
	public ConformationCO2(Space space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
			
		IAtom n1 = atomList.get(SpeciesTraPPECO2.indexC);
		n1.getPosition().E(new double[] {0, 0, 0});
		
		IAtom n2 = atomList.get(SpeciesTraPPECO2.indexOleft);
		n2.getPosition().E(new double[] {-bondlength, 0, 0});
		
		IAtom n3 = atomList.get(SpeciesTraPPECO2.indexOright);
		n3.getPosition().E(new double[] {bondlength, 0, 0});
		
		
	}
	
    public final static double [] Echarge = new double [3];
    static {
    	
    	//add + charge to C, - charge to O, what is the unit of the pt charge?
    	
        ConformationCO2.Echarge[SpeciesTraPPECO2.indexC] = Electron.UNIT.toSim( 0.70);
        
        ConformationCO2.Echarge[SpeciesTraPPECO2.indexOleft] = Electron.UNIT.toSim( -0.75);
        ConformationCO2.Echarge[SpeciesTraPPECO2.indexOright] = Electron.UNIT.toSim( -0.75);

        }
	
	protected final Space space;
	protected static final double bondlength = 1.16;

	
	protected Vector vector;
	
	private static final long serialVersionUID = 1L;
}
