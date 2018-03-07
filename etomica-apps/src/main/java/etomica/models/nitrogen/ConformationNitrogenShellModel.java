/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.units.Electron;

 /**
  *  Conformation for Nitrogen (Shell Model)
  *  
  *	 Reference: Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
  *             phases, JCP 112(15) 6745 (2000)
  * 
  * @author Tai Boon Tan
  *
  */
public class ConformationNitrogenShellModel implements IConformation, java.io.Serializable{
	
	public ConformationNitrogenShellModel(Space space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
		
		IAtom n1 = atomList.get(SpeciesN2ShellModel.indexN1);
		n1.getPosition().E(new double[] {-bondOrigN, 0, 0});
		
		IAtom n2 = atomList.get(SpeciesN2ShellModel.indexN2);
		n2.getPosition().E(new double[] {bondOrigN, 0, 0});
			
		IAtom center = atomList.get(SpeciesN2ShellModel.indexCenter);
		center.getPosition().E(new double[] {0, 0, 0});
		
		IAtom p1Left = atomList.get(SpeciesN2ShellModel.indexP1left);
		p1Left.getPosition().E(new double[] {-bondOrigP, 0, 0});
		
		IAtom p1Right = atomList.get(SpeciesN2ShellModel.indexP1right);
		p1Right.getPosition().E(new double[] {bondOrigP, 0, 0});
		
	}
	
    public final static double [] Echarge = new double [5];
    static {
    
        ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexN1] = Electron.UNIT.toSim( 2.563085);
        ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexN2] = Electron.UNIT.toSim( 2.563085);
    	ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexCenter] = Electron.UNIT.toSim( -2.49737);
        ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexP1left] = Electron.UNIT.toSim(-1.3144);
        ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexP1right] = Electron.UNIT.toSim(-1.3144);

    }
	
	protected final Space space;
	protected static final double bondOrigN = 0.5485;
	protected static final double bondOrigP = 0.81; // 0.83519;
	
	protected Vector vector;
	
	private static final long serialVersionUID = 1L;
}
