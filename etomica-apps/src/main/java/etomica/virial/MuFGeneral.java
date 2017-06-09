/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.api.IPotential;
import etomica.api.IPotentialMolecular;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;

/**
 * @author shu
 *
 * mu_i * mu_j * General Mayer f function
 */
public class MuFGeneral implements MayerFunction, java.io.Serializable {
	
	private double mu;
	private IPotentialMolecular potential;
	public MuFGeneral(double mu, IPotentialMolecular potential) {
		this.mu=mu;
		this.potential=potential;
	}

	public double f(IMoleculeList pair, double r2, double beta) {
	    double x = -beta*potential.energy(pair);//-betaU

		IMolecule molecule1 = pair.getMolecule(0);
		IMolecule molecule2 = pair.getMolecule(1);
		IAtomList atomList1 = molecule1.getChildList();
		IAtomList atomList2 = molecule2.getChildList();
        IAtomOriented atom1 = (IAtomOriented)atomList1.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)atomList2.getAtom(0);
		// should have a loop to loop over all the atoms in the molecules 
        Vector v1 = atom1.getOrientation().getDirection();//dipole1
        Vector v2 = atom2.getOrientation().getDirection();//dipole2
        double cos12= v1.dot(v2);
		double mu_Dot_mu= mu*mu*cos12;
			
    	if ( x== Double.NEGATIVE_INFINITY){
    		return -mu_Dot_mu;
    	}
	    double fFunction = 0.0;
        if (Math.abs(x) < 0.01) {
            fFunction =  x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        } else{
        	fFunction =  Math.exp(x) - 1.0;
        }
        double sum = fFunction*mu_Dot_mu;
		//System.out.println("cos (mu1) and (mu2): " + cos12);
		//System.out.println("(mu) : " + mu);
		return sum;
			
	}

	/* (non-Javadoc)
	 * @see etomica.virial.MayerFunction#getPotential()
	 */  
	public IPotential getPotential() {
		// TODO Auto-generated method stub
		return potential;
	}
	
	public void setBox(Box newBox) {
	    potential.setBox(newBox);
	}
}
