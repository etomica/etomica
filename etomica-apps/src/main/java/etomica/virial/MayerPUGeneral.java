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
import etomica.space.Space;

/**
 * f + u 
 * f: exp(-beta U)
 * u: 0 for r < sigmaHS, dipole energy otherwise
 * a Potential2 object.
 */
public class MayerPUGeneral implements MayerFunction, java.io.Serializable {

    /**
     * Constructor Mayer function using given potential.
     */
	private double mu;
	private IPotentialMolecular potential;
	private double sigmaHS;
	private Vector dr;
    public MayerPUGeneral(Space space, double sigmaHS, double mu, IPotentialMolecular potential) {
    	this.mu=mu;
		this.potential=potential;
		this.sigmaHS=sigmaHS;
		dr=space.makeVector();    }

    public double f(IMoleculeList pair, double r2, double beta) {
        double x = -beta*potential.energy(pair);//-beta U 
        double fbond;
        if (Math.abs(x) < 0.01) {
        	fbond = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        } else {
        	fbond = Math.exp(x)-1;
        }        
        
        double ener = 0;
        
        if (r2 > sigmaHS * sigmaHS){
        
        	IMolecule molecule1 = pair.getMolecule(0);
        	IMolecule molecule2 = pair.getMolecule(1);
        	IAtomList atomList1 = molecule1.getChildList();
        	IAtomList atomList2 = molecule2.getChildList();
        	IAtomOriented atom1 = (IAtomOriented)atomList1.getAtom(0);
        	IAtomOriented atom2 = (IAtomOriented)atomList2.getAtom(0);
        	// should have a loop to loop over all the atoms in the molecules 
        	Vector v1 = atom1.getOrientation().getDirection();//dipole1
        	Vector v2 = atom2.getOrientation().getDirection();//dipole2
        	double cos12= v1.dot(v2);//cos_D1_D2
	
        	dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        	// normalize dr, the vector between the molecules
        	dr.normalize();
	
        	//cos(dipole 1 and r12)
        	double cos_D1_r = v1.dot(dr);
        	//cos(r12 and dipole 2)
        	double cos_r_D2=dr.dot(v2);
        	double r12Magnitude = Math.sqrt(r2);
        	ener  = mu * mu * (cos12 - 3.0  * cos_D1_r * cos_r_D2) / r12Magnitude /r2;
        
        }
        return fbond + beta * ener;
    }

    public IPotential getPotential() {
        return potential;
    }

    public void setBox(Box newBox) {
        potential.setBox(newBox);
    }

}
