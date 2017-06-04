/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * A combination of the basic LEPS potential describing 3 atoms moving in 1-dimension, paired with a harmonic potential.
 * @author msellers
 *
 */
public class P1LepsHarmonic extends Potential1 implements PotentialSoft {
	
    private static final long serialVersionUID = 1L;

    private final Vector[] force;
	
    //LEPS Potential Parameters
    public double alpha = 1.942;
    public double dab = 4.746;
    public double dbc = 4.746;
    public double dac = 3.445;
    public double rzero = 0.742;
    public double rac = 3.742;
    public double a = 0.05, b = 0.80, c = 0.05;
    
    //Harmonic Potential Parameters
    public double kc = 0.2025;
    public double cs = 1.154;
    
    //Gaussian Parameters
    public double A1 = 1.5;
    public double A2 = 6.0;
    public double xz1 = 2.02083;
    public double xz2 = 0.8;
    public double yz1 = -0.172881;
    public double yz2 = 2.0;
    public double sx1 = 0.1;
    public double sx2 = 5.0;
    public double sy1 = 0.35;
    public double sy2 = 0.7;
    
	public P1LepsHarmonic(Space space){
		super(space);
	    force = new Vector[]{space.makeVector()};
	}
	
	public double energy(IAtomList atom) {
		double x = atom.getAtom(0).getPosition().getX(0);
		double y = atom.getAtom(0).getPosition().getX(1);
		
		double energy;
		
		double first = 4.0*A1*Math.exp(-((x-xz1)*(x-xz1))/(2.0*sx1)-((y-yz1)*(y-yz1))/(2.0*sy1)) + 4.0*A2*Math.exp(-((x-xz2)*(x-xz2))/(2.0*sx2)-((y-yz2)*(y-yz2))/(2.0*sy2))
				+ dac*Math.exp(alpha*(-rac+rzero))*(-2.0+3.0*Math.exp(alpha*(-rac+rzero)))/(1.0+c) + dab*Math.exp(alpha*(rzero-x))*(-2.0+3.0*Math.exp(alpha*(rzero-x)))/(1.0+a)
				+ dbc*Math.exp(alpha*(-rac+rzero+x))*(-2.0+3.0*Math.exp(alpha*(-rac+rzero+x)))/(1.0+b);
		
		double s1, s2, s3;	
		
			s1 = dac*dac*Math.exp(2.0*alpha*(-rac+rzero))*(-6.0+Math.exp(alpha*(-rac+rzero)))*(-6.0+Math.exp(alpha*(-rac+rzero)))/(1.0+c)/(1.0+c) 
				+ dab*dab*Math.exp(2.0*alpha*(rzero-x))*(-6.0+Math.exp(alpha*(rzero-x)))*(-6.0+Math.exp(alpha*(rzero-x)))/(1.0+a)/(1.0+a)
				- dab*dac*Math.exp(-2.0*alpha*(rac-(2.0*rzero)+x))*(-1.0+6.0*Math.exp(alpha*(rac-rzero)))*(-1.0+6.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+c);
			
			s2 = dab*dbc*Math.exp(-2.0*alpha*(rac-(2.0*rzero)))*(-1.0+6.0*Math.exp(alpha*(rac-rzero-x)))*(-1.0+6.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+b)
				- dac*dbc*Math.exp(alpha*(-4.0*rac+(2.0*rzero)+x))*(-6.0*Math.exp(alpha*rac)+Math.exp(alpha*rzero))*(-6.0*Math.exp(alpha*rac)+Math.exp(alpha*(rzero+x)))/(1.0+b)/(1.0+c);
			s3 = dbc*dbc*Math.exp(2.0*alpha*(-rac+rzero+x))*(-6.0+Math.exp(alpha*(-rac+rzero+x)))*(-6.0+Math.exp(alpha*(-rac+rzero+x)))/(1.0+b)/(1.0+b);
		
		double second = Math.sqrt(s1-s2+s3);
		
		energy = first - second	+ 8.0*kc*(-rac/2.0 + y/cs + x)*(-rac/2.0 + y/cs + x);
		
		energy =  0.25 * energy;
		
		return energy;
	}
		
	public Vector[] gradient(IAtomList atom){

		double x = atom.getAtom(0).getPosition().getX(0);
		double y = atom.getAtom(0).getPosition().getX(1);
		double dudx;
		double dudy;
		
		double one =  alpha*dab*Math.exp(2*alpha*(rzero-x))*(-3.0+Math.exp(alpha*(-rzero+x)))/(2.0*(1+a)) 
					 + alpha*dbc*Math.exp(alpha*(-rac+rzero+x))*(-1.0+3.0*Math.exp(alpha*(-rac+rzero+x)))/(2.0*(1+b));
		
		double num =  4.0*alpha*dab*dab*Math.exp(2.0*alpha*(rzero-x))*(-6.0+Math.exp(alpha*(rzero-x)))*(-3.0+Math.exp(alpha*(rzero-x)))/(1.0+a)/(1.0+a)       
					- 2.0*alpha*dab*dac*Math.exp(-2.0*alpha*(rac-(2.0*rzero)+x))*(-1.0+6.0*Math.exp(alpha*(rac-rzero)))*(-1.0+3.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+c)
					- 2.0*alpha*dab*dbc*Math.exp(-2.0*alpha*(rac-(2.0*rzero)))*(-1.0+6.0*Math.exp(alpha*(rac-rzero-x)))*(-1.0+3.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+b)
					+ 2.0*alpha*dab*dbc*Math.exp(-2.0*alpha*(rac-(2.0*rzero)))*(-1.0+3.0*Math.exp(alpha*(rac-rzero-x)))*(-1.0+6.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+b)
					+ 2.0*alpha*dac*dbc*Math.exp(alpha*(-4.0*rac+(2.0*rzero)+x))*(-6.0*Math.exp(alpha*rac)+Math.exp(alpha*rzero))*(-3.0*Math.exp(alpha*rac)+Math.exp(alpha*(rzero+x)))/(1.0+b)/(1.0+c)
					- 4.0*alpha*dbc*dbc*Math.exp(2.0*alpha*(-rac+rzero+x))*(-6.0+Math.exp(alpha*(-rac+rzero+x)))*(-3.0+Math.exp(alpha*(-rac+rzero+x)))/(1.0+b)/(1.0+b);
				 
				
		double denom =	dac*dac*Math.exp(2.0*alpha*(-rac+rzero))*(-6.0+Math.exp(alpha*(-rac+rzero)))*(-6.0+Math.exp(alpha*(-rac+rzero)))/(1.0+c)/(1.0+c)
						 + dab*dab*Math.exp(2.0*alpha*(rzero-x))*(-6.0+Math.exp(alpha*(rzero-x)))*(-6.0+Math.exp(alpha*(rzero-x)))/(1.0+a)/(1.0+a)
						 - dab*dac*Math.exp(-2.0*alpha*(rac-(2.0*rzero)+x))*(-1.0+6.0*Math.exp(alpha*(rac-rzero)))*(-1.0+6.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+c)
						 - dab*dbc*Math.exp(-2.0*alpha*(rac-(2.0*rzero)))*(-1.0+6.0*Math.exp(alpha*(rac-rzero-x)))*(-1.0+6.0*Math.exp(alpha*(-rzero+x)))/(1.0+a)/(1.0+b)  
						 - dac*dbc*Math.exp(alpha*(-4.0*rac+(2.0*rzero)+x))*(-6.0*Math.exp(alpha*rac)+Math.exp(alpha*rzero))*(-6.0*Math.exp(alpha*rac)+Math.exp(alpha*(rzero+x)))/(1.0+b)/(1.0+c) 
						 + dbc*dbc*Math.exp(2.0*alpha*(-rac+rzero+x))*(-6.0+Math.exp(alpha*(-rac+rzero+x)))*(-6.0+Math.exp(alpha*(-rac+rzero+x)))/(1.0+b)/(1.0+b);
				
						 
		double three = (A1*Math.exp((-(x-xz1)*(x-xz1)/(2.0*sx1))-((y-yz1)*(y-yz1)/(2.0*sy1)))*(-x+xz1))/sx1 + (A2*Math.exp((-(x-xz2)*(x-xz2)/(2.0*sx2))-((y-yz2)*(y-yz2)/(2.0*sy2)))*(-x+xz2))/sx2
				+ 4.0*kc*((-rac/2.0) + x + y/cs);
		
		dudx = one + (num / (8.0*Math.sqrt(denom))) + three;
		
		dudy =  (-2.0*kc*(cs*(rac-2.0*x)-2.0*y))/cs/cs + (A1*Math.exp(-(x-xz1)*(x-xz1)/(2.0*sx1) - (y-yz1)*(y-yz1)/(2.0*sy1))*(-y+yz1))/sy1 
				+ (A2*Math.exp(-(x-xz2)*(x-xz2)/(2.0*sx2) - (y-yz2)*(y-yz2)/(2.0*sy2))*(-y+yz2))/sy2;
		
		force[0].setX(0,dudx);
		force[0].setX(1,dudy);
		
		return force;
	}
	
	public Vector[] gradient(IAtomList atom, Tensor pressureTensor){
		
		return gradient(atom);
	}
	
	public double virial(IAtomList atom){
		
		return 0.0;
	}

}
