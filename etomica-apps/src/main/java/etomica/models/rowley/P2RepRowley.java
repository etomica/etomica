/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Energy;

/**
 * Purely repulsive potential from Rowley et al (2006) used for interactions between satellite sites, X.
 * These fudge sites are used to represent a region of high electron density belonging to the oxygen atom of simple alcohols.    
 *
 * K.R. Schadel
 * May 2008
 */
public final class P2RepRowley extends Potential2SoftSpherical {

    public P2RepRowley (Space space) {
        this(space, 1.0, 1.0);
    }
    
    public P2RepRowley (Space space, double BXX, double CXX) {
        super(space);
        setBXX(BXX);
        setCXX(CXX);

    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
    	
    	return BXX*Math.exp(-CXX*r);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	
    	double r = Math.sqrt(r2);
    	
    	double first_derivative = -BXX*CXX*Math.exp(-CXX*r);
    	
    	return r*first_derivative;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
    	
    	double r = Math.sqrt(r2);
    	
    	double second_derivative = BXX*CXX*CXX*Math.exp(-CXX*r);
 	
    	return second_derivative*r2;
    }
            
    /**
     * Integral of the potential, used to evaluate corrections for potential truncation.
     * Specifically, this is the integral from rC (the argument) to infinity of
     * u(r) A r^(D-1), where D is the spatial dimension, and A is the area of a unit
     * sphere in D dimensions.  Normally, the long-range potential correction would be obtained
     * by multiplying this quantity by the pair density nPairs/V, where nPairs is the number of pairs of atoms
     * affected by this potential, and V is the volume they occupy.
     */
    public double uInt(double rC) {
    	
    	double A = space.sphereArea(1.0);
    	
        return A*BXX*Math.exp(-CXX*rC)*(rC*rC/CXX + 2*rC/(CXX*CXX) + 2/(CXX*CXX*CXX));  
        
    }


    public double getBXX() {return BXX;}
 
    public final void setBXX(double eps) {
    	BXX = eps;
    }
    public Dimension getBXXDimension() {return Energy.DIMENSION;}
   
    
    public double getCXX() {return CXX;}

    public final void setCXX(double a) {
    	CXX = a;
    }
    
    private static final long serialVersionUID = 1L;
    private double BXX;
    private double CXX;
 
}

