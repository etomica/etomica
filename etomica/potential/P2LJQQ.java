/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.api.IVector;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.exception.MethodNotImplementedException;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Coulomb;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.util.RandomNumberGenerator;

/**
 * Lennard Jones molecule with a quadrupole.
 *
 * This is for naphthalene/CO2 mixture virial coefficients calculation
 * Both species are modeled as single site plus QQ 
 * from Albo
 * modified from Jayant K. Singh
 * @author shu
 * Dec.10,2010
 * QQ is orientation independent 
 * the potential is already sphericalization
 */
public class P2LJQQ extends Potential2SoftSpherical  {

    public P2LJQQ(ISpace space) {
        this(space, 1, 1, 1);
    }

    public P2LJQQ(ISpace space, double sigma, double epsilon,  double momentSquared) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setQuadrupolarMomentSquare(momentSquared);
       
    }
      public void setBox(IBox box) {
        boundary = box.getBoundary();
    }
         // LJ contribution
    public double getSigma() {return sigma;}
    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }
    public double getEpsilon() {return epsilon;}

    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4*epsilon;
    }
//??
    public void setQuadrupolarMomentSquare(double moment){
        Q2=moment;
    }
    
    public double getQuadrupolarMomentSquare() {
        return Q2;
    }
     // LJ contribution
    /**
     * Sets the temperature used for Boltzmann-weighting of the orientational
     * average energy used in u(double) and integral(double)
     */
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    // return a Boltzmann-weighted orientational average of the energy
    // for the given distance
    // this is the part I want
    public double u(double r2) {
        double s2 = sigma2/r2;
        if (s2 > 4  ) return Double.POSITIVE_INFINITY ;
        double s6 = s2*s2*s2;
        double r4 = r2*r2;
        return epsilon4*s6*(s6 - 1.0) - (7.0/(5.0*temperature))*Q2*Q2/(r4*r4*r2);
    }
     
    private static final long serialVersionUID = 1L;
    private double sigma , sigma2;
    private double epsilon, epsilon4;
    private double Q2;
    private IBoundary boundary;
    protected double temperature;
	@Override
	public double du(double r2) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double d2u(double r2) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double uInt(double rC) {
		// TODO Auto-generated method stub
		return 0;
	}
}
