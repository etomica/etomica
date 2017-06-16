/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

/**
 * Provides information about all the normal modes for a periodic system.  The periodicity is
 * described by a set of wave vectors, which are provided via a WaveVectorFactory.  For
 * each wave vector there is a set of coupled degrees of freedom (e.g., xyz motions) that have been 
 * further decomposed into normal modes with corresponding frequencies and eigenvectors.  
 * 
 * In most cases these frequencies/vectors are determined via a simulation and recorded to file; the
 * implmentation of this class then reads the data and provides it via the interface methods.     
 */
public interface NormalModes {

    /**
     * Returns a factory that provides the wave vectors and coefficients for 
     * the periodic system.
     */
    public WaveVectorFactory getWaveVectorFactory();

    /**
     * Returns an array giving the frequencies (squared) corresponding to the normal-mode
     * motions. First index indicates the wave vector, and second index indicates the
     * eigenvector. Length of second index is coordinateDim.
     */
    public double[][] getOmegaSquared();

    /**
     * First index corresponds to the wave vector; second index gives the eigenvector, and
     * third index gives the elements of the vector.  Length of 2nd and 3rd dimensions is 
     * coordinateDim. 
     */
    public double[][][] getEigenvectors();
    
    /**
     * Set the fudge factor applied to frequencies.  The squared-frequencies returned
     * by getOmegaSquared will be the nominal values divided by the given
     * fudge factor.  Thus a smaller value of this fudge factor will make for "tighter"
     * harmonic springs, corresponding to to smaller deviation from the lattice sites.
     */
    public void setHarmonicFudge(double newHarmonicFudge);
    
    /**
     * Set the temperature used for calculating omega-squared.  Omega-squared
     * is kT/eigenvalue.
     */
    public void setTemperature(double newTemperature);
}
