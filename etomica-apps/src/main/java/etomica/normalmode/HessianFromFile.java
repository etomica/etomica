/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


/**
 * MC simulation of FCC soft-sphere model in 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * @author Tai Boon Tan
 */
public class HessianFromFile {

    public static void main(String[] args) {
    	
        // defaults
        int nA =500;
        String filename = "inputSSDB"+nA;//+"_d0962";
        NormalModesFromFile nm = new NormalModesFromFile(filename, 3);
        CalcHarmonicA.doit(nm, 3, 0.001, nA);
        CalcHarmonicA.doit(nm, 3, 1.2147, nA);
        
        for (int i =1; i<17; i++){
        	double temp = i*0.1;
        	CalcHarmonicA.doit(nm, 3, temp, nA);
        }
    }
    
}