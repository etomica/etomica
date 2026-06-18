/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesGeneral;

import java.io.Serializable;

/**
 * Wave vector factory that returns wave vectors for a 1D system.  
 * These wave vectors are given by 2 Pi m / (N a) = 2 Pi m / L = 2 Pi m rho / N,
 * for m = 1, 2, 3,...,N/2.  If N is odd there will be (N-1)/2 vectors, each with 
 * a coefficient of 1.0 (because none are on the Brillouin-zone boundary) while if 
 * N is even there will be N/2, with the last one having a coefficient of 0.5. Wave
 * vectors corresponding to m <= 0 are not included, although they are technically are among the full
 * set of wave vectors for the system.
 *
 * @author Andrew Schultz and David Kofke
 */
public class WaveVectorFactory1D implements WaveVectorFactory, Serializable {

    public void makeWaveVectors(Box box) {

        int nA = box.getMoleculeList().size();
        double L = box.getBoundary().getBoxSize().getX(0);
        
        int mMax = nA/2;

        waveVectors = new Vector1D[mMax+1];
        coefficients = new double[mMax+1];
        
        for(int m = 0; m<=mMax; m++) {
            waveVectors[m] = new Vector1D(2. * Math.PI * m / L);
            coefficients[m] = 1.0;
        }
        coefficients[0] = 0.5;
        
        if(nA % 2 == 0) {
            coefficients[mMax] = 0.5;
        }

    }
    
    public Vector[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getCoefficients() {
        return coefficients;
    }
    
    public static void main(String[] args) {
        int nCells = 6;
        Space sp = Space1D.getInstance();
        Simulation sim = new Simulation(sp);
        SpeciesGeneral species = SpeciesGeneral.monatomic(sp, AtomType.simpleFromSim(sim));
        sim.addSpecies(species);
        Box box = new Box(sim.getSpace());
        sim.addBox(box);
        box.getBoundary().setBoxSize(new Vector1D(nCells));
        box.setNMolecules(species, nCells);
        
        WaveVectorFactory1D foo = new WaveVectorFactory1D();
        foo.makeWaveVectors(box);
        Vector[] waveVectors = foo.getWaveVectors();
        double[] coefficients = foo.getCoefficients();
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(coefficients[i]+" "+waveVectors[i]);
        }
    }

    private static final long serialVersionUID = 1L;

    protected Vector[] waveVectors;
    protected double[] coefficients;
}
