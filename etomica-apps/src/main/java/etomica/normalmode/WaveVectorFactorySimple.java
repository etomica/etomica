/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;

import java.io.Serializable;

/**
 * Wave vector factory that returns wave vectors appropriate for a box with 
 * a single-atom basis.  The box shape need not be rectangular so long as it
 * matches the primitive's shape.
 *
 * @author Andrew Schultz
 */
public class WaveVectorFactorySimple implements WaveVectorFactory, Serializable {

    public WaveVectorFactorySimple(Primitive primitive, Space _space) {
        this.primitive = primitive;
        this.space = _space;
    }
    
    public void makeWaveVectors(Box box) {
        // If we weren't given wave vectors, determine them from the box boundary and primitve
        // assume 1-molecule basis and matchup betwen the box and the primitive
    
        double[] d = primitive.getSize();
        int[] numCells = new int[space.D()];
        Vector[] reciprocals =  primitive.makeReciprocal().vectors();
        Vector[] waveVectorBasis = new Vector[reciprocals.length];
                
        for (int i=0; i<space.D(); i++) {
            waveVectorBasis[i] = space.makeVector();
            waveVectorBasis[i].E(reciprocals[i]);
            numCells[i] = (int)Math.round(Math.sqrt(box.getBoundary().getEdgeVector(i).squared()) / (d[i]));
            waveVectorBasis[i].TE(1.0/numCells[i]);
        }

        int[] kMin = new int[space.D()];
        int[] kMax= new int[space.D()];
        for (int i=0; i<kMax.length; i++) {
            kMin[i] = -(numCells[i]-1)/2;
            kMax[i] = numCells[i]/2;
        }
        int[][][] waveVectorIndices = new int[2*kMax[0]+1][2*kMax[1]+1][2*kMax[2]+1];
        int count = 0;
        int[] idx = new int[3];
        // if we have an odd number of cells, we flip the first non-zero component
        // positive.  If we have an even number of cells, we do so, but ignore any
        // component equal to the max.  And, for even number of cells, if we flip,
        // we re-flip any component that flips to -max. 
        
        boolean [] flip2 = new boolean [3];
        
        for (int i=0; i<3; i++){
        	flip2[i] = numCells[i] %2 == 0;
        }
        
        // this will find N-1 vectors.  Some of them have negatives 
        // within the set others do not.  If its negative is within the set, 
        // exclude the negative, but remember it was there -- they will have 
        // coefficients of '1' while the ones without a negative in the set 
        // will have coefficients of '0.5'.
        int [] k = new int [3];
        
        for ( k[0] = kMin[0]; k[0] < kMax[0]+1; k[0]++) {
            for ( k[1] = kMin[1]; k[1] < kMax[1]+1; k[1]++) {
                for (k[2] = kMin[2]; k[2] < kMax[2]+1; k[2]++) {
//                    if (kx == 0 && ky == 0 && kz == 0) continue;
                    for (int i=0; i<3; i++) {
                        idx[i] = kMax[i];
                    }
                    
                    boolean flip = false;
outer:              for (int i=0; i<3; i++){
                    	for (int j=0; j<i; j++){
                    		if (k[j] > 0 && (!flip2[j] || k[j] < kMax[j] )){
                    			break outer;
                    		}
                    	}
                    	if (k[i] < 0){
                    		flip = true;
                    		break;
                    	}
                    }
                    
                    if (flip) {
                        idx[0] -= k[0];
                        idx[1] -= k[1];
                        idx[2] -= k[2];
                    }
                    else {
                        idx[0] += k[0];
                        idx[1] += k[1];
                        idx[2] += k[2];
                    }
                    
                    if (flip) {
                        for (int i=0; i<3; i++) {
                            if (idx[i] == 0 && flip2[i]) {
                                idx[i] = 2*kMax[i];
                            }
                        }
                    }

                    if (waveVectorIndices[idx[0]][idx[1]][idx[2]] == 0) {
                        // this one was unique
                        count++;
                    }
                    waveVectorIndices[idx[0]][idx[1]][idx[2]]++;
                }
            }
        }
        waveVectors = new Vector3D[count];
        coefficients = new double[count];
        count = 0;
        for (int kx = -kMax[0]; kx < kMax[0]+1; kx++) {
            for (int ky = -kMax[1]; ky < kMax[1]+1; ky++) {
                for (int kz = -kMax[2]; kz < kMax[2]+1; kz++) {
                    if (waveVectorIndices[kx+kMax[0]][ky+kMax[1]][kz+kMax[2]] > 0) {
                        waveVectors[count] = space.makeVector();
                        waveVectors[count].Ea1Tv1(kx, waveVectorBasis[0]);
                        waveVectors[count].PEa1Tv1(ky, waveVectorBasis[1]);
                        waveVectors[count].PEa1Tv1(kz, waveVectorBasis[2]);
                        coefficients[count] = waveVectorIndices[kx+kMax[0]][ky+kMax[1]][kz+kMax[2]]/2.0;
                        count++;
                    }
                }
            }
        }
    }
    
    public Vector[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getCoefficients() {
        return coefficients;
    }
    
    public static void main(String[] args) {
        int [] nCells = new int []{2,3,4};
        Space sp = Space3D.getInstance();
        Simulation sim = new Simulation(sp);
        SpeciesGeneral species = SpeciesGeneral.monatomic(sp, AtomType.simpleFromSim(sim));
        sim.addSpecies(species);
        Box box = new Box(sp);
        sim.addBox(box);
        box.getBoundary().setBoxSize(new Vector3D(nCells[0], nCells[1], nCells[2]));
        box.setNMolecules(species, nCells[0]*nCells[1]*nCells[2]);
        Primitive primitive = new PrimitiveCubic(sim.getSpace(), 1);
        
        WaveVectorFactorySimple foo = new WaveVectorFactorySimple(primitive, sp);
        foo.makeWaveVectors(box);
        Vector[] waveVectors = foo.getWaveVectors();
        double[] coefficients = foo.getCoefficients();
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(coefficients[i]+" "+waveVectors[i]);
        }
    }

    private static final long serialVersionUID = 1L;
    protected final Primitive primitive;
    protected Vector[] waveVectors;
    protected double[] coefficients;
    protected final Space space;
}
