/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceScalar;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.space.ISpace;
import etomica.units.Undefined;

/**
 * Meter for calculation of structure factor of a group of atoms based on a
 * particular wave vector.  The meter is initially constructed using the smallest 
 * reciprocal lattice vector as the wave vector and loops over all atoms in
 * the box.  GetData returns the square of the magnitude of the structure factor.
 *
 * @author Michael Sellers
 */

public class MeterStructureFactor extends DataSourceScalar {
	
    private static final long serialVersionUID = 1L;
    protected BravaisLatticeCrystal lattice;
	protected final ISpace space;
    protected IBox box;
    protected double[] struct;
    protected IVectorMutable [] waveVec;
    protected IMoleculeList moleculeList;

    
    
	
	/**
	 * Creates meter with default to compute the structure factor
	 * for all atoms in the box.
	 * @param parent
	 */
	public MeterStructureFactor(ISpace space, BravaisLatticeCrystal aLattice, IBox aBox){
		super("Structure factor", Undefined.DIMENSION);
		this.space = space;
        
        this.lattice = aLattice;
        this.box = aBox;
        
        waveVec = new IVectorMutable[lattice.getPrimitive().makeReciprocal().vectors().length];
        for(int i=0; i<waveVec.length; i++){
        	waveVec[i] = space.makeVector();
        	waveVec[i].E(lattice.getPrimitive().makeReciprocal().vectors()[i]);
        }
      
        struct = new double[waveVec.length];
        moleculeList = box.getMoleculeList();
	}
	
	public void reset() {
		
		for(int i=0; i<waveVec.length; i++){
        	waveVec[i].setX(0, 0.0);
        	waveVec[i].setX(1, 0.0);
        	waveVec[i].setX(2, 0.0);
        	struct[i] = 0.0;
		}
		
		moleculeList = box.getMoleculeList();
	}
	
	public void actionPerformed() {
		double term1 = 0;
		double term2 = 0;
		double dotprod = 0;
		int numAtoms = moleculeList.getMoleculeCount();
		IVectorMutable workvector = space.makeVector();
		for(int k=0; k<waveVec.length; k++){
			term1 = 0;
			term2 = 0;
			dotprod = 0;
			for(int i=0; i<numAtoms; i++){
				workvector.E(moleculeList.getMolecule(i).getChildList().getAtom(0).getPosition());
				dotprod = waveVec[k].dot(workvector);
				term1 += Math.cos(dotprod); 
				term2 += Math.sin(dotprod);
			}
			struct[k] = ((term1*term1) + (term2*term2))/(numAtoms*numAtoms);
		}
	}
	
	/**
	 * @param waveVec Sets a custom wave vector array.
	 */
	public void setWaveVec(IVectorMutable [] _waveVec){
		 waveVec = space.makeVectorArray(_waveVec.length);
		 struct = new double[waveVec.length];
		for(int i=0; i<_waveVec.length; i++){
			waveVec[i].E(_waveVec[i]);
		}
		
	}
	
	/**
	 * @param atomList Sets the list of atoms for factor calculation.
	 */
	public void setAtoms(IMoleculeList moleculeList){
		this.moleculeList = moleculeList;
	}

	//this is really the square of the magnitude of the structure factor
    public double[] getDataAsArray() {
    	return struct;
    }
    	
	public double getDataAsScalar() {
		for(int i=0; i<struct.length-1; i++){
			struct[struct.length-1] += struct[i];
		}
		return struct[struct.length-1];
	}
	
}
