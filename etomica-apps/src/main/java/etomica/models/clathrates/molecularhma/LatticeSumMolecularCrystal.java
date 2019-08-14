/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.clathrates.molecularhma;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.*;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.spaceNd.TensorND;

/**
 *
 * calculates dynamical matrix for rigid molecular crystal
 */
public class LatticeSumMolecularCrystal {

    public LatticeSumMolecularCrystal(PotentialMaster potentialMaster, Box box, final Space space, int basisDim) {

    	this.potentialMaster = potentialMaster;
    	this.box = box;
    	this.space = space;
    	this.basisDim = basisDim;
    	this.atomPosDef = new MoleculePositionGeometricCenterPBC(space, box.getBoundary());
    	this.com0 = space.makeVector();
    	this.com1 = space.makeVector();
    	int atomsPerMol = box.getLeafList().size() / box.getMoleculeList().size();
    	this.tmpAtomicTensor3 = new Tensor[basisDim*atomsPerMol];//46X4=184 atoms dimensional array of 3D Tensor
    	for (int i=0; i<basisDim*atomsPerMol ;i++){
    		tmpAtomicTensor3[i] = space.makeTensor();
    	}
		pcForce = new PotentialCalculationForceSum();
		atomAgentManager = new AtomLeafAgentManager<>(a -> space.makeVector() , box);
		pcForce.setAgentManager(atomAgentManager);
         IteratorDirective id = new IteratorDirective();
        id.includeLrc = false;
        potentialMaster.calculate(box, id, pcForce);
		tmpDrr1 = space.makeTensor();
    }

	public void reset() {
		IteratorDirective id = new IteratorDirective();
		id.includeLrc = false;
		potentialMaster.calculate(box, id, pcForce);
	}
    public Tensor atomicToMolecularD(AtomicTensorAtomicPair aTensor, IMolecule mol0, IMolecule mol1){
//		System.out.println("I'm in latticesummolecularclass");

    	TensorND D6mol = new TensorND(6);
    	Tensor3D D3tt  = new Tensor3D();	Tensor3D D3tr  = new Tensor3D();
    	Tensor3D D3rt  = new Tensor3D();	Tensor3D D3rr  = new Tensor3D();
    	Tensor3D D3tt_ = new Tensor3D(); Tensor3D D3tr_ = new Tensor3D();
    	Tensor3D Rk = new Tensor3D();  Tensor3D Rk_ = new Tensor3D();
    	Tensor3D Rkp = new Tensor3D(); 
    	Vector Xk = space.makeVector();
    	Vector Xkp = space.makeVector();
    	MoleculePositionCOM com_0 = new MoleculePositionCOM(space);
    	Vector com0 = com_0.position(mol0);
    	MoleculePositionCOM com_1 = new MoleculePositionCOM(space);
    	Vector com1 = com_1.position(mol1);

		int numSites0 = mol0.getChildList().size();
		int numSites1 = mol1.getChildList().size();
    	for (int atomk=0; atomk < numSites0; atomk++){
    		Vector posk = mol0.getChildList().get(atomk).getPosition();
    		Xk.Ev1Mv2(posk, com0);
    		box.getBoundary().nearestImage(Xk);
    		Rk.setComponent(0,0,0.0); Rk.setComponent(1,1,0.0);	Rk.setComponent(2,2,0.0);
    		Rk.setComponent(0,1, Xk.getX(2));  Rk.setComponent(1,0,-Xk.getX(2));   
    		Rk.setComponent(0,2,-Xk.getX(1));  Rk.setComponent(2,0, Xk.getX(1));
    		Rk.setComponent(1,2, Xk.getX(0));  Rk.setComponent(2,1,-Xk.getX(0));
    		for (int atomkp=0; atomkp < numSites1; atomkp++){
    			if(atomk == atomkp && mol0 == mol1) continue;//ADDED:: Non-self
        		Vector poskp = mol1.getChildList().get(atomkp).getPosition();
        		Xkp.Ev1Mv2(poskp, com1);
        		box.getBoundary().nearestImage(Xkp);
        		Rkp.setComponent(0,0,0.0); Rkp.setComponent(1,1,0.0);	Rkp.setComponent(2,2,0.0);
        		Rkp.setComponent(0,1, Xkp.getX(2));  Rkp.setComponent(1,0,-Xkp.getX(2));      
        		Rkp.setComponent(0,2,-Xkp.getX(1));  Rkp.setComponent(2,0, Xkp.getX(1));
        		Rkp.setComponent(1,2, Xkp.getX(0));  Rkp.setComponent(2,1,-Xkp.getX(0));


        		D3tt_.E(aTensor.atomicTensor(mol0.getChildList().get(atomk) , mol1.getChildList().get(atomkp)));
        		D3tt.PE(D3tt_);
        		D3tr_.E(D3tt_); D3tr_.TE(Rkp);  D3tr.PE(D3tr_);
        		Rk_.E(Rk);  Rk_.TE(D3tt_);  D3rt.ME(Rk_);
                Rk_.TE(Rkp); D3rr.ME(Rk_);    
        		tmpAtomicTensor3[mol0.getChildList().get(atomk).getLeafIndex()].ME(D3tt_);//self summation(sum rule)
        	}//atomkp
    		
    		if(mol0 == mol1){//self
    			D3tt_.E(tmpAtomicTensor3[mol0.getChildList().get(atomk).getLeafIndex()]);
    			
        		D3tt.PE(D3tt_);//Transform to molecular

        		D3tr_.E(D3tt_);	//take atomic and multiply by Rk
        		D3tr_.TE(Rk);  
        		D3tr.PE(D3tr_);//add all atomic to Transform to molecular
        		
        		Rk_.E(Rk);
        		Rk_.TE(D3tt_); //D3tt_=  - SUM(k =/= k')
        		D3rt.ME(Rk_);
//Drr Original        		
        		Rk_.TE(Rk);
        		D3rr.ME(Rk_);
//Drr Symmetric Part
 
                Tensor3D tmpDrr2 = new Tensor3D();
	        	IAtom atom = mol0.getChildList().get(atomk);
	        	
	        	Vector fk = atomAgentManager.getAgent(atom);//gradient NOT fk
	        	Xk.TE(-1);
	            tmpDrr1.Ev1v2(Xk, fk);
	            if(true){ //Symmetrize? SUM_over_atoms{tmpDrr1} should = ZERO; i.e. torque=0
		            tmpDrr2.E(tmpDrr1);
		            tmpDrr1.transpose();
		            tmpDrr1.PE(tmpDrr2); 
		            tmpDrr1.TE(0.5); 
	            }
	            tmpDrr1.setComponent(0, 0,  -fk.getX(1) * Xk.getX(1) - fk.getX(2) * Xk.getX(2)); 
	            tmpDrr1.setComponent(1, 1,  -fk.getX(0) * Xk.getX(0) - fk.getX(2) * Xk.getX(2)); 
	            tmpDrr1.setComponent(2, 2,  -fk.getX(0) * Xk.getX(0) - fk.getX(1) * Xk.getX(1)); 
	            D3rr.PE(tmpDrr1);
    		}//mol0 == mol1 ?
    	}//atomk

    	for (int i=0; i<3; i++){
        	for (int j=0; j<3; j++){
        		D6mol.setComponent(i,  j , D3tt.component(i, j));
        		D6mol.setComponent(i,  j+3, D3tr.component(i, j));
        		D6mol.setComponent(i+3,j, D3rt.component(i, j));
        		D6mol.setComponent(i+3,j+3, D3rr.component(i, j));
        	}
    	}
	//	for (int i=0; i<6; i++) {
	//		for (int j = 0; j < 6; j++) {
	//			System.out.println( D6mol.component(4, 4) );
	//		}
	//	}
		return D6mol;
    }
	PotentialCalculationForceSum pcForce;
    private final int basisDim;
    protected final Box box;
    protected final Space space;
    protected IMoleculePositionDefinition atomPosDef;
    protected final Vector com0, com1;
     protected final Tensor[] tmpAtomicTensor3;//46X4=184 atoms dimentional array of 3dim Tensor
	protected PotentialMaster potentialMaster;
	protected AtomLeafAgentManager<Vector> atomAgentManager;
	protected   Tensor tmpDrr1;

    public interface AtomicTensorAtomicPair{
    	public Tensor atomicTensor(IAtom atom0, IAtom atom1);
    }
}
