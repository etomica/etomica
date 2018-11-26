/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.lattice.crystal.Primitive;
import etomica.molecule.*;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.spaceNd.TensorND;

public class LatticeSumMolecularCrystal {

    public LatticeSumMolecularCrystal(PotentialMaster potentialMaster, Box box, final Space space, int basisDim, Primitive primitive) {

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
        kFactory = new WaveVectorFactorySimple(primitive, space);
        kFactory.makeWaveVectors(box);
        double[] kCoefficients = kFactory.getCoefficients(); //kCoefficients=0.5 non-deg.; = 1 degenerate twice!

		PotentialCalculationForceSum pcForce = new PotentialCalculationForceSum();
		atomAgentManager = new AtomLeafAgentManager<>(a -> space.makeVector() , box);
        pcForce.setAgentManager(atomAgentManager);
        IteratorDirective id = new IteratorDirective();
        id.includeLrc = false;
        potentialMaster.calculate(box, id, pcForce);
        tmpDrr1 = space.makeTensor();
    }

    public Tensor[][][][] calculateSum(AtomicTensorAtomicPair atomicTensorAtomicPair) {
    	Vector[] kv = kFactory.getWaveVectors();
//    	System.out.println(Arrays.toString(kv));
    	IMoleculeList molList = box.getMoleculeList();
    	Vector posl0 = space.makeVector();
    	Vector poslp = space.makeVector();
    	Vector dRpR0 = space.makeVector();

    	Tensor[][][] sumR = new Tensor[basisDim][basisDim][kv.length];
    	Tensor[][][] sumI = new Tensor[basisDim][basisDim][kv.length];
        
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                for(int k=0; k<kv.length; k++) { //kv.length = 1
                	sumR[jp][j][k] = new TensorND(6);
                	sumI[jp][j][k] = new TensorND(6);
                }
            }
        }
        
    	posl0.E(molList.get(0).getChildList().get(0).getPosition());

        for(int j=0; j<basisDim; j++) {//1st u.c.
        	System.out.println("j = "+j+" out of "+basisDim );
        	IMolecule moleculej = molList.get(j);
        	int L = box.getMoleculeList().size()/basisDim;
            for(int lp=0; lp<L; lp++) {//basisDim=46
            	if(lp==0){
            		poslp.E(molList.get(0).getChildList().get(0).getPosition());
            		dRpR0.Ev1Mv2(poslp, posl0);
            	}
                for(int jp=0; jp<basisDim; jp++) {//basisDim=46
                	if(jp==j && lp==0) continue; // skip tp next jp iter.: self j=jp
                	
                	IMolecule moleculejp = molList.get(basisDim*lp+jp);
                	if(jp==0){
                		poslp.E(moleculejp.getChildList().get(0).getPosition());
                		dRpR0.Ev1Mv2(poslp, posl0);
                	}
                	Tensor D6jjp = atomicToMolecularD(atomicTensorAtomicPair,moleculej,moleculejp);
                	for(int k=0; k<kv.length; k++){
                		sumR[j][jp][k].PEa1Tt1(Math.cos(kv[k].dot(dRpR0)),D6jjp);
                		sumI[j][jp][k].PEa1Tt1(Math.sin(kv[k].dot(dRpR0)),D6jjp);
                	}
                }//jp
            }//lp
            
    		poslp.E(molList.get(0).getChildList().get(0).getPosition());//acc. ineffic.
    		dRpR0.Ev1Mv2(poslp, posl0);
        	Tensor D6jj = atomicToMolecularD(atomicTensorAtomicPair,moleculej,moleculej);
        	for(int k=0; k<kv.length; k++){
        		sumR[j][j][k].PEa1Tt1(Math.cos(kv[k].dot(dRpR0)),D6jj);
        		sumI[j][j][k].PEa1Tt1(Math.sin(kv[k].dot(dRpR0)),D6jj);
        	}
            
        }//j mol
        return new Tensor[][][][] {sumR, sumI};
    }
    
    
    
    public Tensor atomicToMolecularD(AtomicTensorAtomicPair aTensor, IMolecule mol0, IMolecule mol1){

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

//    	com0.E(mol0.getChildList().getAtom(2).getPosition()); // O (-5.970371160466783, 5.978273273935142, -2.996126942837739)
//    	com1.E(mol1.getChildList().getAtom(2).getPosition()); // O (-6.016203213551466, 6.025148464416224, -2.996521341193713)
//    	com0.PE(4.4);
//    	com1.PE(4.4);
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
        		
        		D3tr_.E(D3tt_);	
        		D3tr_.TE(Rk);  
        		D3tr.PE(D3tr_);
        		
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
    	return D6mol;
    }
    public WaveVectorFactorySimple getWaveVectorFactory(){
    	return kFactory;
    }
    
    private final int basisDim;
    protected final Box box;
    protected final Space space;
    protected IMoleculePositionDefinition atomPosDef;
    protected final Vector com0, com1;
    protected WaveVectorFactorySimple kFactory;
    protected final Tensor[] tmpAtomicTensor3;//46X4=184 atoms dimentional array of 3dim Tensor
	protected PotentialMaster potentialMaster;
	protected AtomLeafAgentManager<Vector> atomAgentManager;
	protected   Tensor tmpDrr1;

    public interface AtomicTensorAtomicPair{
    	public Tensor atomicTensor(IAtom atom0, IAtom atom1);
    }
}
