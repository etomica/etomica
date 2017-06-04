/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.api.*;
import etomica.atom.MoleculePositionCOM;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IMoleculePositionDefinition;
import etomica.potential.IPotentialMolecularSecondDerivative;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Tensor3D;

/** 
 * 3-point potential for water that can calculate gradient and torque (for the
 * center of mass of the water molecule).
 */
public class P2Water4PSoft extends P2Water4P implements IPotentialMolecularSecondDerivative  {

    public P2Water4PSoft(final Space space, double sigma, double epsilon,
                         double chargeH, double rCut, IMoleculePositionDefinition positionDefinition) {
        super(space, sigma, epsilon, chargeH,rCut, positionDefinition);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
		torque = new Vector[2];
		torque[0] = space.makeVector();
		torque[1] = space.makeVector();
		tau = new Vector[2];
		tau[0] = space.makeVector();
		tau[1] = space.makeVector();
		secondDerivative = new Tensor[3];
		this.secondDerivative[0] = space.makeTensor();
		this.secondDerivative[1] = space.makeTensor();
		this.secondDerivative[2] = space.makeTensor();
		tmpSecondD = new Tensor[3];
		this.tmpSecondD[0] = space.makeTensor();
		this.tmpSecondD[1] = space.makeTensor();
		this.tmpSecondD[2] = space.makeTensor();
		
		
		tmpDrr = space.makeTensor();
		
		rotationTensor3D = new RotationTensor3D();//debug only
		
		tmpTensor = new Tensor[4];
		for(int i=0;i<4;i++){
			tmpTensor[i] = space.makeTensor();
		}
		
		fWork = space.makeVector();
		dr1 = space.makeVector();
		dr0 = space.makeVector();
		dr = space.makeVector();
		gradientAndTorque = new Vector[][]{gradient,torque};
		epsilon48 = epsilon*48.0;
		
		//debug only
		a = new Vector[2][4];
		for(int i=0;i<2;i++){
			for(int j=0;j<4;j++){
			a[i][j] = space.makeVector();
			}
		}
		centerMass = space.makeVector();//debug only
    }

    public Vector[][] gradientAndTorque(IMoleculeList pair){
        IMolecule water1 = pair.getMolecule(0);
        IMolecule water2 = pair.getMolecule(1);

        //compute O-O distance to consider truncation	
        Vector O1r = (water1.getChildList().getAtom(2)).getPosition();
        Vector O2r = (water2.getChildList().getAtom(2)).getPosition();

        work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
        boundary.nearestImage(work);
        shift.PE(work);
        double r2 = work.squared();

        if(r2 < 1.6) {
            throw new RuntimeException("waters are overlapped in gradient");
        }
        
        com1.E(positionDefinition.position(water1));
        com2.E(positionDefinition.position(water2));
        dr.Ev1Mv2(com1, com2);
        dr.PE(shift);
        
        if(dr.squared() > rCut*rCut){
        	gradient[0].E(0);
        	gradient[1].E(0);
        	torque[0].E(0);
        	torque[1].E(0);
        	return   	gradientAndTorque;
        }
        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        double du = -epsilon48*s6*(s6 - 0.5);// u = epsilon4*s6*(s6 - 1.0)
	
        gradient[0].Ea1Tv1(du/r2,work);//force vector 

        work.Ev1Mv2(O2r, com2);
        work.XE(gradient[0]);
        torque[1].E(work);
        work.Ev1Mv2(com1, O1r);
        work.XE(gradient[0]);
        torque[0].E(work);
        
		Vector H11r = water1.getChildList().getAtom(0).getPosition();
		Vector H12r = water1.getChildList().getAtom(1).getPosition();
		Vector H21r = water2.getChildList().getAtom(0).getPosition();
		Vector H22r = water2.getChildList().getAtom(1).getPosition();
        Vector M1r = water1.getChildList().getAtom(3).getPosition();
        Vector M2r = water2.getChildList().getAtom(3).getPosition();
        
        // M1-M2
        work.Ev1Mv2(M1r, M2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMM/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(M2r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, M1r);
        work.XE(fWork);
        torque[0].PE(work);

        // M1-H21
        work.Ev1Mv2(M1r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, M1r);
        work.XE(fWork);
        torque[0].PE(work);

        // M1-H22
        work.Ev1Mv2(M1r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, M1r);
        work.XE(fWork);
        torque[0].PE(work);

        // H11-M2
        work.Ev1Mv2(H11r, M2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(M2r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H11r);
        work.XE(fWork);
        torque[0].PE(work);

        // H11-H21
        work.Ev1Mv2(H11r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H11r);
        work.XE(fWork);
        torque[0].PE(work);

        // H11-H22
        work.Ev1Mv2(H11r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H11r);
        work.XE(fWork);
        torque[0].PE(work);

        // H12-M2
        work.Ev1Mv2(H12r, M2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(M2r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H12r);
        work.XE(fWork);
        torque[0].PE(work);

        // H12-H21
        work.Ev1Mv2(H12r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H12r);
        work.XE(fWork);
        torque[0].PE(work);

        // H12-H22
        work.Ev1Mv2(H12r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H12r);
        work.XE(fWork);
        torque[0].PE(work);

        gradient[1].Ea1Tv1(-1, gradient[0]);

		return gradientAndTorque;
	}
    
   public Tensor [] secondDerivative(IMoleculeList pair) {
	   secondDerivative[0].E(0);
	   secondDerivative[1].E(0);
	   secondDerivative[2].E(0);
	   IMolecule water1 = pair.getMolecule(0);
	   IMolecule water2 = pair.getMolecule(1);
	   Vector O1 = water1.getChildList().getAtom(2).getPosition();
	   Vector O2 = water2.getChildList().getAtom(2).getPosition();
	   
	   work.Ev1Mv2(O1, O2);
	   shift.Ea1Tv1(-1,work);
	   boundary.nearestImage(work);
	   shift.PE(work);
	   
	   if(work.squared() < 1.6) {
		   throw new RuntimeException("waters are overlapped in gradient");
	   }
	   
	   com1.E(positionDefinition.position(water1));
	   com2.E(positionDefinition.position(water2));
	   dr.Ev1Mv2(com1, com2);
	   dr.PE(shift);
	   
//        System.out.println("O1 = " + O1);  
//        System.out.println("com1 = " + com1);
//        System.exit(2);
	   
	   if(dr.squared() > rCut*rCut){ 
		   return secondDerivative;
	   }
	   
	   //TODO 
	   Tensor[] t = atomicToMolecularD(water1, water2);
//	   Tensor[] t = atomicToMolecularD(water2, water1);// debug only
	   
	   secondDerivative[0].E(t[0]);
	   secondDerivative[1].E(t[1]);
	   secondDerivative[2].E(t[2]);
	   
	   
	   //debug only
//	   t[1].transpose();
//	   t[2].transpose();
//	   System.out.println("dudidj = \n" + t[0]);
//	   System.out.println("dudidi transpose = \n" + t[1]);
//	   System.out.println("dudjdj transpose = \n" + t[2]);
	   
//	   dr.E(gradientAndTorque(pair)[1][0]);
//	   System.out.println("torque = " + dr);
	   
	   //debug only TODO
	   double prec = 1E-4;
	   if(false){//dtau finite test
	   Vector tau0 = space.makeVector();
	   Vector taui = space.makeVector();
	   for(int i=0;i<3;i++){
	   tau0.E(gradientAndTorque(pair)[1][0]);
	   doRotation(i,water1,1,prec);
	   taui.E(gradientAndTorque(pair)[1][0]);
	   taui.ME(tau0);
	   taui.TE(-1/prec);
	   doRotation(i,water1,-1,prec);
	   System.out.println("dtaudxi" +taui);
	   }
	   System.exit(2);
	   }
	   
	   if(false){//dtau finite test
	   Vector tau0 = space.makeVector();
	   Vector taui = space.makeVector();
	   for(int i=0;i<3;i++){
	   tau0.E(gradientAndTorque(pair)[1][1]);
	   doRotation(i,water2,1,prec);
	   taui.E(gradientAndTorque(pair)[1][1]);
	   taui.ME(tau0);
	   taui.TE(-1/prec);
	   doRotation(i,water2,-1,prec);
	   System.out.println("dtaudxi" +taui);
	   }
	   System.exit(2);
	   }
	   return secondDerivative;
    }
   
   protected Tensor[] atomicToMolecularD(IMolecule mol0, IMolecule mol1){
	   tmpSecondD[0].E(0);
	   tmpSecondD[1].E(0);
	   tmpSecondD[2].E(0);
	   Tensor3D D3rr  = new Tensor3D();
	   Tensor3D D3rri = new Tensor3D();
	   Tensor3D D3rrj = new Tensor3D();
	   Tensor3D D3rr_ = new Tensor3D(); 
	   Tensor3D Rk = new Tensor3D(); 
	   Tensor3D Rk_ = new Tensor3D();
	   Tensor3D Rkp = new Tensor3D(); 
	   Tensor3D Rkp_ = new Tensor3D();
	   Vector Xk = space.makeVector();
	   Vector Xkp = space.makeVector();
	   MoleculePositionCOM com_0 = new MoleculePositionCOM(space);
	   Vector comk = com_0.position(mol0);
	   MoleculePositionCOM com_1 = new MoleculePositionCOM(space);
	   Vector comkp = com_1.position(mol1);
	   int numSites0 = mol0.getChildList().getAtomCount();
	   int numSites1 = mol1.getChildList().getAtomCount();
	   Vector f0 = space.makeVector();
	   Vector f1 = space.makeVector();
//	   IVectorMutable torque0 = space.makeVector();
//	   IVectorMutable torque1 = space.makeVector();
	   Vector fkp = space.makeVector();
	   Vector fk = space.makeVector();
	   
	   D3rr.E(0);
	   D3rri.E(0);
	   D3rrj.E(0);
	   for(int i=0;i<4;i++){
		   tmpTensor[i].E(0);
	   }
	   
	   for (int atomk=0; atomk < numSites0; atomk++){ 
		   fk.E(0);//force on atomk
		   IAtom atom0 = mol0.getChildList().getAtom(atomk);
		   Vector posk = atom0.getPosition();
		   Xk.Ev1Mv2(posk, comk);
		   boundary.nearestImage(Xk);
		   Rk.setComponent(0,0,0.0); Rk.setComponent(1,1,0.0);	Rk.setComponent(2,2,0.0);
		   Rk.setComponent(0,1, Xk.getX(2));  Rk.setComponent(1,0,-Xk.getX(2));   
		   Rk.setComponent(0,2,-Xk.getX(1));  Rk.setComponent(2,0, Xk.getX(1));
		   Rk.setComponent(1,2, Xk.getX(0));  Rk.setComponent(2,1,-Xk.getX(0));
		   
   			for (int atomkp=0; atomkp < numSites1; atomkp++){  
   			IAtom atom1 = mol1.getChildList().getAtom(atomkp);
   			Vector poskp = atom1.getPosition();
   			Xkp.Ev1Mv2(poskp, comkp);
   			boundary.nearestImage(Xkp);

   			Rkp.setComponent(0,0,0.0); Rkp.setComponent(1,1,0.0);	Rkp.setComponent(2,2,0.0);
   			Rkp.setComponent(0,1, Xkp.getX(2));  Rkp.setComponent(1,0,-Xkp.getX(2));      
   			Rkp.setComponent(0,2,-Xkp.getX(1));  Rkp.setComponent(2,0, Xkp.getX(1));
   			Rkp.setComponent(1,2, Xkp.getX(0));  Rkp.setComponent(2,1,-Xkp.getX(0));

   			D3rr_.E(aTensor.atomicTensor(atom0 , atom1));
   			Rk_.E(Rk);  Rk_.TE(D3rr_);         
   			Rk_.TE(Rkp); D3rr.ME(Rk_); 
   			tmpTensor[atomk].ME(D3rr_);
   			
//  			fk.PE(pairForce(atom0, atom1));
//  			fkp.E(pairForce(atom1, atom0));
  			
  			fk.PE(pairForce(atom1, atom0));
  			fkp.E(pairForce(atom0, atom1));
  			
//  			fkp.Ea1Tv1(fkp.dot(Xkp)/Xkp.squared(), Xkp);//TODO
//  			f1.PE(fkp);//TODO
  			
   	   		{//jj begin  
   	   			D3rr_.TE(-1.0);  
   	   			Rkp_.E(Rkp);
   	   			Rkp_.TE(D3rr_); 
   	   			Rkp_.TE(Rkp);
   	   			D3rrj.ME(Rkp_);
   	   			
   	   			Xkp.TE(-1);//be cautious that Xkp change sign here
   	   			tmpDrr.Ev1v2(Xkp, fkp);
   	   			tmpDrr.setComponent(0, 0,  -fkp.getX(1) * Xkp.getX(1) - fkp.getX(2) * Xkp.getX(2)); 
   	   			tmpDrr.setComponent(1, 1,  -fkp.getX(0) * Xkp.getX(0) - fkp.getX(2) * Xkp.getX(2)); 
   	   			tmpDrr.setComponent(2, 2,  -fkp.getX(0) * Xkp.getX(0) - fkp.getX(1) * Xkp.getX(1)); 
   	   			D3rrj.PE(tmpDrr);
   	   		}//jj end
  			
   		}//atomkp
   			
   		{//ii
   			D3rr_.E(tmpTensor[atomk]);
   			Rk_.E(Rk);
   			Rk_.TE(D3rr_); 
   			Rk_.TE(Rk);
   			D3rri.ME(Rk_);
   			
   			Xk.TE(-1);
   			tmpDrr.Ev1v2(Xk, fk);
   			tmpDrr.setComponent(0, 0,  -fk.getX(1) * Xk.getX(1) - fk.getX(2) * Xk.getX(2)); 
   			tmpDrr.setComponent(1, 1,  -fk.getX(0) * Xk.getX(0) - fk.getX(2) * Xk.getX(2)); 
   			tmpDrr.setComponent(2, 2,  -fk.getX(0) * Xk.getX(0) - fk.getX(1) * Xk.getX(1)); 
   			D3rri.PE(tmpDrr);
   		}
	   }//atomk

	   
	   tmpSecondD[0].PE(D3rr);
	   tmpSecondD[1].PE(D3rri);
	   tmpSecondD[2].PE(D3rrj);
	   
	   D3rr.E(0);
	   D3rri.E(0);
	   D3rrj.E(0);
//	   System.out.println("t[0] = "  + tmpSecondD[0]);
//	   System.out.println("t[1] = "  + tmpSecondD[1]);
//	   System.out.println("t[2] = "  + tmpSecondD[2]);
//	   System.exit(0);
	   return tmpSecondD;
   }

   
  	AtomicTensorAtomicPair aTensor = new AtomicTensorAtomicPair() {
  		final Tensor identity = new Tensor3D(new double[][] {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}});
  		Tensor D3ES = space.makeTensor();
  		Tensor tmpD3 = space.makeTensor();

  		//dr == r0 - r1
  		public Tensor atomicTensor(IAtom atom0, IAtom atom1) {
  			D3ES.E(0);
  			IMolecule water1 = atom0.getParentGroup();
  			IMolecule water2 = atom1.getParentGroup();
  			Vector O1 = (water1.getChildList().getAtom(2)).getPosition();
  			Vector O2 = (water2.getChildList().getAtom(2)).getPosition();

  			work.Ev1Mv2(O1, O2);
  			shift.Ea1Tv1(-1,work);
  			boundary.nearestImage(work);
  			shift.PE(work);
  			
  			if(work.squared() < 1.6) {
  				throw new RuntimeException("waters are overlapped in gradient");
  			}
  			
  			com1.E(positionDefinition.position(water1));
  			com2.E(positionDefinition.position(water2));
  			dr.Ev1Mv2(com1, com2);
  			dr.PE(shift);

  			double r2 = dr.squared();
  			
  			if(r2 > rCut*rCut ) return D3ES;
  			if((atom0.getIndex()==2 && atom1.getIndex() != 2) || ( atom0.getIndex()!=2 && atom1.getIndex() == 2)) return D3ES;


  			dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
  			dr.PE(shift);
  			r2= dr.squared();
  			// debug only
//  			System.out.println("atomPairs = ( " + + atom0.getIndex() + ", " + + atom1.getIndex() + ")" ) ;
//  			System.out.println("r = " + Math.sqrt(r2));
//  			System.out.println("sigma = " + sigma + " sigma2 = " + sigma2);
//  			System.out.println("rCut = " + rCut);
  			
  			tmpD3.Ev1v2(dr,dr);//outer product
  			double dW = 0;
  			double d2W = 0;
              
      		if(atom0.getIndex() == 2 &&  atom1.getIndex() == 2){
      			double s2 = sigma2/r2;
      			double s6 = s2*s2*s2;
      			dW = 24.0*epsilon*s6*(1-2*s6);
      			d2W = 24.0*epsilon*s6*(-7 + 26*s6);
      		}else{
      			double r = Math.sqrt(r2);
      			double q1 = (atom0.getIndex() == 3)?chargeM:chargeH;
      			double q2 = (atom1.getIndex() == 3)?chargeM:chargeH;
      			dW = -q1*q2/r;
      			d2W = 2.0*q1*q2/r;
      		}
      		
      		tmpD3.TE(1.0/(r2*r2)*(dW - d2W));
      		tmpD3.PEa1Tt1(-dW/r2,identity);
            D3ES.PE(tmpD3);
  	        
  			return D3ES;
  		}
  	}; 
  	
  	
  	public Vector pairForce(IAtom atom0, IAtom atom1 ){
  		fWork.E(0);
  		IMolecule water1 = atom0.getParentGroup();
  		IMolecule water2 = atom1.getParentGroup();
  		Vector O1 = (water1.getChildList().getAtom(2)).getPosition();
  		Vector O2 = (water2.getChildList().getAtom(2)).getPosition();

  		work.Ev1Mv2(O1, O2);
  		shift.Ea1Tv1(-1,work);
  		boundary.nearestImage(work);
  		shift.PE(work);

  		if(work.squared() < 1.6) {
  			throw new RuntimeException("waters are overlapped in gradient");
  		}
  		
  		com1.E(positionDefinition.position(water1));
  		com2.E(positionDefinition.position(water2));
  		dr.Ev1Mv2(com1, com2);
  		dr.PE(shift);
  		
  		if(dr.squared() > rCut*rCut ) return fWork;
  		if((atom0.getIndex()==2 && atom1.getIndex() != 2) || ( atom0.getIndex()!=2 && atom1.getIndex() == 2)) return fWork;
  		
  		work.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
  		work.PE(shift);
  		double r2 = work.squared();
  		
  		if(atom0.getIndex() == 2 &&  atom1.getIndex() == 2){
  			
  			double s2 = sigma2/r2;
  			double s6 = s2*s2*s2;
  			fWork.Ea1Tv1(24.0*epsilon*s6*(1-2*s6)/r2, work);
  		}else{
  			double q1 = (atom0.getIndex() == 3)?chargeM:chargeH;
  			double q2 = (atom1.getIndex() == 3)?chargeM:chargeH;
  			fWork.Ea1Tv1(-q1*q2/(r2*Math.sqrt(r2)), work);
  		}
  		return fWork;
  				
  	}
    
    public Vector[] gradient(IMoleculeList atoms) {
        // do extra work to calculate torque
        gradientAndTorque(atoms);
        return gradient;
    }
    
    public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
        gradientAndTorque(atoms);
        //FIXME  
        //pressureTensor.PEv1v2(gradient[0],dr);
        return gradient;
    }

    public double virial(IMoleculeList atoms) {
        //FIXME
        return 0;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
	public double getSigma() {return sigma;}
    
	public double getEpsilon() {return epsilon;}
	
	public double getChargeM(){return chargeM;} 
	public double getChargeH(){return chargeH;}
	
	
	
	//TODO debug only
    protected void doTransform(IMolecule molecule) {
    	IAtomList childList =  molecule.getChildList();
            for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
                IAtom atom = childList.getAtom(iChild);
                Vector r = atom.getPosition();
                r.ME(centerMass);
                boundary.nearestImage(r);
                rotationTensor3D.transform(r);
                r.PE(centerMass);
            }
    }
    
    protected void doRotation(int axis,IMolecule mol,double pOrm,double prec){
    	centerMass.E(positionDefinition.position(mol));
    	dr.E(0);
    	dr.setX(axis, pOrm);
    	rotationTensor3D.setRotationAxis(dr,prec);
    	doTransform(mol);
    }
	//debug only  TODO
	
    private static final long serialVersionUID = 1L;
	protected final Vector[] gradient, torque,tau;
	protected final Vector[][] gradientAndTorque;
	protected final Tensor[] secondDerivative, tmpSecondD ;	
	protected double epsilon48;
	protected final Vector fWork,dr,dr0,dr1;
	protected  final Tensor tmpDrr;
	protected final Tensor[] tmpTensor;
	
	
	
	
	protected transient RotationTensor3D rotationTensor3D;//debug only
	protected Vector centerMass;//debug only
	protected final Vector[][] a;//debug only
	
//	protected AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent> atomAgentManager;
	
	public interface AtomicTensorAtomicPair{
    	public Tensor atomicTensor(IAtom atom0, IAtom atom1);
    }
	
}
