/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * EAM (Embedded Atom Potential) potential with a lattice sum
 * 
 * @author Sabry Moustafa
 */
public class PotentialEAM_LS extends PotentialN implements PotentialSoft{

    protected double n, m, eps, a, c, rC1, rC2;
    protected Boundary boundary;
    protected final Vector dr, dR;
    protected Vector[] gradient;
    protected Vector[] rhograd;
    protected double [][] secondder;
    protected final Vector Lxyz, drtmp, dRtmp, R1,R2,R3;
    protected final int[] nShells;
    protected final Vector[] a0;
    protected final CoordinateDefinition coordinateDefinition;

    
    public PotentialEAM_LS(CoordinateDefinition coordinateDefinition, Space space, double n, double m, double eps, double a, double c, double rC, Vector[] a0) {
        super(space);
        
        this.n = n;
        this.m = m;
        this.eps = eps;
        this.a = a;
        this.c = c;
        this.coordinateDefinition = coordinateDefinition;
        rC1 = rC;
        rC2 = rC1;
        dr=space.makeVector();
        dR=space.makeVector();
        gradient=new Vector[0];
        rhograd=new Vector[0];
		dRtmp = space.makeVector();
		drtmp = space.makeVector();
        R1 = space.makeVector();
        R2 = space.makeVector();
        R3 = space.makeVector();

    	Lxyz = space.makeVector();
        this.a0 = a0;
		nShells = new int[] {(int) Math.ceil(rC/a0[0].getX(0) - 0.49999), (int) Math.ceil(rC/a0[1].getX(1) - 0.49999), (int) Math.ceil(rC/a0[2].getX(2) - 0.49999)};
		System.out.println("nShells = " + nShells[0] + " "+ nShells[1] + " "+ nShells[2] );
    }
    
    public double getRange() {
        return rC1;
    }
    public void setRange(double rC) {
        rC1 = rC;
        rC2 = rC1;
    }
    //coordina
    public double energy(IAtomList atoms) {
      double sumV=0;
      double rhoi=0;
      double rij, Rij, Lij;
      Vector ipos=atoms.getAtom(0).getPosition();
      Vector Ri = coordinateDefinition.getLatticePosition(atoms.getAtom(0));
      Vector shiftR = space.makeVector();

      
      for(int j=1;j<atoms.getAtomCount();j++){
        Vector jpos=atoms.getAtom(j).getPosition();
        Vector Rj = coordinateDefinition.getLatticePosition(atoms.getAtom(j));
        dr.Ev1Mv2(ipos, jpos);
        dR.Ev1Mv2(Ri, Rj);
        shiftR.E(dR);
        boundary.nearestImage(dR);
        shiftR.ME(dR);
        dr.ME(shiftR);
        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
          R1.setX(0, nx*a0[0].getX(0)); R1.setX(1, nx*a0[0].getX(1)); R1.setX(2, nx*a0[0].getX(2));
          Lxyz.E(R1);
	      for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
	        R2.setX(0, ny*a0[1].getX(0)); R2.setX(1, ny*a0[1].getX(1)); R2.setX(2, ny*a0[1].getX(2));
	        Lxyz.Ev1Pv2(R1, R2);
	        for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
	        R3.setX(0, nz*a0[2].getX(0)); R3.setX(1, nz*a0[2].getX(1)); R3.setX(2, nz*a0[2].getX(2));
	        Lxyz.Ev1Pv2(R1, R2); Lxyz.PE(R3);

  			  drtmp.Ev1Pv2(dr, Lxyz);
  			  dRtmp.Ev1Pv2(dR, Lxyz);
			  Lij = Math.sqrt(Lxyz.squared());
			  rij = Math.sqrt(drtmp.squared());
			  Rij = Math.sqrt(dRtmp.squared());
			  //pair pot.
		      if(j==1 && Lij<=rC1 && Lij > 0){ //self with 1/2 (Not needed for n-body!)
			    sumV += 0.5* eps *Math.pow(a /Lij , n);
		      }
		      if(Rij<=rC1 && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
			    sumV += eps *Math.pow(a /rij , n);
		      }
		      //n-body pot.
		      if(j==1 && Lij<=rC2  && Lij > 0){ //self
		        rhoi += Math.pow(a /Lij , m);
		      }
		      if(Rij<=rC2){
		        rhoi += Math.pow(a /rij , m);
		      }
	        }
	      }
	    }
      }
      double frho = -eps * c *Math.sqrt(rhoi);
      return sumV + frho;
    }
    
    public void setBox(Box box) {
        boundary=box.getBoundary();
    }

    public double virial(IAtomList atoms) {
      Vector ipos=atoms.getAtom(0).getPosition();
      Vector Ri = coordinateDefinition.getLatticePosition(atoms.getAtom(0));
      Vector gij2b = space.makeVector();
      Vector gijnb = space.makeVector();
      Vector shiftR = space.makeVector();

      double rhoi = 0;
      double vir2b = 0;
      double virnb = 0;
      double dvdr, Lij;
      double rij, Rij;
      double drhodr;

      for(int j=1;j<atoms.getAtomCount();j++){
        Vector jpos =atoms.getAtom(j).getPosition();
        Vector Rj = coordinateDefinition.getLatticePosition(atoms.getAtom(j));
        dr.Ev1Mv2(ipos, jpos);
        dR.Ev1Mv2(Ri, Rj);
        shiftR.E(dR);
        boundary.nearestImage(dR);
        shiftR.ME(dR);
        dr.ME(shiftR);

        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
            R1.setX(0, nx*a0[0].getX(0)); R1.setX(1, nx*a0[0].getX(1)); R1.setX(2, nx*a0[0].getX(2));
            Lxyz.E(R1);
      	  for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
  	        R2.setX(0, ny*a0[1].getX(0)); R2.setX(1, ny*a0[1].getX(1)); R2.setX(2, ny*a0[1].getX(2));
  	        Lxyz.Ev1Pv2(R1, R2);
      	    for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
    	        R3.setX(0, nz*a0[2].getX(0)); R3.setX(1, nz*a0[2].getX(1)); R3.setX(2, nz*a0[2].getX(2));
    	        Lxyz.Ev1Pv2(R1, R2); Lxyz.PE(R3);
        	  drtmp.Ev1Pv2(dr, Lxyz);
        	  dRtmp.Ev1Pv2(dR, Lxyz);
			  Lij = Math.sqrt(Lxyz.squared());
        	  rij = Math.sqrt(drtmp.squared());
        	  Rij = Math.sqrt(dRtmp.squared());
        	  if(j==1 && Lij<=rC1 && Lij > 0){
          	    dvdr =  -eps * n / a *Math.pow(a /Lij , n +1.0) ;
    			gij2b.Ea1Tv1(dvdr/Lij, Lxyz);
    			vir2b += 0.5*gij2b.dot(Lxyz); // Note the 1/2 here
          	  }
        	  if(Rij<=rC1 && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
        	    dvdr =  -eps * n / a *Math.pow(a /rij , n +1.0) ;
  			    gij2b.Ea1Tv1(dvdr/rij, drtmp);
  			    vir2b += gij2b.dot(drtmp);
        	  }

		      if(j==1 && Lij<=rC2  && Lij > 0){ //self
	            rhoi += Math.pow(a /Lij , m);
          	    drhodr= -m / a *Math.pow(a /Lij, m +1);
	        	gijnb.Ea1Tv1(drhodr/Lij, Lxyz);
	  			virnb += gijnb.dot(Lxyz);//WHY no 1/2? Bcs Fij== fij-fji = 2fij and 1/2*Fij=fij (which we compute)
			  }
			  if(Rij<=rC2){
			    rhoi += Math.pow(a /rij , m);
          	    drhodr= -m / a *Math.pow(a /rij, m +1);
	        	gijnb.Ea1Tv1(drhodr/rij, drtmp);
	  			virnb += gijnb.dot(drtmp);
			  }
  		    }
      	  }
   	    }
      }                  	  
      double f=Math.sqrt(rhoi);
      virnb *= -eps * c /2.0/f;
      double virial = vir2b + virnb;
      return virial;
    }

    
    
    public Vector[] gradient(IAtomList atoms) {
      if(gradient.length<atoms.getAtomCount()){
        rhograd=new Vector[atoms.getAtomCount()];
        gradient=new Vector[atoms.getAtomCount()];
        for(int j=0;j<atoms.getAtomCount();j++){
          gradient[j]=space.makeVector();
          rhograd[j]=space.makeVector();
        }
      }
      gradient[0].E(0);
      Vector ipos=atoms.getAtom(0).getPosition();
      Vector Ri = coordinateDefinition.getLatticePosition(atoms.getAtom(0));
        Vector shiftR = space.makeVector();
        double rhoi=0;
    	double dvdr;
        double rij, Rij, Lij;
        double drhodr;
        //S: Do NOT start from 0 as it will be -SUM(j>0); see below
        for(int j=1;j<atoms.getAtomCount();j++){
          gradient[j].E(0);
          rhograd[j].E(0);
          Vector jpos=atoms.getAtom(j).getPosition();
          Vector Rj = coordinateDefinition.getLatticePosition(atoms.getAtom(j));
          dr.Ev1Mv2(ipos, jpos);
          dR.Ev1Mv2(Ri, Rj);
          shiftR.E(dR);
          boundary.nearestImage(dR);
          shiftR.ME(dR);
          dr.ME(shiftR);
    	  for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
              R1.setX(0, nx*a0[0].getX(0)); R1.setX(1, nx*a0[0].getX(1)); R1.setX(2, nx*a0[0].getX(2));
              Lxyz.E(R1);
    	    for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
    	        R2.setX(0, ny*a0[1].getX(0)); R2.setX(1, ny*a0[1].getX(1)); R2.setX(2, ny*a0[1].getX(2));
    	        Lxyz.Ev1Pv2(R1, R2);
      	   	  for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
      	        R3.setX(0, nz*a0[2].getX(0)); R3.setX(1, nz*a0[2].getX(1)); R3.setX(2, nz*a0[2].getX(2));
    	        Lxyz.Ev1Pv2(R1, R2); Lxyz.PE(R3);
      			drtmp.Ev1Pv2(dr, Lxyz);
      			dRtmp.Ev1Pv2(dR, Lxyz);
  			    Lij = Math.sqrt(Lxyz.squared());
    			rij = Math.sqrt(drtmp.squared());
    			Rij = Math.sqrt(dRtmp.squared());
    			
    			//2-body
	            if(Rij<=rC1 && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
			      dvdr =  -eps * n / a *Math.pow(a /rij , n +1.0) ;
			      gradient[j].PEa1Tv1(-dvdr/rij, drtmp);
		        }
	            //n-body
		        if(j==1 && Lij<=rC2  && Lij > 0){ //self
			      rhoi += Math.pow(a /Lij , m);
			    }
		        if(Rij<=rC2){
		          rhoi += Math.pow(a /rij , m);
		          drhodr= -m / a *Math.pow(a /rij, m +1);
		          rhograd[j].PEa1Tv1(-drhodr/rij, drtmp);
		        }
    	      }
    	    }
    	  }          
        }//End j
        double f=Math.sqrt(rhoi);
        for (int j=1;j<atoms.getAtomCount();j++){
            gradient[j].PEa1Tv1(-eps * c /2.0/f,rhograd[j]);//Adds the n-body to the 2-body for j
            gradient[0].ME(gradient[j]);
        }
        return gradient;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
 }
