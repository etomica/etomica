/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.ElectronVolt;

public class PotentialCuLREP extends PotentialN implements PotentialSoft{

    protected double c0, c1, c2, c3, c4, mLREP, nLREP;
    protected double A, rC1, rC2, B, r0;
    protected Vector dr, drij, drik, drjk;
    protected IBoundary boundary; 
    protected Vector[] gradient;
    protected Vector[] rhograd;
    protected double [][] secondder;
    protected int test = 1;
    
    public PotentialCuLREP(Space space) {
        super(space);
        
                this.mLREP   = 4.0;
                this.nLREP   = 6.0;
                this.rC1 = 6.1;
                this.rC2 = 7.8;
                this.c0  = ElectronVolt.UNIT.toSim(0.123554);
                this.c1  = ElectronVolt.UNIT.toSim(-0.134361);
                this.c2  = ElectronVolt.UNIT.toSim(0.0543818);
                this.c3  = ElectronVolt.UNIT.toSim(-0.981194E-2);
                this.c4  = ElectronVolt.UNIT.toSim(0.675816E-3);
                this.A   = ElectronVolt.UNIT.toSim(Math.sqrt(0.656618E-4));
                this.A = A*A;
                this.B   = 1.836569;
                this.r0  = 2.552655;
                System.out.println(" LREP parameters: A: "+A+" B: "+B+" rC1: "+rC1+" rC2: "+rC2+" c0: "+c0+" c1: "+c1+" c2: "+c2+" c3: "+c3+" c4: "+c4 + " A: " + A + " B: " + B + " r0: " + r0);

        dr=space.makeVector();
        drij=space.makeVector();
        drik=space.makeVector();
        drjk=space.makeVector();
        gradient=new Vector[0];
        rhograd=new Vector[0];
    }
    
    public double getRange() {
        return rC2; // rC2=7.8>rC1=6.1
    }
    public void setRange(double rC1, double rC2) {
    	this.rC1 = rC1; //d>c
    	this.rC2 = rC2; //d>c
    }

    
    public double energy(IAtomList atoms) {
        double sumV=0;
        double rhoi=0;
        double rij_c, rij_d , v1, v2;
        int ng = atoms.getAtomCount();
        if(test==1){
        	System.out.println(" ng: " + ng);
        	test = 0;
        }
        Vector ipos=atoms.getAtom(0).getPosition();
        for(int j=1;j<ng;j++){
          Vector jpos=atoms.getAtom(j).getPosition();
          dr.Ev1Mv2(ipos, jpos);
          boundary.nearestImage(dr);
          double rij=Math.sqrt(dr.squared());
          if(rij<=rC1 && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
            rij_c = rij - rC1 ;
            v1  = Math.pow(rij_c, mLREP);
            v2  = c0 + rij*(c1 + rij*(c2 + rij*(c3 + c4*rij))) ;
            sumV += v1*v2 ;
       	  }
    	  if(rij<=rC2){
   		    rij_d=rij-rC2;
		    rhoi += Math.pow(rij_d, nLREP)*Math.exp(-B*(rij/r0-1.0));
   		  }
        }//j
        rhoi = A*rhoi;
        
        double Utot=sumV-Math.sqrt(rhoi);
        return Utot;//getLeafIndex
        
    }
    
    public void setBox(Box box) {
        boundary=box.getBoundary();
    }

    public double virial(IAtomList atoms) {
        double virial=0;
        int ng = atoms.getAtomCount();
        gradient(atoms);
        Vector ipos=atoms.getAtom(0).getPosition();
        for(int j=1;j<ng;j++){
            Vector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            virial -= gradient[j].dot(dr);
        }
        return virial;
    }

    public Vector[] gradient(IAtomList atoms) {
        double v1, v2, dv1, dv2, dvdr, drhodr, dudrho, rij, rij_c, rij_d;
        int ng = atoms.getAtomCount();
        if(gradient.length<ng){
            rhograd=new Vector[ng];
            gradient=new Vector[ng];
            
            for(int j=0;j<ng;j++){
                gradient[j]=space.makeVector();
                rhograd[j]=space.makeVector();
            }
        }
        
        Vector ipos= atoms.getAtom(0).getPosition();
        int iAtom   = atoms.getAtom(0).getLeafIndex();
        double rhoi = 0;
        gradient[0].E(0);
        //S: Does NOT start from 0 as it will be -SUM(j>0); see below
        for(int j=1;j<ng;j++){
            gradient[j].E(0);
            rhograd[j].E(0);
            Vector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            rij=Math.sqrt(dr.squared());
            //2-body
            if(rij<=rC1 && iAtom < atoms.getAtom(j).getLeafIndex()){
                rij_c = rij - rC1 ;
                v1  = Math.pow(rij_c, mLREP);
                dv1 =   mLREP*v1/rij_c;
                v2  = c0 + rij*(c1 + rij*(c2 + rij*(c3 + c4*rij))) ;
                dv2 = c1 + rij*(2*c2 + rij*(3*c3 + 4*c4*rij)) ;
                dvdr = dv1*v2 + v1*dv2 ; 
                gradient[j].Ea1Tv1(-dvdr/rij, dr);
            }
            //n-body
            if(rij<=rC2){
                rij_d  = rij - rC2;
                double expB = Math.exp(-B*(rij/r0-1.0));
                double rij_d_n = Math.pow(rij_d, nLREP);
		        rhoi  += rij_d_n*expB;
                drhodr = expB*rij_d_n/rij_d * (nLREP - B/r0*rij_d);
                rhograd[j].Ea1Tv1(drhodr/rij, dr);
            }
        }//end j loop
        //ui = -sqrt(rhoi)
        dudrho = -1.0/2.0/Math.sqrt(A*rhoi);
        for (int j=1; j<ng; j++){
          gradient[j].PEa1Tv1(-dudrho * A , rhograd[j]);//-ve because we need Fj, not Fi.
          gradient[0].ME(gradient[j]);
        }
        return gradient;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
    
    
    
    public double [][] secondder(IAtomList atoms){
		int ng = atoms.getAtomCount();
    	secondder = new double[3*ng][3*ng];
        Vector ipos=atoms.getAtom(0).getPosition();
        double v1, v2, dv1, dv2, d2v1, d2v2;
    	double dvdr, d2vdr2, g1, g2;
        double rij, rij2 , rij3 , rij4;
        double rik;
        double dphidrij=0,d2phidrij2=0, dudrho, d2udrho2;
        double dphidrik;
        double tmpD, tmpD1;
        double rij_d, rij_c;
		double phi1, dphi1, d2phi1;
		double phi2, dphi2, d2phi2;
		//rhoi
        double rhoi = 0;
        for(int j=1;j<ng;j++){
            Vector jpos=atoms.getAtom(j).getPosition();
            drij.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(drij);
            rij=Math.sqrt(drij.squared());
            if(rij<=rC2){
       		    rij_d=rij-rC2;
		        rhoi += Math.pow(rij_d, nLREP)*Math.exp(-B*(rij/r0-1.0));
            }
        }
        rhoi     = A*rhoi;
        dudrho   = -1.0/2.0/Math.pow(rhoi, 0.5);
        d2udrho2 =  1.0/4.0/Math.pow(rhoi, 1.5);
        
        int iAtom = atoms.getAtom(0).getLeafIndex();
        
        for(int j=1;j<ng;j++){
          Vector jpos=atoms.getAtom(j).getPosition();
          drij.Ev1Mv2(ipos, jpos);
          boundary.nearestImage(drij);

          rij=Math.sqrt(drij.squared()); 
          
          rij2 = rij*rij;   rij3 = rij*rij2;   rij4 = rij2*rij2;

          
          /**pairwise*/
          if(rij<=rC1 && iAtom < atoms.getAtom(j).getLeafIndex()){
            rij_c = rij - rC1;
            
            v1  = Math.pow(rij_c, mLREP);
            dv1 =   mLREP*v1/rij_c;
            d2v1 = (mLREP-1)*dv1/rij_c;
            
            v2   = c0 + c1*rij +   c2*rij2 +   c3*rij3 +    c4*rij4; 
            dv2  =      c1     + 2*c2*rij  + 3*c3*rij2 +  4*c4*rij3; 
            d2v2 =               2*c2      + 6*c3*rij  + 12*c4*rij2; 
            
            dvdr = dv1*v2 + v1*dv2 ; 
            g1 = rij*dvdr;
            
            d2vdr2 =  2*dv1*dv2 + v1*d2v2 + d2v1*v2;
            g2 = rij2*d2vdr2;
            
            for (int c=0; c<3; c++){
              tmpD = -g1/rij2 + (g1-g2)/rij4*drij.getX(c)*drij.getX(c);
              secondder[c][3*j+c] += tmpD;//ij: xx
              secondder[3*j+c][c] += tmpD;//ji: xx
              for (int r=c+1; r<3; r++){
            	tmpD = (g1-g2)/rij4*drij.getX(r)*drij.getX(c);
            	secondder[r][3*j+c] += tmpD;//ij: xy
            	secondder[c][3*j+r] += tmpD;//ij: yx
            	
            	secondder[3*j+r][c] += tmpD;//ji: xy
            	secondder[3*j+c][r] += tmpD;//ji: yx
              }
            }
          }//if rij<rC1
          
           /**EAM: 1st derivative terms*/
          
           if(rij > rC2) continue;
           
           rij_d   = rij-rC2;
           
           phi1   = Math.pow(rij_d, nLREP);
           dphi1  = nLREP*phi1/rij_d;
           d2phi1 = (nLREP-1)*dphi1/rij_d;
           
           phi2   = Math.exp(-B*(rij/r0-1.0));
           dphi2  = -B/r0*phi2; 
           d2phi2 = -B/r0*dphi2; 
           
	       dphidrij = phi1*dphi2 + dphi1*phi2;
	       dphidrij = A*dphidrij;
           
           d2phidrij2 = 2*dphi1*dphi2 + d2phi1*phi2 + phi1*d2phi2;
           d2phidrij2 = A*d2phidrij2;
           for (int c=0; c<3; c++){
             tmpD = -1/rij*dudrho*dphidrij + 1/rij3*dudrho*dphidrij*drij.getX(c)*drij.getX(c);
             secondder[c][3*j+c] += tmpD; // ij delta_ab: xx
             secondder[3*j+c][c] += tmpD; // ji delta_ab: xx
             for (int r=c+1; r<3; r++){
               tmpD = 1/rij3*dudrho*dphidrij*drij.getX(r)*drij.getX(c);
               secondder[r][3*j+c] += tmpD;//ij: xy
               secondder[c][3*j+r] += tmpD;//ij: yx
               
               secondder[3*j+r][c] += tmpD;//ji: xy
               secondder[3*j+c][r] += tmpD;//ji: yx
             }
           }
            
           
            
            /**EAM: 1st&2nd derivative terms*/
            //k-loop:             
            for(int k=1;k<ng;k++){
              Vector kpos=atoms.getAtom(k).getPosition();
              drik.Ev1Mv2(ipos, kpos);
              boundary.nearestImage(drik);
              rik=Math.sqrt(drik.squared());

              if(rik<=rC2){
              	double rik_d = rik - rC2;
            	
                phi1   = Math.pow(rik_d, nLREP);
                dphi1  = nLREP*phi1/rik_d;
                d2phi1 = (nLREP-1)*dphi1/rik_d;
                
                phi2   = Math.exp(-B*(rik/r0-1.0));
                dphi2  = -B/r0*phi2; 
                d2phi2 = -B/r0*dphi2; 
                
                dphidrik = phi1*dphi2 + dphi1*phi2;
     	        dphidrik = A*dphidrik;

            	
                for (int c=0; c<3; c++){
                //xx
                  //ij
                  tmpD = -d2udrho2*dphidrij*dphidrik*drij.getX(c)*drik.getX(c)/rij/rik;
                  secondder[c][3*j+c] += tmpD; //ij Direct: xx 
                  secondder[3*j+c][c] += tmpD; //ji Direct: xx
                  if(j==k){//k can equal j for Dij:Direct
                 	tmpD = -dudrho*d2phidrij2*drij.getX(c)*drij.getX(c)/rij/rij;
                    secondder[c][3*j+c] += tmpD;// ij Direct (j=k): xx
                    secondder[3*j+c][c] += tmpD;// ji Direct (j=k): xx
                  }
                  //jk
                  if(j != k){
                  	tmpD = d2udrho2*dphidrij*dphidrik*drij.getX(c)*drik.getX(c)/rij/rik;
                    secondder[3*j+c][3*k+c] =  tmpD; // jk Indirect
                  }
                  
                  
                  
                  
              //xy
                  for (int r=c+1; r<3; r++){
                  //ij
                  	tmpD  = -d2udrho2*dphidrij*dphidrik*drij.getX(c)*drik.getX(r)/rij/rik;
                	tmpD1 = -d2udrho2*dphidrij*dphidrik*drij.getX(r)*drik.getX(c)/rij/rik;
                   	secondder[r][3*j+c]  += tmpD;  // ij Direct: xy
                    secondder[c][3*j+r]  += tmpD1; // ij Direct: yx
                	secondder[3*j+r][c]  += tmpD1; // ji Direct: xy
                	secondder[3*j+c][r]  += tmpD;  // ji Direct: yx
                    if(j==k){//k can equal j for Dij:Direct
                    	tmpD = -dudrho*d2phidrij2*drij.getX(c)*drij.getX(r)/rij/rij;
                        secondder[r][3*j+c]   += tmpD;// ij direct: xy
                        secondder[c][3*j+r]   += tmpD;// ij direct: yx
                        secondder[3*j+r][c]   += tmpD;// ji direct: xy
                        secondder[3*j+c][r]   += tmpD;// ji direct: yx
                    }
                    
                  //jk
                    if(j != k){
                     tmpD  = d2udrho2*dphidrij*dphidrik*drij.getX(r)*drik.getX(c)/rij/rik;
                     tmpD1 = d2udrho2*dphidrij*dphidrik*drij.getX(c)*drik.getX(r)/rij/rik;
                      secondder[3*j+r][3*k+c] =  tmpD;  // jk Indirect: xy
                      secondder[3*j+c][3*k+r] =  tmpD1; // jk Indirect: yx
                    }
                  }
                }
              }//if rik<rc2
            }//k  
            
            
        }//j
        
        //Self
        for(int i=0;i<ng;i++){
          for(int j=0;j<ng;j++){
          	if(i != j){
              for (int c=0; c<3; c++){
                for (int r=0; r<3; r++){
                  secondder[3*i+r][3*i+c] -=  secondder[3*i+r][3*j+c];              		
                }
              }          		
          	}
          }//j
        }//i
                return secondder;
    }
}
