/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * EFS (Extended Finnis-Sinclair) potential
 * 
 * @author Joe Kromph
 */
public class PotentialEFS extends PotentialN implements PotentialSoft{

    protected double rCc, c0, c1, c2, c3, c4;
    protected double A, rCd, B;
    protected Vector dr, drij, drik, drjk;
    protected Boundary boundary;
    protected Vector[] gradient;
    protected Vector[] rhograd;
    protected double [][] secondder;
    
    public PotentialEFS(Space space, double A, double B, double c, double d, double c0, double c1,
                        double c2, double c3, double c4) {
        super(space);
        
        this.rCc=c;
        this.rCd=d;
        this.c0=c0;
        this.c1=c1;
        this.c2=c2;
        this.c3=c3;
        this.c4=c4;
        this.A=A;
        this.B=B;
        dr=space.makeVector();
        drij=space.makeVector();
        drik=space.makeVector();
        drjk=space.makeVector();
        gradient=new Vector[0];
        rhograd=new Vector[0];
    }
    
    public double getRange() {
        return rCd; // d>c
    }
    public void setRange(double rC) {
    	rCd = rC; //d>c
    }

    
    public double energy(IAtomList atoms) {
        double sumV=0;
        double rhoi=0;
        Vector ipos=atoms.getAtom(0).getPosition();
        
        for(int j=1;j<atoms.getAtomCount();j++){
          Vector jpos=atoms.getAtom(j).getPosition();
          dr.Ev1Mv2(ipos, jpos);
          boundary.nearestImage(dr);
          double rij=Math.sqrt(dr.squared());
          double rij_c = rij - rCc;
          if(rij<=rCc && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
            sumV+=rij_c*rij_c*(c0+c1*rij+c2*rij*rij+c3*rij*rij*rij+c4*rij*rij*rij*rij);
       	  }
    	  if(rij<=rCd){
   		    double rij_d=rij-rCd;
		    rhoi += rij_d*rij_d*(1+B*B*rij_d*rij_d);
   		  }
        }//j
        
        double frho=A*Math.sqrt(rhoi);
        double Utot=sumV-frho;
        return Utot;
        
    }
    
    public void setBox(Box box) {
        boundary=box.getBoundary();
    }

    public double virial(IAtomList atoms) {
        double virial=0;
        gradient(atoms);
        Vector ipos=atoms.getAtom(0).getPosition();
        for(int j=1;j<atoms.getAtomCount();j++){
            Vector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            virial -= gradient[j].dot(dr);
        }
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
        
        double rhoi=0;
        //S: Does NOT start from 0 as it will be -SUM(j>0); see below
        for(int j=1;j<atoms.getAtomCount();j++){
            gradient[j].E(0);
            rhograd[j].E(0);
            Vector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            double rij_c = rij - rCc;
            if(rij<=rCc && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
            	double dvdr;
                dvdr=rij_c*(2*(c0+rij*(c1+rij*(c2+rij*(c3+c4*rij))))+(rij-rCc)*(c1+rij*(2*c2+rij*(3*c3+4*c4*rij))));
                gradient[j].Ea1Tv1(-dvdr/rij, dr);
            }
            if(rij<=rCd){
                double rij_d = rij-rCd;
                double drhodr;
        		rhoi += A*A*rij_d*rij_d*(1+B*B*rij_d*rij_d);
                drhodr= A*A*(2*rij_d+4*B*B*rij_d*rij_d*rij_d);
                rhograd[j].Ea1Tv1(-drhodr/rij, dr);
            }
        }
        double f=Math.sqrt(rhoi);
        for (int j=1;j<atoms.getAtomCount();j++){
            gradient[j].PEa1Tv1(-1.0/2.0/f,rhograd[j]);
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
        double rhoi=0;
    	double dvdr, d2vdr2, g1, g2;
        double rij, rij2 , rij3 , rij4, sij;
        double rik, sik;
        double dphidrij=0,d2phidrij2=0, dudrho, d2udrho2;
        double dphidrik;
        double tmpD, tmpD1;
		
		//Get rhoi
        for(int j=1;j<ng;j++){
            Vector jpos=atoms.getAtom(j).getPosition();
            drij.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(drij);
            rij=Math.sqrt(drij.squared());
            if(rij<=rCd){
       		    double rij_d=rij-rCd;
    		    rhoi += rij_d*rij_d*(1+B*B*rij_d*rij_d);
            }
        }
        rhoi = A*A*rhoi;
        dudrho   = -1.0/2.0/Math.pow(rhoi, 0.5);
        d2udrho2 =  1.0/4.0/Math.pow(rhoi, 1.5);
        
        for(int j=1;j<ng;j++){
          Vector jpos=atoms.getAtom(j).getPosition();
          drij.Ev1Mv2(ipos, jpos);
          boundary.nearestImage(drij);

          rij=Math.sqrt(drij.squared()); 
          rij2 = rij*rij;   rij3 = rij*rij2;   rij4 = rij2*rij2;
          
          /**pairwise*/
          if(rij<=rCc && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
            double rij_c = rij - rCc;
            dvdr=rij_c*(2*(c0+rij*(c1+rij*(c2+rij*(c3+c4*rij))))+(rij-rCc)*(c1+rij*(2*c2+rij*(3*c3+4*c4*rij))));
            g1 = rij*dvdr;
            d2vdr2 =  2.0*(c0+c1*rij+c2*rij2+c3*rij3+c4*rij4)
            		+ 4.0*rij_c*(c1+2*c2*rij+3*c3*rij2+4*c4*rij3)
            		+ rij_c*rij_c*(2*c2+6*c3*rij+12*c4*rij2);
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
          }
           /**EAM: 1st derivative terms*/
           if(rij > rCd) continue;
           double rij_d = rij-rCd;
           double rij_d2 = rij_d*rij_d; 
           double rij_d3 = rij_d*rij_d*rij_d; 
           dphidrij   = A*A*(2.0*rij_d + 4*B*B*rij_d3);
           d2phidrij2 = A*A*(2.0 + 12.0*B*B*rij_d2);
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

              if(rik<=rCd){
              	double rik_d = rik - rCd;
            	double rik_d3 = rik_d*rik_d*rik_d;
            	dphidrik   = A*A*(2.0*rik_d + 4*B*B*rik_d3);
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
                    /**No need to do ji as we have jk and kj*/
//                  secondder[3*k+c][3*j+c] =  tmpD; // kj Indirect
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
                      /**No need to do ji as we have jk and kj*/
//                    secondder[3*k+r][3*j+c] =  tmpD1; // kj Indirect: xy
//                    secondder[3*k+c][3*j+r] =  tmpD;  // kj Indirect: yx
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
