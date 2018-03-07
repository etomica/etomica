/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * EAM (Embedded Atom Method) potential
 * 
 * @author Joe Kromph
 * @author Sabry Moustafa
 */
public class PotentialEAM extends PotentialN implements PotentialSoft{

    protected double n , m , eps , a , C , rC1, rC2;
    protected Boundary boundary;
    protected final Vector dr, drij, drik, drjk;
    protected Vector[] gradient;
    protected Vector[] rhograd;
    protected double [][] secondder;
    
    public PotentialEAM(Space space, double n , double  m , double  eps , double  a , double  C, double rC) {
        super(space);
        
        this.n=n;
        this.m=m;
        this.eps = eps;
        this.a = a;
        this.C = C;
        rC1 = rC;
        rC2 = rC1;
        dr=space.makeVector();
        drij=space.makeVector();
        drik=space.makeVector();
        drjk=space.makeVector();
        gradient=new Vector[0];
        rhograd=new Vector[0];
    }
    
    public double getRange() {
        return rC1;
    }
    public void setRange(double rC) {
        rC1 = rC;
        rC2 = rC1;
    }
    
    public double energy(IAtomList atoms) {
        double sumV=0;
        double rhoi=0;
        Vector ipos=atoms.get(0).getPosition();
        
        for(int j = 1; j<atoms.size(); j++){
            Vector jpos=atoms.get(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            if(rij<=rC1 && atoms.get(0).getLeafIndex() < atoms.get(j).getLeafIndex()){
              sumV += eps*Math.pow(a/rij , n);            		
            }
        	if(rij<=rC2){
       			rhoi += Math.pow(a/rij , m);
            }
        }
        double frho = -eps*C*Math.sqrt(rhoi);
        
        
        return sumV + frho;
        
        
    }
    
    public void setBox(Box box) {
        boundary=box.getBoundary();
    }

    public double virial(IAtomList atoms) {
        double virial=0;
        gradient(atoms);
        Vector ipos=atoms.get(0).getPosition();
        for(int j = 1; j<atoms.size(); j++){
            Vector jpos=atoms.get(j).getPosition();
            dr.Ev1Mv2(jpos,ipos);
            boundary.nearestImage(dr);
            virial += gradient[j].dot(dr);
        }
        return virial;
    }

    public Vector[] gradient(IAtomList atoms) {
        
        if(gradient.length<atoms.size()){
            rhograd=new Vector[atoms.size()];
            gradient=new Vector[atoms.size()];
            
            for(int j = 0; j<atoms.size(); j++){
                rhograd[j] = space.makeVector();
                gradient[j]= space.makeVector();
            }
        }
        
        gradient[0].E(0);
        Vector ipos=atoms.get(0).getPosition();
        
        double rhoi=0;
    	double dvdr;

        //S: Does NOT start from 0 as it will be -SUM(j>0); see below
        for(int j = 1; j<atoms.size(); j++){
            gradient[j].E(0);
            rhograd[j].E(0);
            Vector jpos=atoms.get(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            if(rij<=rC1 && atoms.get(0).getLeafIndex() < atoms.get(j).getLeafIndex()){
	            dvdr =  -eps*n/a*Math.pow(a/rij , n+1) ;
	            gradient[j].Ea1Tv1(-dvdr/rij, dr);//fji <= -fij
            }
            if(rij<=rC2){
                double drhodr;
                rhoi += Math.pow(a/rij , m);
                drhodr = -m/a*Math.pow(a/rij, m+1);
                rhograd[j].Ea1Tv1(-drhodr/rij, dr);//fji <= -fij
            }          
        }
        double f=Math.sqrt(rhoi);
        for (int j = 1; j<atoms.size(); j++){
            gradient[j].PEa1Tv1(-eps*C/2.0/f,rhograd[j]);//Adds the n-body to the 2-body for j
            gradient[0].ME(gradient[j]);
        }
        return gradient;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
    
    
    
    
    public double [][] secondder(IAtomList atoms){
		int ng = atoms.size();
    	secondder = new double[3*ng][3*ng];
        Vector ipos=atoms.get(0).getPosition();
        double rhoi=0;
    	double dvdr, d2vdr2, g1, g2;
        double rij, rij2 , rij3 , rij4, sij;
        double rik, sik;
        double dphidrij=0,d2phidrij2=0, dudrho, d2udrho2;
        double dphidrik;
        double tmpD, tmpD1;
		
		//Get rhoi
        for(int j=1;j<ng;j++){
            Vector jpos=atoms.get(j).getPosition();
            drij.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(drij);
            rij=Math.sqrt(drij.squared());
            sij=a/rij;
            if(rij<=rC2){
                rhoi += Math.pow(sij , m);
            }
        }
        dudrho   = -eps*C/2.0/Math.pow(rhoi, 0.5);
        d2udrho2 =  eps*C/4.0/Math.pow(rhoi, 1.5);
        
        
        for(int j=1;j<ng;j++){
          Vector jpos=atoms.get(j).getPosition();
          drij.Ev1Mv2(ipos, jpos);
          boundary.nearestImage(drij);

          rij=Math.sqrt(drij.squared()); 
          sij=a/rij;
          rij2 = rij*rij;   rij3 = rij*rij2;   rij4 = rij2*rij2;
          
          /**pairwise*/
          if(rij<=rC1 && atoms.get(0).getLeafIndex() < atoms.get(j).getLeafIndex()){
            dvdr =  -eps*n/a*Math.pow(sij , n+1) ;
            g1 = rij*dvdr;
            d2vdr2 =  eps*n*(n+1.0)/a/a*Math.pow(sij , n+2) ;
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
           if(rij>rC2) continue;
            
           dphidrij   =        -m/a*Math.pow(sij, m+1);
           d2phidrij2 = m*(m+1)/a/a*Math.pow(sij, m+2);
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
              Vector kpos=atoms.get(k).getPosition();
              drik.Ev1Mv2(ipos, kpos);
              boundary.nearestImage(drik);
              rik=Math.sqrt(drik.squared());
              sik=a/rik;

              if(rik<=rC2){
                dphidrik   = -m/a*Math.pow(sik, m+1);
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
        
        
/**DEBUG: F.D.*/   
//      System.out.println();
//      System.out.println();
//      System.out.println();
//      System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++");
//
//      int i = 11, j = 0;
//      for(int r=0;r<3;r++){
//        for(int c=0;c<3;c++){
//          System.out.print(secondder[3*i+r][3*j+c] + "  ");
//        }
//        System.out.println();
//      }
//
//      double U0 = energy(atoms);
//      IVectorMutable dx = space.makeVector();
//      IVectorMutable dy = space.makeVector();
//      double d = 0.001;
//      dx.setX(0, d);        dx.setX(1, 0);        dx.setX(2, 0);
//      dy.setX(0, 0);        dy.setX(1, d);        dy.setX(2, 0);
//      atoms.getAtom(i).getPosition().PE(dx);
//      atoms.getAtom(j).getPosition().PE(dy);
//      double Uipjp = energy(atoms);
//      atoms.getAtom(i).getPosition().ME(dx);
//      atoms.getAtom(i).getPosition().ME(dx);
//      double Uimjp = energy(atoms);
//      atoms.getAtom(j).getPosition().ME(dy);
//      atoms.getAtom(j).getPosition().ME(dy);
//      double Uimjm = energy(atoms);
//      atoms.getAtom(i).getPosition().PE(dx);
//      atoms.getAtom(i).getPosition().PE(dx);
//      double Uipjm = energy(atoms);
//      atoms.getAtom(i).getPosition().ME(dx);
//      atoms.getAtom(j).getPosition().PE(dy);
//      double Dij = 1/4.0/d/d*(Uipjp-Uipjm-Uimjp+Uimjm);
//      System.out.println(Dij);
//      System.out.println(secondder[3*i+0][3*j+1] - Dij);
//      System.out.println((secondder[3*i+0][3*j+1] - Dij)/Dij);
//      System.exit(0);

    return secondder;
    }
}

