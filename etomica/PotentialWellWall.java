package simulate;

import java.beans.*;
import java.awt.*;
import java.beans.Beans;

// Works only for horizontal wall

public class PotentialWellWall extends simulate.Potential
{

  private double coreDiameter, coreDiameterSquared;
  private double wellDiameter, wellDiameterSquared;
  private double lambda; //wellDiameter = coreDiameter * lambda
  private double epsilon;
  private transient final double[] r12 = new double[Space.D];
  private transient final double[] v12 = new double[Space.D];

  double collisionDiameter, collisionRadius, sig2;

    public PotentialWellWall(double d, double lambda, double epsilon) 
    {
        setCoreDiameter(d);
        setLambda(lambda);
        setEpsilon(epsilon);
    }

    public double collisionTime(Atom atom1, Atom atom2) {
   
        Atom disk;
        AtomWall wall;
        if(atom2 instanceof AtomWall) {
           disk = atom1;
           wall = (AtomWall)atom2;
        }
        else {
           disk = atom2;
           wall = (AtomWall)atom1;
        }
       int i=1;
       
       double dr = wall.r[i] - disk.r[i];
       double dv = 0.0-disk.p[i]*disk.rm;
       
       double bij = dr*dv;    //r12 . v12

       double adr = Math.abs(dr);
    double adv = Math.abs(dv);
    if(adr < wellDiameter) {  // Already inside wells

      if(adr < coreDiameter) {   // Inside core; collision now if approaching, at well if separating
        return (bij < 0) ? 0.0 : (wellDiameter-adr)/adv;
      }

      if(bij < 0.0) {    // Check for hard-core collision
        return  (adr-coreDiameter)/adv;
      }
      else {           // Moving away from each other, wells collide next
	    return (wellDiameter-adr)/adv;
      }
      
    }
    else {              // Outside wells; look for collision at well
        return (bij < 0) ? (adr-wellDiameter)/adv : Double.MAX_VALUE;
    }
    }
    
    public void bump(Atom atom1, Atom atom2)  //this needs updating to check for isStationary
    {
        Atom disk;
        AtomWall wall;
        if(atom2 instanceof AtomWall) {
           disk = atom1;
           wall = (AtomWall)atom2;
        }
        else {
           disk = atom2;
           wall = (AtomWall)atom1;
        }
        
 /*     if(wall.isVertical()) {
           disk.p[0] *= -1;
      }
        if(wall.isHorizontal()) {
            disk.p[1] *= -1;
        }*/
       
    double eps = 1.0e-6;
    int i=1;
       double dr = wall.r[i] - disk.r[i];
       double dv = -disk.p[i]*disk.rm;
       
       double bij = dr*dv;                           //r12 . v12

       double adr = Math.abs(dr);
    double adv = Math.abs(dv);
//    space.uEr1Mr2(r12,atom2.r,atom1.r);  //use instance method   //r2-r1
//    Space.uEa1Tv1Ma2Tv2(v12,atom2.rm,atom2.p,atom1.rm,atom1.p);  //v2-v1
    double r2 = dr*dr;
    // ke is kinetic energy due to components of velocity
//    double reduced_m = 1.0/((1.0/atom1.rm+ 1.0/atom2.rm)*atom1.rm*atom2.rm);
    double reduced_m = disk.mass;
    double ke = bij*bij*reduced_m/(2.0*r2);
    double s, r2New;
       if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
         s = 2.0*reduced_m*bij/r2;
         r2New = r2;
         disk.p[1] *= -1;
         return;
       }
       else {    // Well collision
        if(bij > 0.0) {         // Separating
	      if(ke < epsilon) {     // Not enough kinetic energy to escape
	         s = 2.0*reduced_m*bij/r2;
	         r2New = (1-eps)*wellDiameterSquared;
	         disk.p[1] *= -1;
	         
	       }
	    else {                 // Escape
//	  s = (0.5*bij/r - Math.sqrt(0.5*(ke - epsilon)))/r;
	         s = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m))/r2;
//             disk.p[1] += s*r12[1];
             disk.p[1] += s*dr;
	         r2New = (1+eps)*wellDiameterSquared;
             }
        }
        else {                  // Approaching
//	s = (0.5*bij/r + Math.sqrt(0.5*(ke + epsilon)))/r;
	         s = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m))/r2;
//             disk.p[1] += s*r12[1];
             disk.p[1] += s*dr;
    	     r2New = (1-eps)*wellDiameterSquared; 
          }
        }    
        
        if(r2New == r2) return;
        if(dr > 0) {
            disk.r[1] = wall.r[1] - Math.sqrt(r2New);
        }
        else {
            disk.r[1] = wall.r[1] + Math.sqrt(r2New);
        }
     }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
        sig2 = c*c;
        setCoreDiameter(c);
    }
    
      public double getCoreDiameter() {return coreDiameter;}
    public void setCoreDiameter(double c) {
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }

    public double getLambda() {return lambda;}
    public void setLambda(double lam) {
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    
    public double getEpsilon() {return epsilon;}
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
}