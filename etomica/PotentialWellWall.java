package simulate;

import java.beans.*;
import java.awt.*;
import java.beans.Beans;

// Square-well potential between a disk and a wall
// Not written to handle moving wall
// No momentum or energy accumulation in wall is recorded

public class PotentialWellWall extends Potential implements PotentialHard
{

  private double coreDiameter, coreDiameterSquared;
  private double wellDiameter, wellDiameterSquared;
  private double lambda; //wellDiameter = coreDiameter * lambda
  private double epsilon; //well depth
  private double elasticity;  //see set and get methods for description

    public PotentialWellWall(double d, double lambda, double epsilon, double elasticity) 
    {
        setCoreDiameter(d);
        setLambda(lambda);
        setEpsilon(epsilon);
        setElasticity(elasticity);
    }

    public double collisionTime(AtomHard atom1, AtomHard atom2) {   
        AtomHardDisk disk;
        AtomHardWall wall;
        if(atom2 instanceof AtomHardWall) {
           disk = (AtomHardDisk)atom1;
           wall = (AtomHardWall)atom2;
        }
        else {
           disk = (AtomHardDisk)atom2;
           wall = (AtomHardWall)atom1;
        }
        
        int i = (wall.isHorizontal()) ? 1 : 0;

        double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
        double dv = 0.0-disk.p[i]*disk.rm;   //stationary wall
       
        double bij = dr*dv;

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
    
    public void bump(AtomHard atom1, AtomHard atom2)  //this needs updating to check for isStationary
    {
        AtomHardDisk disk;
        AtomHardWall wall;
        if(atom2 instanceof AtomHardWall) {
           disk = (AtomHardDisk)atom1;
           wall = (AtomHardWall)atom2;
        }
        else {
           disk = (AtomHardDisk)atom2;
           wall = (AtomHardWall)atom1;
        }
               
        int i = (wall.isHorizontal()) ? 1 : 0;
        double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
        double dv = -disk.p[i]*disk.rm;
       
        double bij = dr*dv;     //r12 . v12

        double adr = Math.abs(dr);
        double adv = Math.abs(dv);
        double r2 = dr*dr;
        
    // ke is kinetic energy due to components of velocity
//    double reduced_m = 1.0/((1.0/atom1.rm+ 1.0/atom2.rm)*atom1.rm*atom2.rm);
        double reduced_m = disk.mass;
        double ke = bij*bij*reduced_m/(2.0*r2);
        double s;
        int rSign;
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
            s = 2.0*reduced_m*bij/r2;
            disk.p[i] *= -1;
            return;
        }
        else {    // Well collision
            if(bij > 0.0) {         // Separating
	            if(ke < epsilon) {     // Not enough kinetic energy to escape
               //	                s = 2.0*reduced_m*bij/r2;
	                disk.p[i] *= -1;
	                rSign = -1;
	            }
	            else {                 // Escape
               //	  s = (0.5*bij/r - Math.sqrt(0.5*(ke - epsilon)))/r;
	                s = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m))/r2;
                    disk.p[i] /= elasticity;
                    disk.p[i] += s*dr;       //should check more carefully if this line comes before or after the previous one.
	                rSign = +1;
                }
            }
            else {                  // Approaching
               //	s = (0.5*bij/r + Math.sqrt(0.5*(ke + epsilon)))/r;
	            s = reduced_m*(bij + Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m))/r2;
                disk.p[i] += s*dr;
                disk.p[i] *= elasticity;
    	        rSign = -1; 
            }
            double eps = 1.0e-6;
            double offset = (1.0+rSign*eps)*wellDiameter;
            if(dr > 0) {
                disk.r[i] = wall.r[i] - offset;
            }
            else {
                disk.r[i] = wall.r[i] + offset;
            }
        }    
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

// To aid in getting molecule to stick on collision with wall, after passing through well the
// disk's energy is reduced by the factor "elasticity".  The energy is returned if the particle
// later manages to escape.  A perfectly elastic collision is obtained with elasticity = 1.0.
    public double getElasticity() {return elasticity;}
    public void setElasticity(double e) {
        elasticity = e;
    }
}