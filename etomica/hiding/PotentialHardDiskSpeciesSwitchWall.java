package simulate;

import java.beans.*;
import java.awt.*;

public class PotentialHardDiskSpeciesSwitchWall extends Potential implements PotentialHard
{
    SpeciesDisks changeSpecies;
    
    double collisionDiameter, collisionRadius, sig2;

    public PotentialHardDiskSpeciesSwitchWall() {
        setCollisionDiameter(0.0);
    }

    public void setChangeSpecies(Species s) {
        changeSpecies = (SpeciesDisks)s;
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
        
        int i;
        
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}
        
        if(parentPhase.noGravity || i==0 || !(wall.isStationary() || disk.isStationary())) {
            double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
            double dv = wall.p[i]*wall.rm-disk.p[i]*disk.rm;
            double time = -dr/dv;
            return (time > 0.0) ? time : Double.MAX_VALUE;
        }
        else {  //This could be greatly simplified given that collisionRadius = 0
            double time = Double.MAX_VALUE;
            double dr, dv;
            dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
            dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
            if(Math.abs(dr) < collisionRadius) {   //this may still need some work
                return (dr*dv > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
            dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
            double a = wall.isStationary() ? -(parentPhase.getG()) : parentPhase.getG();
            double discrim = dv*dv - 2*a*dr;
            if(discrim > 0) {
                boolean adr = (a*dr > 0);
                boolean adv = (a*dv > 0);
                int aSign = (a > 0) ? +1 : -1;
                if(adr && adv) {time = Double.MAX_VALUE;}
                else if(adr) {time = (-dv - aSign*Math.sqrt(discrim))/a;}
                else if(-a*dr/(dv*dv) < 1.e-7) {if(dr*dv<0) time = -dr/dv*(1+0.5*dr*a/(dv*dv));} //series expansion for small acceleration
                else {time = (-dv + aSign*Math.sqrt(discrim))/a;}
            }
            return time;
        }
    }
    
// not suited for multiatomic molecules; need to work on IntegratorHard (advanceToCollision method) to make ready
    
    public void bump(AtomHard atom1, AtomHard atom2)  //this needs updating to check for isStationary
    {
        double eps = 1.0e-6;
        
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
        
       Molecule m = disk.parentMolecule;
       Species oldSpecies = m.parentSpecies;
       oldSpecies.deleteMolecule(m);
       changeSpecies.addMolecule(m);
       for(Atom a=m.firstAtom(); a!=m.lastAtom.getNextAtom(); a=a.getNextAtom()) {
          ((AtomDisk)a).setDiameter(changeSpecies.getDiameter());
       }
       
       //Ensure wall and atom are separating
        int i;
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}
        double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
        double dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
        if(dr*dv > 0.0) { //OK, they are separating
            return;}  
        else {            //otherwise, put atom just on other side of wall
            int sign = (dv > 0.0) ? -1 : +1;
            disk.r[i] = wall.r[i] + sign*eps;
        }
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
        sig2 = c*c;
    }
}