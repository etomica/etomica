package simulate;

import java.beans.*;
import java.awt.*;

public class PotentialHardDiskWall extends simulate.Potential
{
    double collisionDiameter, sig2;

    public PotentialHardDiskWall(double d) {
        setCollisionDiameter(d);
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
        
        double time = Double.MAX_VALUE;
        double dr, t;
        double dtdr;
        int i;
        
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}
        
        dr = wall.r[i] - disk.r[i];
        dtdr = 1.0/(disk.p[i]*disk.rm);
        t = dr*dtdr;

        if(t > 0.0) {time = Math.max(0.0,t-collisionDiameter/2.0*Math.abs(dtdr));}
        return time;
    }
    
    public void bump(Atom atom1, Atom atom2)
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
        
        if(wall.isVertical()) {
            disk.p[0] *= -1;
        }
        if(wall.isHorizontal()) {
            disk.p[1] *= -1;
        }
        
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        sig2 = c*c;
    }
}