package simulate;

import java.beans.*;
import java.awt.*;

public class PotentialHardDiskSpeciesSwitchWall extends PotentialHardDiskWall implements PotentialHard
{
    SpeciesDisks changeSpecies;
    
    double collisionDiameter, collisionRadius, sig2;

    public PotentialHardDiskSpeciesSwitchWall() {
        setCollisionDiameter(0.0);
    }

    public void setChangeSpecies(Species s) {
        changeSpecies = (SpeciesDisks)s;
    }
    
    
// not suited for multiatomic molecules; need to work on IntegratorHard (advanceToCollision method) to make ready
    
    public void bump(AtomPair pair) {
   
        double eps = 1.e-6;
        
        Atom disk;
        Atom wall;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
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
}