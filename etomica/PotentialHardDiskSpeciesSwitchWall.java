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
        boolean wallFirst;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
           wallFirst = false;
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
           wallFirst = true;
        }
        
       Molecule m = disk.parentMolecule;
       Space.Coordinate coord = m.coordinate();
       Phase phase = m.parentPhase();
       phase.deleteMolecule(m);
       m = changeSpecies.getMolecule();
       phase.addMolecule(m);
       m.translateTo(coord.position());
       m.accelerateTo(coord.momentum());
       disk = m.firstAtom();
       if(wallFirst) pair.reset(wall,disk);  //not suited to multiatomics
       else pair.reset(disk,wall);

       //Ensure wall and atom are separating
        int i = (((AtomType.Wall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate
//        double dr = wall.position(i) - disk.position(i);   //no PBC
//        double dv = wall.momentum(i)*wall.rm()-disk.momentum(i)*disk.rm();
        double dr = pair.dr(i);
        double dv = pair.dv(i);
        if(dr*dv > 0.0) { //OK, they are separating
            return;}  
        else {            //otherwise, put atom just on other side of wall
            int sign = (dv > 0.0) ? -1 : +1;
//            disk.r[i] = wall.r[i] + sign*eps;
            disk.position().PE(i,sign*eps);
        }
    }
}