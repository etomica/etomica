package etomica;

/**
 * An unphysical potential that has the sole action of swapping the species identity of
 * any molecule of a certain species that passes a plane (a wall).  
 * Momentum of molecule is unchanged by interaction with wall.  
 * One of species1Index and species2Index should correspond to a wall type species, and the 
 * other a disk type.
 * Presently this potential is suitable only for monatomic molecules from speciesDisks
 */
public class PotentialHardDiskSpeciesSwitchWall extends PotentialHardDiskWall
{
  /**
   * Molecules crossing the wall are changed to a molecule from this species
   */
    SpeciesDisks changeSpecies;
    
    double collisionDiameter, collisionRadius, sig2;

    public PotentialHardDiskSpeciesSwitchWall() {
        setCollisionDiameter(0.0);
    }
 /**
  * Always returns zero
  */
  public double energy(AtomPair pair) {return 0.0;}
 /**
  * Always returns zero
  */
  public double energyLRC(int n1, int n2, double V) {return 0.0;}
  
 /**
  * Always returns false
  */
  public boolean overlap(AtomPair pair) {return false;}
 
  /**
   * Molecules interacting with wall are changed to a molecule from this species
   */
    public void setChangeSpecies(Species s) {
        changeSpecies = (SpeciesDisks)s;
    }
    
   /**
    * Implements "collision dynamics" between wall and molecule.  Has net effect
    * of changing replacing the entering molecule with a new molecule of type
    * given from changeSpecies
    */
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