package simulate;
import java.util.Enumeration;

public final class IntegratorHardNbrlist extends IntegratorHard {

  private Atom nextCollider;
  private transient final double[] dr = new double[Space.D];
    
  public IntegratorHardNbrlist() {
    super();
  }

  
// --------- advances all atom coordinates by tStep --------

  protected void advanceAcrossTimeStep(double tStep) {
    
    for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        a.decrementCollisionTime(tStep);
        Space.uEa1Tv1(dr,tStep*a.rm,a.p);
        a.translate(dr);         //needs modification for nonspherical atom
    }
    
    for(Molecule m=firstPhase.firstMolecule; m!=null; m=m.getNextMolecule()) {
        if(m.needNeighborUpdate()) {
          m.updateNeighborList();
          downList(m.firstAtom);   //change this for multiatomic molecule
          upList(m.firstAtom);
        }
    }
    findNextCollider();
  }
  
//--------------------------------------------------------------

  protected void upList(Atom atom) {
    
    Atom nextMoleculeAtom = atom.getMolecule().lastAtom.getNextAtom();  //first atom on next molecule
    double minCollisionTime = Double.MAX_VALUE;
    
    int atomSpeciesIndex = atom.getSpeciesIndex();
    atom.setCollision(Double.MAX_VALUE, null, null);
    
    //Loop through remaining uplist atoms in this atom's molecule
    Potential1 p1 = firstPhase.potential1[atomSpeciesIndex];
    for(Atom a=atom.getNextAtom(); a!=nextMoleculeAtom; a=a.getNextAtom()) {
        Potential potential = p1.getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < minCollisionTime) {
            minCollisionTime = time;
            atom.setCollision(time,a,potential);
        }
    }
    
    //Loop through remaining uplist atoms in firstPhase
    Potential2[] p2 = firstPhase.potential2[atomSpeciesIndex];
    Molecule molecule = atom.getMolecule();
//    if(molecule.needNeighborUpdate) {molecule.updateNeighborList();}
    for(Enumeration enum=molecule.getNeighborList(); enum.hasMoreElements();) {
       Molecule m = (Molecule)enum.nextElement();
       for(Atom a=m.firstAtom; a!=m.lastAtom.getNextAtom(); a=a.getNextAtom()) {
          Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
          double time = potential.collisionTime(atom,a);
          if(time < minCollisionTime) {
             minCollisionTime = time;
             atom.setCollision(time,a,potential);
          } 
        }
    }
  }

//--------------------------------------------------------------

  protected void downList(Atom atom) {
    
    Atom previousMoleculeAtom = atom.getMolecule().firstAtom.getPreviousAtom();
    
    int atomSpeciesIndex = atom.getSpeciesIndex();
    
    //Loop through remaining downlist atoms in this atom's molecule
    Potential1 p1 = firstPhase.potential1[atomSpeciesIndex];
    for(Atom a=atom.getPreviousAtom(); a!=previousMoleculeAtom; a=a.getPreviousAtom()) {
        Potential potential = p1.getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < a.getCollisionTime()) {
            a.setCollision(time,atom,potential);
        }
    }
    
    //Loop through remaining downlist atoms in firstPhase
    Potential2[] p2 = firstPhase.potential2[atomSpeciesIndex];
    Molecule molecule = atom.getMolecule();
    for(Molecule m=molecule.getPreviousMolecule(); m!=null; m=m.getPreviousMolecule()) {
//       if(m.needNeighborUpdate) {m.updateNeighborList();}
       if(m.hasNeighbor(molecule)) {
         for(Atom a=m.firstAtom; a!=m.lastAtom.getNextAtom(); a=a.getNextAtom()) {
           Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
           double time = potential.collisionTime(atom,a);
           if(time < a.getCollisionTime()) {a.setCollision(time,atom,potential);}
         }
       }
    }
  }
  
//--------------------------------------------------------------

  public void initialize() {
    
    for(Molecule m=firstPhase.firstMolecule; m!=null; m=m.getNextMolecule()) {
        m.updateNeighborList();
    }
    for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        upList(a);
    }
    findNextCollider();
  }
        
}

