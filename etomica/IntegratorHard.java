package simulate;
import java.util.Observable;

public class IntegratorHard extends Integrator {

  private Atom nextCollider;
  protected transient final double[] dr = new double[Space.D];
    
  public IntegratorHard() {
    super();
  }

  
//--------------------------------------------------------------
// steps all particles across time interval tStep, handling
// all intervening collisions

  public void doStep(double tStep) {

    if(tStep < nextCollider.getCollisionTime()) {
      advanceAcrossTimeStep(tStep);
      if(isothermal) {
        scaleMomenta(Math.sqrt(this.temperature/firstPhase.getKineticTemperature()));
      }
      return;
    }

    double tStepNew = tStep - nextCollider.getCollisionTime();
    advanceToCollision();
    doStep(tStepNew);
    return;
  }

//--------------------------------------------------------------
 
  protected void findNextCollider() {
    //find next collision pair by looking for minimum collisionTime
    double minCollisionTime = Double.MAX_VALUE;
    for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
      if(a.getCollisionTime() < minCollisionTime) {
        minCollisionTime = a.getCollisionTime();
        nextCollider = a;}
    }
}
//--------------------------------------------------------------
// advances to next collision, applies collision dynamics to involved 
// species and updates collision time/partners

  protected void advanceToCollision() {
    
    advanceAcrossTimeStep(nextCollider.getCollisionTime());
    Atom partner = nextCollider.getCollisionPartner();
    nextCollider.getCollisionPotential().bump(nextCollider,partner);
    for(Atom a=firstPhase.firstAtom; a!=partner; a=a.getNextAtom()) {
        if(a.getCollisionPartner()==nextCollider || a.getCollisionPartner()==partner || a==nextCollider) {
            upList(a);
        }
    }
    upList(partner);
    downList(nextCollider);
    downList(partner);

    findNextCollider();
    
    firstPhase.updatedKineticEnergy = false;
    firstPhase.updatedPotentialEnergy = false;
}


// --------- advances all atom coordinates by tStep --------

  protected void advanceAcrossTimeStep(double tStep) {
    
    if(firstPhase.noGravity) {
        for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
            a.decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            Space.uEa1Tv1(dr,tStep*a.rm,a.p);
            a.translate(dr);         //needs modification for nonspherical atom
        }
    }
    else {
        double t2 = tStep*tStep;
        for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
            a.decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            Space.uEa1Tv1Pa2Tv2(dr,tStep*a.rm,a.p,t2,firstPhase.gravity.gVector);
            a.translate(dr);         //needs modification for nonspherical atom
            Space.uEa1Tv1(dr,tStep*a.mass,firstPhase.gravity.gVector);
            a.accelerate(dr);
        }
    }
  }

 /**
  * Update of collision list when gravity is changed.
  */
  public void update(Observable o, Object arg) {
    if(o instanceof Gravity) {
      for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        if(a.isStationary() || a.getCollisionPartner().isStationary()) {
            upList(a);
        }
        if(a.isStationary()) {downList(a);}
      }
      findNextCollider();
    }
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
//    Potential2[] p2 = firstPhase.potential2[atomSpeciesIndex];
    for(Atom a=nextMoleculeAtom; a!=null; a=a.getNextAtom()) {
//        Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
        Potential potential = firstPhase.potential2[a.getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < minCollisionTime) {
            minCollisionTime = time;
            atom.setCollision(time,a,potential);
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
        Potential potential = p1.getPotential(a,atom);
        double time = potential.collisionTime(a,atom);
        if(time < a.getCollisionTime()) {
            a.setCollision(time,atom,potential);
        }
    }
    
    //Loop through remaining downlist atoms in firstPhase
 //   Potential2[] p2 = firstPhase.potential2[atomSpeciesIndex];
    for(Atom a=previousMoleculeAtom; a!=null; a=a.getPreviousAtom()) {
        Potential potential = firstPhase.potential2[a.getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,a);
 //       Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < a.getCollisionTime()) {
            a.setCollision(time,atom,potential);
        }
    }
  }
  
//--------------------------------------------------------------

  public void initialize() {
    for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        upList(a);
    }
    findNextCollider();
  }
  
//--------------------------------------------------------------

    public void scaleMomenta(double s) {
      double rs = 1.0/s;
      for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        Space.uTEa1(a.p,s);
        a.collisionTime *= rs;
      }
    }
}

