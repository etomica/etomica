package simulate;
import java.util.Observable;

public class IntegratorHard extends Integrator {

  private Agent nextCollider;
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
    for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=null; a=a.getNextAtomC()) {
      Agent ia = (Agent)a.ia;
      double ct = ia.getCollisionTime();
      if( ct < minCollisionTime) {
        minCollisionTime = ct;
        nextCollider = ia;
      }
    }
}
//--------------------------------------------------------------
// advances to next collision, applies collision dynamics to involved 
// species and updates collision time/partners

  protected void advanceToCollision() {
    
    advanceAcrossTimeStep(nextCollider.getCollisionTime());
    AtomC partner = nextCollider.getCollisionPartner();
    if(partner == null) {
        upList(nextCollider.atom);
        downList(nextCollider.atom);
    }
    else {
//        Atom partnerNextAtom = partner.nextMoleculeFirstAtom();  //put this back in for multiatomic speciesSwitch; also need to do more work with loop below
        AtomC partnerNextAtom = partner.getNextAtomC();
        nextCollider.getCollisionPotential().bump(nextCollider.atom,partner);
        
        boolean upListedN = false;
        boolean upListedP = false;
        for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=partnerNextAtom; a=a.getNextAtomC()) {  //note that nextCollider's or partner's position in linked-list may have been moved by the bump method
            AtomC aPartner = ((Agent)a.ia).getCollisionPartner();
            if(aPartner==nextCollider.atom || aPartner==partner) {
                upList(a);
                if(a == nextCollider.atom) {upListedN = true;}
                else if(a == partner) {upListedP = true;}
            }
        }
        if(!upListedN) {upList(nextCollider.atom);}
        if(!upListedP) {upList(partner);}
        downList(nextCollider.atom);
        downList(partner);
    }

    findNextCollider();
}


// --------- advances all atom coordinates by tStep --------

  protected void advanceAcrossTimeStep(double tStep) {
    
    if(firstPhase.noGravity) {
        for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=null; a=a.getNextAtomC()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            Space.uEa1Tv1(dr,tStep*a.rm,a.p);
            a.translate(dr);         //needs modification for nonspherical atom
        }
    }
    else {
        double t2 = tStep*tStep;
        for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=null; a=a.getNextAtomC()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
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
      for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=null; a=a.getNextAtomC()) {
        if(a.isStationary() || ((Agent)a.ia).getCollisionPartner().isStationary()) {
            upList(a);
        }
        if(a.isStationary()) {downList(a);}
      }
      findNextCollider();
    }
  }
  
//--------------------------------------------------------------

  protected void upList(AtomC atom) {
    
    double minCollisionTime = Double.MAX_VALUE;
    Agent aia = (Agent)atom.ia;
    if(!atom.isStationary() && atom instanceof AtomDisk) {  //if mobile, set collision time to time atom takes to move half a box edge
        for(int i=Space.D-1; i>=0; i--) {
            double tnew = Math.abs((atom.parentMolecule.parentSpecies.parentPhase.space.getDimensions(i)-1.0001*((AtomDisk)atom).getDiameter())/atom.p[i]);  //assumes range of potential is .le. diameter
            minCollisionTime = (tnew < minCollisionTime) ? tnew : minCollisionTime;
        }
        minCollisionTime *= 0.5*atom.mass;
    }
    aia.setCollision(minCollisionTime, null, null);
    
    Atom nextMoleculeAtom = atom.nextMoleculeFirstAtom();  //first atom on next molecule
    int atomSpeciesIndex = atom.getSpeciesIndex();
    
    //Loop through remaining uplist atoms in this atom's molecule
    Potential1 p1 = firstPhase.potential1[atomSpeciesIndex];
    for(AtomC a=(AtomC)atom.getNextAtom(); a!=nextMoleculeAtom; a=a.getNextAtomC()) {
        PotentialHard potential = (PotentialHard)p1.getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < minCollisionTime) {
            minCollisionTime = time;
            aia.setCollision(time,a,potential);
        }
    }
    
    //Loop through remaining uplist atoms in firstPhase
//    Potential2[] p2 = firstPhase.potential2[atomSpeciesIndex];
    for(AtomC a=(AtomC)nextMoleculeAtom; a!=null; a=a.getNextAtomC()) {
//        Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
        PotentialHard potential = (PotentialHard)firstPhase.potential2[a.getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < minCollisionTime) {
            minCollisionTime = time;
            aia.setCollision(time,a,potential);
        }
    }
  }

//--------------------------------------------------------------

  protected void downList(AtomC atom) {
    
    AtomC previousMoleculeAtom = (AtomC)atom.getMolecule().firstAtom().getPreviousAtom();
    
    int atomSpeciesIndex = atom.getSpeciesIndex();
    
    //Loop through remaining downlist atoms in this atom's molecule
    Potential1 p1 = firstPhase.potential1[atomSpeciesIndex];
    for(AtomC a=(AtomC)atom.getPreviousAtom(); a!=previousMoleculeAtom; a=a.getPreviousAtomC()) {
        Agent aia = (Agent)a.ia;
        PotentialHard potential = (PotentialHard)p1.getPotential(a,atom);
        double time = potential.collisionTime(a,atom);
        if(time < aia.getCollisionTime()) {
            aia.setCollision(time,atom,potential);
        }
    }
    
    //Loop through remaining downlist atoms in firstPhase
 //   Potential2[] p2 = firstPhase.potential2[atomSpeciesIndex];
    for(AtomC a=previousMoleculeAtom; a!=null; a=a.getPreviousAtomC()) {
        Agent aia = (Agent)a.ia;
        PotentialHard potential = (PotentialHard)firstPhase.potential2[a.getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,a);
 //       Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < aia.getCollisionTime()) {
            aia.setCollision(time,atom,potential);
        }
    }
  }
  
//--------------------------------------------------------------

  public void initialize() {
    deployAgents();
    for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=null; a=a.getNextAtomC()) {
        upList(a);
    }
    findNextCollider();
  }
  
//--------------------------------------------------------------

    public void scaleMomenta(double s) {
      double rs = 1.0/s;
      for(AtomC a=(AtomC)firstPhase.firstAtom(); a!=null; a=a.getNextAtomC()) {
        Space.uTEa1(a.p,s);
        ((Agent)a.ia).collisionTime *= rs;
      }
    }
    
    public IntegratorAgent makeAgent(Atom a) {
        return new Agent((AtomC)a);
    }
    
    private class Agent implements IntegratorAgent {
      public AtomC atom;
      public double time0;  //time since last collision
      private double collisionTime = Double.MAX_VALUE; //time to next collision
      private AtomC collisionPartner;  //next atom scheduled for collision by atom containing this Agent
      private PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent

      public Agent(AtomC a) {atom = a;}
      
      public final void setCollision(double time, AtomC partner, PotentialHard p) {
        collisionTime = time;
        collisionPartner = partner;
        collisionPotential = p;
      }

      public final void decrementCollisionTime(double interval) {
        collisionTime -= interval;
        time0 += interval;
      }
      
      public final double getCollisionTime() {return collisionTime;}
      public final AtomC getCollisionPartner() {return collisionPartner;}
      public final PotentialHard getCollisionPotential() {return collisionPotential;}
  
    }
}

