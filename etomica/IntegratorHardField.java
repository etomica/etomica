package simulate;

/** 
 * Integration of molecules interacting via hard potentials, with a constant
 * force acting on each atom
 */

public class IntegratorHardField extends IntegratorHard {

  protected transient final double[] dr = new double[Space.D];
  
  public IntegratorHardField() {
    super();
  }

// --------- advances all atom coordinates and velocities by tStep --------

  protected void advanceAcrossTimeStep(double tStep) {
    
    double t2 = 0.5*tStep*tStep;
    for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        a.decrementCollisionTime(tStep);
        if(a.isForceFree()) {
            Space.uEa1Tv1(dr,tStep*a.rm,a.p);
            a.translate(dr);         //needs modification for nonspherical atom
        }
        else {
            Space.uEa1Tv1Pa2Tv2(dr,tStep*a.rm,a.p,t2*a.rm,a.f);
            a.translate(dr);         //needs modification for nonspherical atom
            Space.uEa1Tv1(dr,tStep,a.f);
            a.accelerate(dr);
  //          System.out.println(a.rm+" "+a.r[1]+" "+a.p[1]+" "+a.f[1]+" "+a.diameter);
        }
    }
  }
  
    public final void scaleMomenta(double s) {
      double rs = 1.0/s;
      for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        Space.uTEa1(a.p,s);
      }
      for(Atom a=firstPhase.firstAtom; a!=null; a=a.getNextAtom()) {
        if(a.isForceFree() && a.getCollisionPartner().isForceFree()) {
            a.collisionTime *= rs;
        }
        else {  //could improve efficiency by checking to see if collision time had decreased
            upList(a);
        }
        if(!a.isForceFree()) {downList(a);}
      }
      findNextCollider();
    }
  
}

