package simulate;
import java.util.Observable;
import java.awt.Color;

public final class IntegratorHard extends Integrator {

private Agent nextCollider;
private AtomPair.Iterator.A upPairIterator, downPairIterator;
private Atom.Iterator upAtomIterator, downAtomIterator;
private AtomPair atomPair;
private PotentialHard spacePotential;

//boolean bb = true;  used in debugging
            
public IntegratorHard() {
    super();
}

public void registerPhase(Phase p) {
    super.registerPhase(p);
    upPairIterator = p.iterator.makeAtomPairIteratorUp();
    downPairIterator = p.iterator.makeAtomPairIteratorDown();
    upAtomIterator = p.iterator.makeAtomIteratorUp();
    downAtomIterator = p.iterator.makeAtomIteratorDown();
    atomPair = new AtomPair(p);
    spacePotential = (PotentialHard)p.space().makePotential();
}

          
//--------------------------------------------------------------
// steps all particles across time interval tStep, handling
// all intervening collisions

public void doStep(double tStep) {

    if(tStep < nextCollider.getCollisionTime()) {
        advanceAcrossTimeStep(tStep);
        if(isothermal) {
            scaleMomenta(Math.sqrt(this.temperature/firstPhase.kineticTemperature()));
        }
        debugMethod();
        return;
    }

    double tStepNew = tStep - nextCollider.getCollisionTime();
    advanceToCollision();
    doStep(tStepNew);
    return;
}

//debugging tools
// Colors atoms according to whether they are in the up or down neighborlist of the first atom
    private void debugMethod() {
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            a.setColor(Color.black);
        }
        Atom central = firstPhase.firstAtom();
        central.setColor(((Space2DCell.Coordinate)central.coordinate).cell.color);
        upPairIterator.reset(central,Iterator.INTRA);
        while(upPairIterator.hasNext()) {
            AtomPair pair = upPairIterator.next();
            pair.atom2.setColor(Color.blue);
        }
        downPairIterator.reset(central,Iterator.INTRA);
        while(downPairIterator.hasNext()) {
            AtomPair pair = downPairIterator.next();
            pair.atom2.setColor(Color.green);
        }
        /*
        
        upAtomIterator.reset(central);
        while(upAtomIterator.hasNext()) {
            Atom atom = upAtomIterator.next();
            atom.setColor(Color.blue);
        }
        downAtomIterator.reset(central);
        while(downAtomIterator.hasNext()) {
            Atom atom = downAtomIterator.next();
            atom.setColor(Color.green);
        }        
        */
        central.setColor(((Space2DCell.Coordinate)central.coordinate).cell.color);
        central.setColor(Color.red);
    }

//--------------------------------------------------------------
         
protected void findNextCollider() {
    //find next collision pair by looking for minimum collisionTime
    double minCollisionTime = Double.MAX_VALUE;
    for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
        Agent ia = (Agent)a.ia;
        double ct = ia.getCollisionTime();
        if( ct < minCollisionTime) {
            minCollisionTime = ct;
            nextCollider = ia;
        }
    }
    
    //Debugging code.  Colors some atoms according to whether they are next collider
//    if(minCollisionTime <= 0.0) {
//        AtomPair pair = firstPhase.makeAtomPair(nextCollider.atom,nextCollider.collisionPartner);
//        System.out.println("collision time "+minCollisionTime);
//    }
//    if(nextCollider.getCollisionPartner()!=null) {nextCollider.atom.setColor(Color.green);}
//    else if(bb) {nextCollider.atom.setColor(Color.blue); bb=!bb;}
//    else {nextCollider.atom.setColor(Color.red); bb=!bb;}
}
//--------------------------------------------------------------
// advances to next collision, applies collision dynamics to involved 
// species and updates collision time/partners

protected void advanceToCollision() {
            
    advanceAcrossTimeStep(nextCollider.getCollisionTime());
    Atom partner = nextCollider.getCollisionPartner();
    if(partner == nextCollider.atom) {  //"self-collision" is code for "collision" with space
        Atom a = nextCollider.atom;
        atomPair.reset(a,a);
        nextCollider.getCollisionPotential().bump(atomPair);
        upList(a);
        downList(a);
        
        //reset collision partners of atoms that are now up from this atom but still list it as their
        //collision partner.  Assumes this atom was moved down list, but this won't always be the case
        //This bit could be made more efficient
        upPairIterator.reset(a,Iterator.INTRA);
        while(upPairIterator.hasNext()) {
            AtomPair pair = upPairIterator.next();
            if(((Agent)pair.atom2().ia).getCollisionPartner() == a) {  //upList atom could have atom as collision partner if atom was just moved down list
                upList(pair.atom2());
            }
        }
        
        //to keep collision lists perfect, should do an upList on atoms that had this
        //atom on its neighbor list, but no longer do because it has moved away
    }
    else {
//        Atom partnerNextAtom = partner.nextMoleculeFirstAtom();  //put this back in for multiatomic speciesSwitch; also need to do more work with loop below
//        Atom partnerNextAtom = partner.coordinate.nextNeighbor().atom();
        Atom partnerNextAtom = null;  //remove this -- temporary
        atomPair.reset(nextCollider.atom,partner);
        nextCollider.getCollisionPotential().bump(atomPair);
//        nextCollider.getCollisionPotential().bump(nextCollider.atom,partner);
                
        boolean upListedN = false;
        boolean upListedP = false;
//        for(Atom a=firstPhase.firstAtom(); a!=partnerNextAtom; a=a.nextAtom()) {  //note that nextCollider's or partner's position in linked-list may have been moved by the bump method
//        upAtomIterator.reset(firstPhase.firstAtom());
        upAtomIterator.reset();  //first atom in first cell
        while(upAtomIterator.hasNext()) {
            Atom a = upAtomIterator.next();
            if(a == partnerNextAtom) break;
            Atom aPartner = ((Agent)a.ia).getCollisionPartner();
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
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            a.translateBy(tStep*a.rm(),a.momentum());
        }
    }
    else {
        double t2 = tStep*tStep;
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            a.translateBy(tStep*a.rm(),a.momentum());
            a.translateBy(t2,firstPhase.gravity.gVector);
            a.accelerateBy(tStep*a.mass(),firstPhase.gravity.gVector);
        }
    }
}

/**
* Update of collision list when gravity is changed.
*/
public void update(Observable o, Object arg) {
    if(o instanceof Gravity) {
    for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
        if(a.isStationary() || ((Agent)a.ia).getCollisionPartner().isStationary()) {
            upList(a);
        }
        if(a.isStationary()) {downList(a);}
    }
    findNextCollider();
    }
}
          
//--------------------------------------------------------------

protected void upList(Atom atom) {  
            
    double minCollisionTime = Double.MAX_VALUE;
    Agent aia = (Agent)atom.ia;
    aia.setCollision(minCollisionTime, null, null);
            
    upPairIterator.reset(atom,Iterator.INTRA);
    while(upPairIterator.hasNext()) {
        AtomPair pair = upPairIterator.next();
        PotentialHard potential = (PotentialHard)simulation().getPotential(pair);
        double time = potential.collisionTime(pair);
        if(time < minCollisionTime) {
            minCollisionTime = time;
            aia.setCollision(time,pair.atom2(),potential);  //atom2 is inner loop
        }
    }
        
    //Examine interaction with space
    atomPair.reset(atom,atom);
    double time = spacePotential.collisionTime(atomPair);
    if(time < minCollisionTime) {
        minCollisionTime = time;
        aia.setCollision(time,atom,spacePotential);  //"self-collision" is code for collision with space
    }
}

//--------------------------------------------------------------

protected void downList(Atom atom) {
            
    downPairIterator.reset(atom,Iterator.INTRA);
    while(downPairIterator.hasNext()) {
        AtomPair pair = downPairIterator.next();
        Agent aia = (Agent)pair.atom2().ia;  //atom2 is inner loop
        PotentialHard potential = (PotentialHard)simulation().getPotential(pair);
        double time = potential.collisionTime(pair);
        if(time < aia.getCollisionTime()) {
            aia.setCollision(time,atom,potential);
        }
    }
}
          
//--------------------------------------------------------------

public void initialize() {
    deployAgents();
    Atom.Iterator iterator = firstPhase.iterator.makeAtomIteratorUp();
    iterator.reset(firstPhase.firstAtom());
    while(iterator.hasNext()) {upList(iterator.next());}
    findNextCollider();
}
          
//--------------------------------------------------------------

    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            a.coordinate.scaleMomentum(s);
            ((Agent)a.ia).collisionTime *= rs;
        }
    }
            
    public final Integrator.Agent makeAgent(Atom a) {
        return new Agent(a);
    }
            
    public final static class Agent implements Integrator.Agent {  //need public so to use with instanceof
        public Atom atom;
        public double time0;  //time since last collision
        private double collisionTime = Double.MAX_VALUE; //time to next collision
        private Atom collisionPartner;  //next atom scheduled for collision by atom containing this Agent
        private PotentialHard collisionPotential;  //potential governing interaction between collisionPartner and atom containing this Agent
        public double pAccumulator;  //accumulated momentum, for calculation of pressure
        public double qAccumulator;  //accumulated energy, for calculation of heat transfer

        public Agent(Atom a) {
            atom = a;
            zeroAccumulators();
        }
                  
        public final void setCollision(double time, Atom partner, PotentialHard p) {
            collisionTime = time;
            collisionPartner = partner;
            collisionPotential = p;
        }

        public final void decrementCollisionTime(double interval) {
            collisionTime -= interval;
            time0 += interval;
        }
                  
        public final double getCollisionTime() {return collisionTime;}
        public final Atom getCollisionPartner() {return collisionPartner;}
        public final PotentialHard getCollisionPotential() {return collisionPotential;}

        public final void zeroAccumulators() {
            qAccumulator = 0.0;
            pAccumulator = 0.0;
        }
        
    }
}

