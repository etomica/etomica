package simulate;
import java.util.Observable;

public class IntegratorHard extends Integrator {

private Agent nextCollider;
private AtomPair.Iterator.A iterator;
            
public IntegratorHard() {
    super();
}

public void registerPhaseSpace(PhaseSpace p) {
    super.registerPhaseSpace(p);
    iterator = p.makePairIteratorFull();
}

          
//--------------------------------------------------------------
// steps all particles across time interval tStep, handling
// all intervening collisions

public void doStep(double tStep) {

    if(tStep < nextCollider.getCollisionTime()) {
        advanceAcrossTimeStep(tStep);
        if(isothermal) {
            scaleMomenta(Math.sqrt(this.temperature/firstPhaseSpace.kineticTemperature()));
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
    for(Atom a=firstPhaseSpace.firstAtom(); a!=null; a=a.nextAtom()) {
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
    Atom partner = nextCollider.getCollisionPartner();
    if(partner == null) {
        upList(nextCollider.atom);
        downList(nextCollider.atom);
    }
    else {
//        Atom partnerNextAtom = partner.nextMoleculeFirstAtom();  //put this back in for multiatomic speciesSwitch; also need to do more work with loop below
        Atom partnerNextAtom = partner.nextAtom();
        nextCollider.getCollisionPotential().bump(firstPhaseSpace.makeAtomPair(nextCollider.atom,partner));  //should change this to avoid making pair every time
//        nextCollider.getCollisionPotential().bump(nextCollider.atom,partner);
                
        boolean upListedN = false;
        boolean upListedP = false;
        for(Atom a=firstPhaseSpace.firstAtom(); a!=partnerNextAtom; a=a.nextAtom()) {  //note that nextCollider's or partner's position in linked-list may have been moved by the bump method
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
            
    if(firstPhaseSpace.noGravity) {
        for(Atom a=firstPhaseSpace.firstAtom(); a!=null; a=a.nextAtom()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
//            Space.uEa1Tv1(dr,tStep*a.rm,a.p);
//            a.translate(dr);         //needs modification for nonspherical atom
            a.coordinate.translateToward(a.coordinate.momentum(),tStep*a.rm());
        }
    }
    else {
        double t2 = tStep*tStep;
        for(Atom a=firstPhaseSpace.firstAtom(); a!=null; a=a.nextAtom()) {
            ((Agent)a.ia).decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
//            Space.uEa1Tv1Pa2Tv2(dr,tStep*a.rm,a.p,t2,firstPhaseSpace.gravity.gVector);
            a.coordinate.translateToward(a.coordinate.momentum(),tStep*a.rm());
            a.coordinate.translateToward(firstPhaseSpace.gravity.gVector,t2);
//            Space.uEa1Tv1(dr,tStep*a.mass,firstPhaseSpace.gravity.gVector);
            a.coordinate.accelerateToward(firstPhaseSpace.gravity.gVector,tStep*a.mass());
        }
    }
}

/**
* Update of collision list when gravity is changed.
*/
public void update(Observable o, Object arg) {
    if(o instanceof Gravity) {
    for(Atom a=firstPhaseSpace.firstAtom(); a!=null; a=a.nextAtom()) {
        if(a.isStationary() || ((Agent)a.ia).getCollisionPartner().isStationary()) {
            upList(a);
        }
        if(a.isStationary()) {downList(a);}
    }
    findNextCollider();
    }
}
          
//--------------------------------------------------------------

protected void upList(Atom atom) {  //specific to 2D
            
    double minCollisionTime = Double.MAX_VALUE;
    Agent aia = (Agent)atom.ia;
    if(!atom.isStationary() && atom.type instanceof AtomType.Disk) {  //if mobile, set collision time to time atom takes to move half a box edge
        for(int i=Simulation.D-1; i>=0; i--) {
            double tnew = Math.abs((atom.phaseSpace().dimensions().component(i)-1.0001*((AtomType.Disk)atom.type).diameter())/atom.coordinate.momentum(i));  //assumes range of potential is .le. diameter
            minCollisionTime = (tnew < minCollisionTime) ? tnew : minCollisionTime;
        }
        minCollisionTime *= 0.5*atom.mass();
    }
    aia.setCollision(minCollisionTime, null, null);
            
    Atom nextMoleculeAtom = atom.nextMoleculeFirstAtom();  //first atom on next molecule
    int atomSpeciesIndex = atom.getSpeciesIndex();
            
    //Loop through remaining uplist atoms in this atom's molecule
    if(atom != atom.parentMolecule.lastAtom) {
        Potential1 p1 = firstPhaseSpace.potential1[atomSpeciesIndex];
        iterator.reset(atom.nextAtom(),atom.parentMolecule.lastAtom,atom,atom);
        while(iterator.hasNext()) {
//            for(Atom a=atom.nextAtom(); a!=nextMoleculeAtom; a=a.nextAtom()) {
//                PotentialHard potential = (PotentialHard)p1.getPotential(atom,a);
//                double time = potential.collisionTime(atom,a);
            AtomPair pair = iterator.next();
            PotentialHard potential = (PotentialHard)p1.getPotential(pair.atom1(),pair.atom2());
            double time = potential.collisionTime(pair);
            if(time < minCollisionTime) {
                minCollisionTime = time;
                aia.setCollision(time,pair.atom2(),potential);
            }
        }
    }
            
    //Loop through remaining uplist atoms in firstPhaseSpace
    if(atom.parentMolecule != atom.phaseSpace().lastMolecule()) {
        iterator.reset(atom.nextMoleculeFirstAtom(),atom.phaseSpace().lastAtom(),atom,atom);
        while(iterator.hasNext()) {
            AtomPair pair = iterator.next();
            PotentialHard potential = (PotentialHard)firstPhaseSpace.potential2[pair.atom2().getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,pair.atom2());
            double time = potential.collisionTime(pair);
            if(time < minCollisionTime) {
                minCollisionTime = time;
                aia.setCollision(time,pair.atom2(),potential);  //atom2 is inner loop
            }
        }
    }
//    Potential2[] p2 = firstPhaseSpace.potential2[atomSpeciesIndex];
//    for(Atom a=nextMoleculeAtom; a!=null; a=a.nextAtom()) {
//        Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
//        PotentialHard potential = (PotentialHard)firstPhaseSpace.potential2[a.getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,a);
//        double time = potential.collisionTime(atom,a);
//        if(time < minCollisionTime) {
//            minCollisionTime = time;
//            aia.setCollision(time,a,potential);
//        }
//    }
}

//--------------------------------------------------------------

protected void downList(Atom atom) {
            
    Atom previousMoleculeAtom = atom.getMolecule().firstAtom().previousAtom();
            
    int atomSpeciesIndex = atom.getSpeciesIndex();
            
    //Loop through remaining downlist atoms in this atom's molecule
    if(atom != atom.parentMolecule.firstAtom) {
        iterator.reset(atom.parentMolecule.firstAtom,atom.previousAtom(),atom,atom);
        Potential1 p1 = firstPhaseSpace.potential1[atomSpeciesIndex];
        while(iterator.hasNext()) {
            AtomPair pair = iterator.next();
            Agent aia = (Agent)pair.atom2().ia;  //atom2 is inner loop
            PotentialHard potential = (PotentialHard)p1.getPotential(pair.atom2(),atom);
            double time = potential.collisionTime(pair);
            if(time < aia.getCollisionTime()) {
                aia.setCollision(time,atom,potential);
            }
        }
    }
/*        Potential1 p1 = firstPhaseSpace.potential1[atomSpeciesIndex];
        for(Atom a=atom.previousAtom(); a!=previousMoleculeAtom; a=a.previousAtom()) {
            Agent aia = (Agent)a.ia;
            PotentialHard potential = (PotentialHard)p1.getPotential(a,atom);
            double time = potential.collisionTime(a,atom);
            if(time < aia.getCollisionTime()) {
                aia.setCollision(time,atom,potential);
            }
        }*/
            
    //Loop through remaining downlist atoms in firstPhaseSpace
//   Potential2[] p2 = firstPhaseSpace.potential2[atomSpeciesIndex];
    if(atom.parentMolecule != atom.phaseSpace().firstMolecule()) {
        iterator.reset(atom.phaseSpace().firstAtom(),atom.previousMoleculeLastAtom(),atom,atom);
        while(iterator.hasNext()) {
            AtomPair pair = iterator.next();
            Agent aia = (Agent)pair.atom2().ia;  //atom2 is inner loop
            PotentialHard potential = (PotentialHard)firstPhaseSpace.potential2[pair.atom2().getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,pair.atom2());
    //       Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
            double time = potential.collisionTime(pair);
            if(time < aia.getCollisionTime()) {
                aia.setCollision(time,atom,potential);
            }
        }
    }
 /*   for(Atom a=previousMoleculeAtom; a!=null; a=a.previousAtom()) {
        Agent aia = (Agent)a.ia;
        PotentialHard potential = (PotentialHard)firstPhaseSpace.potential2[a.getSpeciesIndex()][atomSpeciesIndex].getPotential(atom,a);
//       Potential potential = p2[a.getSpeciesIndex()].getPotential(atom,a);
        double time = potential.collisionTime(atom,a);
        if(time < aia.getCollisionTime()) {
            aia.setCollision(time,atom,potential);
        }
    }*/
    
}
          
//--------------------------------------------------------------

public void initialize() {
    deployAgents();
    for(Atom a=firstPhaseSpace.firstAtom(); a!=null; a=a.nextAtom()) {
        upList(a);
    }
    findNextCollider();
}
          
//--------------------------------------------------------------

    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        for(Atom a=firstPhaseSpace.firstAtom(); a!=null; a=a.nextAtom()) {
            a.coordinate.scaleMomentum(s);
            ((Agent)a.ia).collisionTime *= rs;
        }
    }
            
    public Integrator.Agent makeAgent(Atom a) {
        return new Agent(a);
    }
            
    public static class Agent implements Integrator.Agent {  //need public so to use with instanceof
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

