package simulate;
import java.beans.*;

public class IntegratorHardSoft extends Integrator {

  SpeciesElement nextCollider;
  double forceStep;
  int colCount=0;
  private transient final double dr[] = new double[Space.D];
  private transient final double dp[] = new double[Space.D];
  private transient final double v12[] = new double[Space.D];
 
  public IntegratorHardSoft() {
    super();
    setForceStep(0.001);
  }

//--------------------------------------------------------------

  public double getForceStep() {return forceStep;}
  public void setForceStep(double step) {forceStep = step;}
  
//--------------------------------------------------------------
// steps all particles across time interval tStep, handling
// all intervening collisions

  public void doStep(double tStep) {

    if(tStep < nextCollider.getCollisionTime()) {
      advanceAcrossTimeStep(tStep);
      return;
    }

    double tStepNew = tStep - nextCollider.getCollisionTime();
    advanceToCollision();
    doStep(tStepNew);
    return;
  }

//--------------------------------------------------------------
 
  private void findNextCollider() {
    //find next collision pair by looking for minimum collisionTime
    double minCollisionTime = Double.MAX_VALUE;
    for(SpeciesElement e=phase.firstElement; e!=null; e=e.getNext()) {
      if(e.getCollisionTime() < minCollisionTime) {
        minCollisionTime = e.getCollisionTime();
        nextCollider = e;}
    }
}
//--------------------------------------------------------------
// advances to next collision, applies collision dynamics to involved 
// species and updates collision time/partners

  private void advanceToCollision() {
    advanceAcrossTimeStep(nextCollider.getCollisionTime());
    SpeciesElement partner = nextCollider.getCollisionPartner();
    this.bump(nextCollider,partner);
    for(SpeciesElement e=phase.firstElement; e!=partner; e=e.getNext()) {
        if(e.getCollisionPartner()==nextCollider || e.getCollisionPartner()==partner || e==nextCollider) {
            upList(e);
        }
    }
    upList(partner);
    downList(nextCollider);
    downList(partner);

    findNextCollider();
    
    phase.updatedKineticEnergy = false;
    phase.updatedPotentialEnergy = false;
}


// --------- advances all particle coordinates and velocities by tStep --------

  private void advanceAcrossTimeStep(double tStep) {  //modified from integratorHard method
    for(SpeciesElement e=phase.firstElement; e!=null; e=e.getNext()) {e.decrementCollisionTime(tStep);}
    double halfTStep = 0.5*tStep;
    for(int j = 0; j < phase.speciesVector.size(); j++) {
        Species species = (Species)phase.speciesVector.elementAt(j);
        if(species instanceof SpeciesDisk) {
            SpeciesDisk speciesDisk = (SpeciesDisk)species;
            double rm = speciesDisk.rm;
            double a1 = speciesDisk.rm*tStep;    // rm*tStep
            double a2 = a1*halfTStep;            // 0.5*rm*tStep^2
            for (int i=speciesDisk.nElements; --i>=0; ) {
               MoleculeAtomic m = speciesDisk.molecule[i];
               Space.uEa1Tv1Pa2Tv2(dr, a1, m.p, a2, m.f);
               Space.uEa1Tv1(dp, tStep, m.f);
               m.translate(dr);   //dr = rm*(tStep*p + 0.5*tStep^2*f
               m.accelerate(dp);  //dp = tStep*f
            }
        }
    }
  }
  
//--------------------------------------------------------------

  private void upList(SpeciesElement element) {

    if(element==phase.lastElement) {return;}    //lastElement has no collision partner

    SpeciesElement partner=element.getCollisionPartner();
    double minCollisionTime = collisionTime(element,partner);  //min collision time likely to be with previous partner
    for(SpeciesElement e=element.getNext(); e!=null; e=e.getNext()) {
//        if(wontCollide(element,e,minCollisionTime)) {continue;}
        
        double pairCollisionTime = collisionTime(element,e);
        if(pairCollisionTime < minCollisionTime) {
            minCollisionTime = pairCollisionTime;
            partner = e;
        }
    }
    element.setCollisionTime(minCollisionTime);
    
    if(partner != element.getCollisionPartner()) {          //found a new collision partner
        updateCollisionPartner(element,partner);
    }
  }
  
//--------------------------------------------------------------

  private void updateCollisionPartner(SpeciesElement collider, SpeciesElement partner) {
     if(collider != nextCollider) {this.bump(collider,collider.getCollisionPartner());}   //advance velocities with old partner
     collider.setCollisionPartner(partner);               //set new partner
     double[] force = phase.potential[collider.getSpeciesIndex()][partner.getSpeciesIndex()].force(collider,partner);
     MoleculeAtomic disk1 = (MoleculeAtomic)collider;     //can remove this if speciesElement updated
     disk1.setPartnerForce(force);
  }

//--------------------------------------------------------------

  private void downList(SpeciesElement element) {
    
//    for(SpeciesElement e=phase.firstElement; e!=element; e=e.getNext()) {
    for(SpeciesElement e=element.getPrevious(); e!=null; e=e.getPrevious()) {
//        if(wontCollide(e,element,e.getCollisionTime()) {continue;}
        
        double pairCollisionTime = collisionTime(element,e);
        if(pairCollisionTime < e.getCollisionTime()) {
            e.setCollisionTime(pairCollisionTime);
            updateCollisionPartner(e,element);
        }
    }
  }

//--------------------------------------------------------------

//  private boolean wontCollide(SpeciesElement e1, SpeciesElement e2, double time) {
    
    


//--------------------------------------------------------------
  
      public double collisionTime(SpeciesElement element1, SpeciesElement element2) {
        MoleculeAtomic disk1 = (MoleculeAtomic)element1;
        MoleculeAtomic disk2 = (MoleculeAtomic)element2;
        PairInteraction pair = phase.potential[element1.getSpeciesIndex()][element2.getSpeciesIndex()].computePairInteraction(element1,element2);
        double fTot = Math.sqrt(Space.v1Dv2(pair.force,pair.force)) + forceStep;
        double delr = Math.abs(Math.sqrt(pair.rSquared) - Math.pow(12.0/1.e12/fTot,1.0/13.0));
        Space.uEa1Tv1Ma2Tv2(v12,disk2.rm,disk2.p,disk1.rm,disk1.p);  //v2-v1 = (p/m)2 - (p/m)1
 //       double r2 = phase.space.rr(disk1,disk2);
 //       double delr2 = forceStep/phase.potential[element1.getSpeciesIndex()][element2.getSpeciesIndex()].dfdr(r2);  
 //       delr2 = delr2*delr2;
        double speed = Math.sqrt(Space.v1Dv2(v12,v12));
        double time = delr/speed;
        double vDotf = disk1.rm*(pair.force[0]*(disk2.p[0]-disk1.p[0])+pair.force[1]*(disk2.p[1]-disk1.p[1]));
        speed = Math.sqrt(speed*speed + 0.5*(vDotf*time + fTot*fTot*time*time));
        time = delr/speed;
        return time;
    }
    
//--------------------------------------------------------------
    
    public void bump(SpeciesElement e1, SpeciesElement e2)
    {
        MoleculeAtomic collider = (MoleculeAtomic)e1;
        MoleculeAtomic partner = (MoleculeAtomic)e2;
        
        double oldForce[] = collider.getPartnerForce();
        e2.subtractForce(oldForce);    //subtract force of 1 on 2
        e1.addForce(oldForce);         //subtract minus force of 1 on 2
        PairInteraction pair = phase.potential[e1.getSpeciesIndex()][e2.getSpeciesIndex()].computePairInteraction(e1,e2);
        e2.addForce(pair.force);       //add force of 1 on 2
        e1.subtractForce(pair.force);  //add minus force
        Space.uEa1T_v1Mv2_(dp, 0.5*collider.time0, oldForce, pair.force);
        collider.accelerate(dp);
        partner.decelerate(dp);
        collider.setPartnerForce(pair.force);  //in case collider has same partner won't have to compute force again
//        colCount++;
//        System.out.println(colCount+" "+pair.rSquared);
    }                                          //also sets collider.time0 to zero

//--------------------------------------------------------------
  
  public void initialize() {
    phase.updateNeighbors();
    phase.updateForces();
    for(SpeciesElement e=phase.firstElement; e!=null; e=e.getNext()) {
        upListInitialize(e);
    }
    findNextCollider();
  }
  
//--------------------------------------------------------------

  private void upListInitialize(SpeciesElement element) {
    SpeciesElement part = null;
    double pairCollisionTime;
    
    int elementSpeciesIndex = element.getSpeciesIndex();

    double minCollisionTime = Double.MAX_VALUE;
    for(SpeciesElement e=element.getNext(); e!=null; e=e.getNext()) {
//        if(wontCollide(element,e,minCollisionTime)) {continue;}
        
        pairCollisionTime = collisionTime(element,e);
        if(pairCollisionTime < minCollisionTime) {
            minCollisionTime = pairCollisionTime;
            part = e;
        }
    }
    element.setCollisionTime(minCollisionTime);
    if(part!=null) {
    element.setCollisionPartner(part);
    double[] force = phase.potential[element.getSpeciesIndex()][part.getSpeciesIndex()].force(element,part);
    MoleculeAtomic disk1 = (MoleculeAtomic)element;     //can remove this if speciesElement updated
    disk1.setPartnerForce(force);
    }
  }
  
}
//--------------------------------------------------------------

