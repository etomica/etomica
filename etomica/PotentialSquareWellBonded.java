package etomica;
import etomica.units.Dimension;

/**
 * Similar to square-well potential, but considers and alters bonding states with collisions.
 * Each atom may bind, in the form of the square-well attraction, to one other atom.  Two atoms
 * approaching each other with this potential may interact as follows:
 * <ul>
 *   <li>If neither is bound, or if they are bound to each other, 
 *       they will interact as square-well atoms.
 *   <li>If one or both is bound, each to another atom, they will bind (form square-well
 *       attraction) with each other, only if the following are satisfied:
 *       <ol>
 *         <li>ke + epsilon > epsilonA + epsilonB
 *         <li>ke > epsilonA + epsilonB + barrierA + barrierB
 *       </ol>
 * </ul>
 * where <i>ke</i> is the kinetic energy of the two atoms for motion directed toward
 * each other; <i>epsilon</i> is the square-well depth for this pair; <i>epsilonA</i>
 * and <i>epsilonB</i> are the square-well depths for each atom with its current
 * binding partner (if it has one); and <i>barrierA</i> and <i>barrierB</i> are the
 * extra barrier potentials (zero by default) for unbinding of each atom with its current partner,
 * if it has one.  Otherwise they repel as hard spheres of the same diameter as the diameter
 * of the square well.  If they bind to each other, their interactions with their current
 * binding partners (if any) become pure hard-sphere repulsion.<br>
 * The first criterion ensures that the total energy after the reaction is sufficient
 * to result in a positive kinetic energy, while the second criterion ensures that the
 * kinetic energy is sufficient to surmount the binding potential plus added activation energy.
 *
 * @author David Kofke
 */
public class PotentialSquareWellBonded extends PotentialSquareWell implements Potential.Reactive {
  
  public String getVersion() {return "PotentialSquareWellBonded:01.05.25/"+super.getVersion();}

  private Potential.Reactive.BondChangeData[] bondData;
  private double barrier;
  
  public PotentialSquareWellBonded() {
    super();
    makeBondData();
  }
  public PotentialSquareWellBonded(Simulation sim) {
    super(sim);
    makeBondData();
  }
  public PotentialSquareWellBonded(double coreDiameter, double lambda, double epsilon) {
    super(coreDiameter, lambda, epsilon);
    makeBondData();
  }
  
  public PotentialSquareWellBonded(Simulation sim, double coreDiameter, double lambda, double epsilon) {
    super(sim, coreDiameter, lambda, epsilon);
    makeBondData();
  }
  
  //constructs 2 BondChangeData classes, one for each atom in possible bondchange event
  private void makeBondData() {
    bondData = new Potential.Reactive.BondChangeData[2];
    for(int i=0; i<2; i++) {
        bondData[i] = new Potential.Reactive.BondChangeData();
        bondData[i].oldPartners = new Atom[1];
        bondData[i].newPartners = new Atom[1];
    }
  }
  
  public void setBarrier(double b) {barrier = b;}
  public double getBarrier() {return barrier;}
  public Dimension getBarrierDimension() {return Dimension.ENERGY;}
  
  public Potential.Reactive.BondChangeData[] getBondChangeData() {return bondData;}

/**
 * Computes next time of collision of the two atoms, assuming free-flight kinematics.
 */
  public double collisionTime(AtomPair pair) {

    Atom.Linker a1L = pair.atom1().atomLink[0];
    Atom a1Partner = (a1L != null) ? a1L.atom() : null;
    
    //inside well but not mutually bonded;  collide now if approaching
    if(a1Partner != pair.atom2() && pair.r2() < wellDiameterSquared) 
        return (pair.vDotr() < 0.0) ? 0.0 : Double.MAX_VALUE;
    //mutually bonded, or outside well; collide as SW
    else return super.collisionTime(pair);
  }

/**
 * Implements collision dynamics between the atom pair.
 */
  public void bump(AtomPair pair) {
    double eps = 1.0e-6;
    double r2 = pair.r2();
    double bij = pair.vDotr();
    // ke is kinetic energy due to components of velocity
    double reduced_m = 1.0/(pair.atom1().rm() + pair.atom2().rm());
    dr.E(pair.dr());
    double ke = bij*bij*reduced_m/(2.0*r2);
    double s, r2New;
    
    Atom a1 = pair.atom1();
    Atom a2 = pair.atom2();
    Atom.Linker a1L = a1.atomLink[0];
    Atom.Linker a2L = a2.atomLink[0];
    Atom a1Partner = (a1L != null) ? a1L.atom() : null;
    Atom a2Partner = (a2L != null) ? a2L.atom() : null;
    
    if(a1Partner == a2) //atoms are bonded to each other
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
            lastCollisionVirial = 2.0*reduced_m*bij;
            r2New = r2;
            bondData[0].atom = null;
        }
        else {    // Well collision; assume separating because mutually bonded
	        if(ke < epsilon) {     // Not enough kinetic energy to escape
	            lastCollisionVirial = 2.0*reduced_m*bij;
	            r2New = (1-eps)*wellDiameterSquared; 
	            bondData[0].atom = null;
	        }
	        else {                 // Escape
	            lastCollisionVirial = reduced_m*bij - Math.sqrt(2.*reduced_m*r2*(ke-epsilon));
	            r2New = (1+eps)*wellDiameterSquared;
	            bondData[0].atom = a1;
	            bondData[0].oldPartners[0] = a2;
	            bondData[0].newPartners[0] = null;
	            bondData[1].atom = a2;
	            bondData[1].oldPartners[0] = a1;
	            bondData[1].newPartners[0] = null;
	   //         bondData[2].atom = null;
	            a1.atomLink[0] = null;
	            a2.atomLink[0] = null;
	        }//end if(ke < epsilon)
        }//end if(2*r2...
    else {                  //not bonded to each other
        if(r2 < (1-eps)*wellDiameterSquared) {   //already inside well; push away
            lastCollisionVirial = 2.0*reduced_m*bij;
            r2New = r2;
            bondData[0].atom = null;
        }
        else { //well collision; decide whether to bond or have hard repulsion
            double deltaE = ke + epsilon;
            double barrier = 0.0;
            if(a1Partner != null) {
                PotentialSquareWellBonded potl = ((PotentialSquareWellBonded)parentSimulation().getPotential(a1,a1Partner));
                deltaE -= potl.getEpsilon();
            //    barrier += potl.getBarrier();
                barrier += potl.getEpsilon() + potl.getBarrier();
            }
            if(a2Partner != null) {
                PotentialSquareWellBonded potl = ((PotentialSquareWellBonded)parentSimulation().getPotential(a2,a2Partner));
                deltaE -= potl.getEpsilon();
            //    barrier += potl.getBarrier();
                barrier += potl.getEpsilon() + potl.getBarrier();
            }
            if(deltaE < 0.0 || ke < barrier) { //insufficient energy to break other bonds; repel
                lastCollisionVirial = 2.0*reduced_m*bij;
                r2New = (1-eps)*wellDiameterSquared;
                bondData[0].atom = null;
            }
            else { //sufficient energy to break other bonds
                lastCollisionVirial = reduced_m*bij + Math.sqrt(2.*reduced_m*r2*deltaE);
                r2New = (1-eps)*wellDiameterSquared;
        
                bondData[0].atom = a1;
                bondData[0].newPartners[0] = a2;
                bondData[0].oldPartners[0] = a1Partner;

                bondData[1].atom = a2;
                bondData[1].newPartners[0] = a1;
                bondData[1].oldPartners[0] = a2Partner;
            
              //redundant bonddata commented out
            //    int idx = 2;
            //    bondData[2].atom = null;
            //    bondData[3].atom = null;
                if(a1Partner != null) {
                    a1Partner.atomLink[0] = null;
            //        bondData[2].atom = a1Partner;
            //       bondData[2].oldPartners[0] = a1;
            //        bondData[2].newPartners[0] = null;
            //        idx++;
                }
                if(a2Partner != null) {
                    a2Partner.atomLink[0] = null;
            //        bondData[idx].atom = a2Partner;
            //        bondData[idx].oldPartners[0] = a2;
            //        bondData[idx].newPartners[0] = null;
                } 
                a1.atomLink[0] = new Atom.Linker(a2);
                a2.atomLink[0] = new Atom.Linker(a1);
            }//end if(deltaE < 0.0) else
        }//end if(r2 < wellDiameterSquared) else
    } //end if(a1Partner == a2) else
    lastCollisionVirialr2 = lastCollisionVirial/r2;
    pair.cPair.push(lastCollisionVirialr2);
    if(r2New != r2) pair.cPair.setSeparation(r2New);
  }//end of bump
}//end of class
    