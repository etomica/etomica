/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.*;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.Debug;

/**
 * Basic square-well potential.
 * Energy is infinite if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of hard core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential. 
 * 
 * This implementation keeps track of which atom pairs are in the well,
 * avoiding numerical roundoff issues after a collision and avoiding any
 * need to bump atoms apart or together after capture and escape.  If a
 * simulation component moves atoms around externally, this class will become
 * confused and fail.  If such movement is necessary, consider using
 * P2SquareWell instead.
 */
public class P2SquareWellRobust extends Potential2HardSpherical implements AtomLeafAgentManager.AgentSource<AtomArrayList> {

    protected final boolean ignoreOverlap;
    protected final BoxAgentManager<AtomLeafAgentManager<AtomArrayList>> boxWellManager;
    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected Tensor lastCollisionVirialTensor;
    protected double lastEnergyChange;
    protected Vector dv;
    protected AtomLeafAgentManager<AtomArrayList> wellManager;
    protected PotentialMasterList potentialMaster;
    protected AtomType[] atomTypes;
    protected boolean suppressMakeAgent = false;

    public P2SquareWellRobust(Space space) {
        this(space, 1.0, 2.0, 1.0, false);
    }

    public P2SquareWellRobust(Space space, double coreDiameter, double lambda, double epsilon, boolean ignoreOverlap) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
        this.ignoreOverlap = ignoreOverlap;
        boxWellManager = new BoxAgentManager<AtomLeafAgentManager<AtomArrayList>>(null, AtomLeafAgentManager.class);
    }
    
    /**
     * Informs the square-well class what the potential master is (if neighbor
     * listing is used) and the types of atoms interacting with this potential.
     * This helps the potential construct a list of atom pairs that are in the
     * well (atoms that do not interact with this potential will be excluded
     * from the list).
     * 
     * If neighbor listing is not used, null should be passed.  If this method
     * is not called (or if null is passed for both parameters), then the
     * potential will construct a list including all atom pairs (even those
     * that do not interact with this potential).  Also, for a large system,
     * the potential must check all pairs (even pairs with large distances). 
     * 
     * @param potentialMaster
     * @param types
     */
    public void setPotentialMaster(PotentialMasterList potentialMaster, AtomType[] types) {
        this.potentialMaster = potentialMaster;
        this.atomTypes = types;
    }

    public double getRange() {
        return wellDiameter;
    }
    
    /**
     * Implements collision dynamics between two square-well atoms.
     * Includes all possibilities involving collision of hard cores, and collision of wells
     * both approaching and diverging
     */
    public void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double reduced_m = 1.0/(rm0+rm1);
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared) || lambda == 1) {   // Hard-core collision
            if (Debug.ON && !ignoreOverlap && Math.abs(r2 - coreDiameterSquared)/coreDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+coreDiameterSquared);
            }
            if (bij > 0) {
                throw new RuntimeException("really?");
//                lastCollisionVirial = 0;
            }
            else {
                lastCollisionVirial = 2.0*reduced_m*bij;
            }
            lastEnergyChange = 0.0;
        }
        else {    // Well collision
            if (Debug.ON && Math.abs(r2 - wellDiameterSquared)/wellDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+wellDiameterSquared);
            }
            // ke is kinetic energy due to components of velocity
            double ke = bij*bij*reduced_m/(2.0*r2);
            if(wellManager.getAgent(atom0).contains(atom1)) {         // Separating
                if(ke < epsilon) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
                    lastEnergyChange = epsilon;
                    AtomArrayList iList = wellManager.getAgent(atom0);
                    iList.removeAndReplace(iList.indexOf(atom1));
                    AtomArrayList jList = wellManager.getAgent(atom1);
                    jList.removeAndReplace(jList.indexOf(atom0));
                }
            }
            else if(ke > -epsilon) {   // Approach/capture
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m));
                lastEnergyChange = -epsilon;
                wellManager.getAgent(atom0).add(atom1);
                wellManager.getAgent(atom1).add(atom0);
            }
            else {                     // Not enough kinetic energy to overcome square-shoulder
                lastCollisionVirial = 2.0*reduced_m*bij;
                lastEnergyChange = 0.0;
            }
        }
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PEa1Tv1( rm0,dv);
        atom1.getVelocity().PEa1Tv1(-rm1,dv);
        atom0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        atom1.getPosition().PEa1Tv1( falseTime*rm1,dv);

    }//end of bump method

    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }

    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }

    /**
     * Computes next time of collision of two square-well atoms, assuming free-flight kinematics.
     * Collision may occur when cores collides, or when wells first encounter each other on
     * approach, or when they edge of the wells are reached as atoms diverge.
     */
    public double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic coord0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic coord1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        double time = Double.POSITIVE_INFINITY;

        if(wellManager.getAgent(coord0).contains(coord1)) {  // Already inside wells

            if(bij < 0.0) {    // Check for hard-core collision
                if(ignoreOverlap && r2 < coreDiameterSquared) {   // Inside core; collision now
                    return falseTime+0.001*Math.sqrt(dr.squared())/Math.sqrt(v2);
                }

                double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                if(discr > 0) {  // Hard cores collide next
                    time = (-bij - Math.sqrt(discr))/v2;
                }
                else {           // Moving toward each other, but wells collide next
                    discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                    time = (-bij + Math.sqrt(discr))/v2;
                }
            }
            else {           // Moving away from each other, wells collide next
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
                time = (-bij + Math.sqrt(discr))/v2;
            }
        }
        else {              // Outside wells; look for collision at well
            if(bij < 0.0) {
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                if(discr > 0) {
                    time = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.allAtoms(pair)) || time < 0)) {
            System.out.println(pair+" r2 "+r2+" bij "+bij+" time "+(time+falseTime));
        }
        return time + falseTime;
    }

  /**
   * Returns infinity if overlapping, -epsilon if otherwise less than well diameter, or zero if neither.
   */
    public double u(double r2) {
        if (r2 > wellDiameterSquared) return 0.0;
        if (r2 > coreDiameterSquared) return -epsilon;
        return Double.POSITIVE_INFINITY;
    }

    public double energyChange() {return lastEnergyChange;}

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {return coreDiameter;}
    /**
     * Accessor method for core diameter.
     * Well diameter is defined as a multiple (lambda) of this, and is updated when core diameter is changed
     */
    public void setCoreDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        if (coreDiameter==0) wellDiameter = lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getCoreDiameterDimension() {return Length.DIMENSION;}

    /**
     * Accessor method for well-diameter multiplier.
     */
    public double getLambda() {return lambda;}
    /**
     * Accessor method for well-diameter multiplier.
     * Well diameter is defined as this multiple of core diameter, and is updated when 
     * this is changed
     */
    public void setLambda(double lam) {
        if (lam < 1.0) throw new IllegalArgumentException("Square-well lambda must be greater than 1.0");
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        if (coreDiameter==0) wellDiameter = lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getLambdaDimension() {return Null.DIMENSION;}

    /**
     * Accessor method for depth of well
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Accessor method for depth of well
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
    
    public void setBox(Box box) {
        super.setBox(box);
        wellManager = boxWellManager.getAgent(box);
        if (wellManager == null) {
            wellManager = makeAgent(box);
            boxWellManager.setAgent(box, wellManager);
        }
    }
    
    protected void handleNewAtomNbr(IAtom a, AtomLeafAgentManager<AtomArrayList> agentManager, NeighborListManager nbrListManager, AtomArrayList iList) {
        if (iList == null) {
            iList = agentManager.getAgent(a);
        }
        PotentialArray potentialArray = potentialMaster.getRangedPotentials(a.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        int u;
        for (u=0; u<potentials.length; u++) {
            if (potentials[u] == this) break;
        }
        if (u == potentials.length) return;
        IAtomList upList = nbrListManager.getUpList(a)[u];
        for (int jj=0; jj<upList.getAtomCount(); jj++) {
            IAtom jAtom = upList.getAtom(jj);

            dr.Ev1Mv2(a.getPosition(), jAtom.getPosition());
            boundary.nearestImage(dr);

            double r2 = dr.squared();
            if (r2 < wellDiameterSquared) {
                iList.add(jAtom);
                agentManager.getAgent(jAtom).add(a);
            }
        }
    }
    
    protected void handleNewAtomNoNbr(IAtom a, AtomLeafAgentManager<AtomArrayList> agentManager, IAtomList leafList, AtomArrayList iList) {
        if (iList == null) {
            iList = agentManager.getAgent(a);
        }
        if (atomTypes != null && a.getType() != atomTypes[0] && a.getType() != atomTypes[1]) return;
        for (int j=a.getLeafIndex()+1; j<leafList.getAtomCount(); j++) {
            IAtom jAtom = leafList.getAtom(j);
            if (atomTypes != null && jAtom.getType() != atomTypes[0] && jAtom.getType() != atomTypes[1]) continue;
            
            dr.Ev1Mv2(a.getPosition(), jAtom.getPosition());
            boundary.nearestImage(dr);

            double r2 = dr.squared();
            if (r2 < wellDiameterSquared) {
                iList.add(jAtom);
                agentManager.getAgent(jAtom).add(a);
            }
        }
    }
    
    public AtomLeafAgentManager<AtomArrayList> makeAgent(Box box) {
        suppressMakeAgent = true;
        AtomLeafAgentManager<AtomArrayList> foo = new AtomLeafAgentManager<AtomArrayList>(this, box, AtomArrayList.class);
        IAtomList leafList = box.getLeafList();
        if (potentialMaster != null) {
            NeighborListManager nbrListManager = potentialMaster.getNeighborManager(box);
            
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtom a = leafList.getAtom(i);
                handleNewAtomNbr(a, foo, nbrListManager, null);
            }
        }
        else {
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtom a = leafList.getAtom(i);
                handleNewAtomNoNbr(a, foo, leafList, null);
            }
        }
        suppressMakeAgent = false;
        return foo;
    }
    
    public AtomArrayList makeAgent(IAtom atom, Box agentBox) {
        if (suppressMakeAgent) {
            // allow [box] makeAgent to find all interacting pairs
            return new AtomArrayList();
        }
        AtomArrayList rv = new AtomArrayList();
        if (potentialMaster != null) {
            handleNewAtomNbr(atom, boxWellManager.getAgent(agentBox), potentialMaster.getNeighborManager(agentBox), rv);
        }
        else {
            handleNewAtomNoNbr(atom, boxWellManager.getAgent(agentBox), agentBox.getLeafList(), rv);
        }
        return rv;
    }

    public void releaseAgent(AtomArrayList iList, IAtom atom, Box agentBox) {
        // atom is going away.  remove it from all of its neighbor's lists
        AtomLeafAgentManager<AtomArrayList> agentManager = boxWellManager.getAgent(agentBox);
        for (int j=0; j<iList.getAtomCount(); j++) {
            AtomArrayList jList = agentManager.getAgent(iList.getAtom(j));
            jList.removeAndReplace(jList.indexOf(atom));
        }
        iList.clear();
    }
}
  
