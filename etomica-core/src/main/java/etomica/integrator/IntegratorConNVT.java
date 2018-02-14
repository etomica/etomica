/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Constant NVT Molecular Dynamics Integrator-Constraint Method
 * 
 * Uses a modified version of the Leap Frog Algorithm. 
 * ke is calculated from the unconstrained velocities at time T.
 * A ratio of the setTemperature to the unconstrained temp (as solved from ke),
 * is used to calculate the new constrained velocities at T+Dt/2.  
 * The positions at T+Dt are solved for from the constrained velocities calculated at T+Dt/2.
 *  
 * @author Chris Iacovella
 * @author David Kofke
 */
public final class IntegratorConNVT extends IntegratorMD {

    private static final long serialVersionUID = 1L;
    public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms;
    Vector work, work1, work2, work3, work4;
    double halfTime, mass;

    protected AtomLeafAgentManager<Vector> agentManager;

    public IntegratorConNVT(PotentialMaster potentialMaster, IRandom random,
                            double timeStep, double temperature, Box box) {
        super(potentialMaster,random,timeStep,temperature, box);
        forceSum = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        allAtoms.setIncludeLrc(false);
        work = this.space.makeVector();
        work1 = this.space.makeVector();
        work2 = this.space.makeVector();
        work3 = this.space.makeVector();
       	work4 = this.space.makeVector();
        agentManager = new AtomLeafAgentManager<>(a -> this.space.makeVector(), box);
        forceSum.setAgentManager(agentManager);
    }

  	public final void setTimeStep(double t) {
    	super.setTimeStep(t);
    	halfTime = timeStep/2.0;
  	}
  	
  	private double Temper;
  	public void setTemp(double temperature) {
  	    Temper=temperature;
	}
          
//--------------------------------------------------------------
// steps all particles across time interval tStep

    protected void doStepInternal() {
        super.doStepInternal();

        double dim = space.D();  //get the dimension

  	    //Compute forces on each atom
  	    forceSum.reset();
  	    potentialMaster.calculate(box, allAtoms, forceSum);

  	    //MoveA
        //Advance velocities from T-Dt/2 to T without constraint
        double Free=0.0;
        //degrees of freedom
        Free=((box.getMoleculeList().getMoleculeCount()-1)*dim); 

        double k=0.0;
        double chi;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Vector v = a.getVelocity();

            work1.E(v); //work1 = v
            work2.E(agentManager.getAgent(a));	//work2=F
            work1.PEa1Tv1(halfTime* a.getType().rm(),work2); //work1= p/m + F*Dt2/m = v + F*Dt2/m

            k+=work1.squared();
        }   

        //calculate scaling factor chi
        chi= Math.sqrt(Temper*Free/(mass*k));

        //calculate constrained velbox.getSpace()ocities at T+Dt/2
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Vector force = agentManager.getAgent(a);
            Vector v = a.getVelocity();

            double scale = (2.0*chi-1.0); 
            work3.Ea1Tv1(scale,v); 
            work4.Ea1Tv1(chi* a.getType().rm(),force);
            work4.TE(timeStep);
            work3.PE(work4);
            v.E(work3);
        } 

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Vector r = a.getPosition();
            Vector v = a.getVelocity();

            work.Ea1Tv1(timeStep,v);
            work.PE(r);
            r.E(work);
        }
  	}
}
