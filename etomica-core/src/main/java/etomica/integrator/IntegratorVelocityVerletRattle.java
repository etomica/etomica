/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Integrator implementing RATTLE algorithm.
 *
 * @author Andrew Schultz
 */
public class IntegratorVelocityVerletRattle extends IntegratorVelocityVerletShake {

    private static final long serialVersionUID = 1L;
    protected final Vector dv;

    public IntegratorVelocityVerletRattle(Simulation sim, PotentialMaster potentialMaster, Space _space) {
        this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }
    
    

	protected void randomizeMomenta() {
		super.randomizeMomenta();
		
		if(true)return;//make the initial translation or rotation velocity zero
		IMoleculeList molecules = box.getMoleculeList();
		  for (int i = 0; i<molecules.getMoleculeCount(); i++){
	            IMolecule molecule = molecules.getMolecule(i);
	        	IAtomList leafList = molecule.getChildList();
	        	Vector totalMom =  space.makeVector();
	        	Vector totalV = space.makeVector();
	        	double totalMass = 0.0;
//	        	for(int j = 0; j < leafList.getAtomCount() -1; j++){
//	        		IAtomKinetic a = (IAtomKinetic)leafList.getAtom(j);
//	        		IVectorMutable v = a.getVelocity();
//	        		double mass = a.getType().getMass();
//	        		totalMass += mass;
//		        	totalMom.PEa1Tv1(mass, v);
//    			}
	//	        totalV.E(totalMom);
	//	        double totalMassInverse = 1.0/totalMass;
	//	        totalV.TE(totalMassInverse);
//	            for(int j = 0; j < leafList.getAtomCount()-1; j++){
//	        	   IAtomKinetic a = (IAtomKinetic)leafList.getAtom(j);
//	        	   a.getVelocity().ME(totalV);//make it with no translation velocity, only rotation velocity 
//	           }
	           
	           //make the initial only two hydrogen have the initial rotation velocity
//	           ((IAtomKinetic)leafList.getAtom(2)).getVelocity().E(0);
//	           ((IAtomKinetic)leafList.getAtom(3)).getVelocity().E(0);
//	            IVectorMutable h1 = leafList.getAtom(0).getPosition();
//           		IVectorMutable h2 = leafList.getAtom(1).getPosition();
//           		IVectorMutable o = leafList.getAtom(2).getPosition();
//				IVectorMutable m = leafList.getAtom(3).getPosition();
//				IVectorMutable h1velocity = ((IAtomKinetic)leafList.getAtom(0)).getVelocity();
//				IVectorMutable h2velocity = ((IAtomKinetic)leafList.getAtom(1)).getVelocity();
//				IVectorMutable h1h2 = space.makeVector();
//	        	IVectorMutable om = space.makeVector();
//				h1h2.Ev1Mv2(h2, h1);
//				om.Ev1Mv2(m, o);
//				h1h2.normalize();
//				om.normalize();
//				h1velocity.PEa1Tv1(-1.0*h1velocity.dot(om),om);
//				h1velocity.PEa1Tv1(-1.0*h1velocity.dot(h1h2), h1h2);
//				h2velocity.PEa1Tv1(-1.0*h2velocity.dot(om),om);
//				h2velocity.PEa1Tv1(-1.0*h2velocity.dot(h1h2), h1h2);
//				totalV.E(h1velocity);
//				totalV.PE(h2velocity);
//				h1velocity.PEa1Tv1(-0.5, totalV);
//				h2velocity.PEa1Tv1(-0.5, totalV);
//				System.out.println("h1velocity = " + h1velocity );
//				System.out.println("h2velocity = " + h2velocity );
	           
	           //make the initial velocity only in oh1h2 plane
//	            IVectorMutable h1 = leafList.getAtom(0).getPosition();
//           		IVectorMutable h2 = leafList.getAtom(1).getPosition();
//           		IVectorMutable o = leafList.getAtom(2).getPosition();
//				IVectorMutable m = leafList.getAtom(3).getPosition();
//				IVectorMutable h1velocity = ((IAtomKinetic)leafList.getAtom(0)).getVelocity();
//				IVectorMutable h2velocity = ((IAtomKinetic)leafList.getAtom(1)).getVelocity();
//				IVectorMutable ovelocity = ((IAtomKinetic)leafList.getAtom(2)).getVelocity();
//				IVectorMutable mvelocity = ((IAtomKinetic)leafList.getAtom(3)).getVelocity();
//				IVectorMutable h1h2 = space.makeVector();
//	        	IVectorMutable om = space.makeVector();
//	        	h1h2.Ev1Mv2(h2,h1);
//	        	om.Ev1Mv2(m, o);
//	        	IVectorMutable ph1h2 = space.makeVector();//ph1h2 is the vector perpendicular to oh1h2 plane
//	        	ph1h2.E(h1h2);
//	        	ph1h2.XE(om);
//	        	ph1h2.normalize();
//	        	h1velocity.PEa1Tv1(-1.0*h1velocity.dot(ph1h2),ph1h2);
//	        	h2velocity.PEa1Tv1(-1.0*h2velocity.dot(ph1h2),ph1h2);
//	        	ovelocity.PEa1Tv1(-1.0*ovelocity.dot(ph1h2),ph1h2);
//	        	mvelocity.PEa1Tv1(-1.0*mvelocity.dot(ph1h2),ph1h2);
	        	
	        	for(int j = 0; j < leafList.getAtomCount() -1; j++){
	        		IAtomKinetic a = (IAtomKinetic)leafList.getAtom(j);
	        		Vector v = a.getVelocity();
	        		double mass = a.getType().getMass();
	        		totalMass += mass;
		        	totalMom.PEa1Tv1(mass, v);
		    	}
		        totalV.E(totalMom);
		        double totalMassInverse = 1.0/totalMass;
		        totalV.TE(totalMassInverse);
		        
		             for(int j = 0; j < leafList.getAtomCount()-1; j++){
		        	   IAtomKinetic a = (IAtomKinetic)leafList.getAtom(j);
//		        	   a.getVelocity().E(totalV);//make it translation velocity 
		        	   a.getVelocity().ME(totalV);//make it  rotation velocity 
		        }
		  }
	}

    public IntegratorVelocityVerletRattle(Simulation sim, PotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature, Space _space) {
        super(sim, potentialMaster,random,timeStep,temperature, _space);
        dv = space.makeVector();
    }

    protected void doStepInternal() {
        currentTime += timeStep;

        IMoleculeList molecules = box.getMoleculeList();

        // RATTLE
        int numBondedMolecules = 0;
        int numIterations = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints != null) {
                numBondedMolecules++;
                IAtomList childList = molecule.getChildList();
                Boundary boundary = box.getBoundary();

                if (drOld.length < bondConstraints.bondedAtoms.length) {
                    Vector[] newDrOld = new Vector[bondConstraints.bondedAtoms.length];
                    System.arraycopy(drOld, 0, newDrOld, 0, drOld.length);
                    for (int j=drOld.length; j<newDrOld.length; j++) {
                        newDrOld[j] = space.makeVector();
                    }
                    drOld = newDrOld;
                }

                for (int j=0; j<bondConstraints.bondedAtoms.length; j++) {
                    IAtom atom0 = childList.getAtom(bondConstraints.bondedAtoms[j][0]);
                    IAtom atom1 = childList.getAtom(bondConstraints.bondedAtoms[j][1]);
                    drOld[j].Ev1Mv2(atom1.getPosition(), atom0.getPosition());
                    boundary.nearestImage(drOld[j]);
                }
            }
            
            IAtomList leafList = molecule.getChildList();
            int nLeaf = leafList.getAtomCount();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                MyAgent agent = agentManager.getAgent(a);
                Vector r = a.getPosition();
                Vector v = a.getVelocity();
                if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                    System.out.println("first "+a+" r="+r+", v="+v+", f="+agent.force);
                }
                if  (a.getType().getMass() != 0) {
                    v.PEa1Tv1(0.5*timeStep*a.getType().rm(),agent.force);  // p += f(old)*dt/2
                }
                r.PEa1Tv1(timeStep,v);         // r += p*dt/m
            }

            IAtomList childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            Boundary boundary = box.getBoundary();
            double[] bondLengths = bondConstraints.bondLengths;

            if (childList.getAtomCount() > moved.length) {
                moved = new boolean[2][childList.getAtomCount()];
            }
            for (int j=0; j<childList.getAtomCount(); j++) {
                moved[1][j] = true;
            }
            
            for (int iter = 0; iter<maxIterations; iter++) {
                numIterations++;
                boolean success = true;
                for (int j=0; j<childList.getAtomCount(); j++) {
                    moved[0][j] = moved[1][j];
                    moved[1][j] = false;
                }
                for (int j=0; j<bondConstraints.bondedAtoms.length; j++) {
                    int iAtom1 = bondedAtoms[j][0];
                    int iAtom2 = bondedAtoms[j][1];
                    if (!moved[0][iAtom1] && !moved[0][iAtom2]) {
                        continue;
                    }
                    IAtomKinetic atom1 = (IAtomKinetic)childList.getAtom(iAtom1);
                    IAtomKinetic atom2 = (IAtomKinetic)childList.getAtom(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
//                    System.out.println(iter+" "+j+" "+atom1.getVelocity());
//                    System.out.println(iter+" "+j+" "+atom2.getVelocity());
                    double dr2 = dr.squared();
                    double bl2 = bondLengths[j]*bondLengths[j];
//                    System.out.println(Math.sqrt(dr2)+" "+Math.sqrt(bl2));
                    double diffSq = bl2 - dr2;
                    if (Math.abs(diffSq/bl2) > shakeTol) {
                        double mass1 = atom1.getType().getMass();
                        double mass2 = atom2.getType().getMass();
                        double rMass = 1.0/mass1 + 1.0/mass2;
                        double drDotDrOld = dr.dot(drOld[j]);
                        if  (drDotDrOld / bl2 < 0.1) {
                            System.out.println("molecule "+i);
                            System.out.println("dr "+dr);
                            System.out.println("drOld "+drOld[j]);
                            System.out.println("drDotDrOld "+drDotDrOld);
                            throw new RuntimeException("oops");
                        }
                        double gab = diffSq / (2.0 * rMass * drDotDrOld);
                        atom2.getPosition().PEa1Tv1( gab/mass2, drOld[j]);
                        atom1.getPosition().PEa1Tv1(-gab/mass1, drOld[j]);

                        gab /= timeStep;
                        
                        atom2.getVelocity().PEa1Tv1( gab/mass2, drOld[j]);
                        atom1.getVelocity().PEa1Tv1(-gab/mass1, drOld[j]);
                        
//                        System.out.println("new "+atom1.getVelocity());
//                        System.out.println("new "+atom2.getVelocity());
                        
                        moved[1][iAtom1] = true;
                        moved[1][iAtom2] = true;
                        success = false;
                    }
                }
                if (success) {
                    break;
                }
                if (iter == maxIterations-1) {
                    System.err.println("failed to converge in shake for molecule "+i);
                }
            }
        }

        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            bondConstraints.relaxMolecule(molecule);
        }

        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, allAtoms, forceSum);

        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            bondConstraints.redistributeForces(molecule, agentManager);
        }

        //Finish integration step
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            Vector velocity = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("second "+a+" v="+velocity+", f="+agentManager.getAgent(a).force);
            }
            if (a.getType().getMass() != 0) {
                velocity.PEa1Tv1(0.5*timeStep*a.getType().rm(),agentManager.getAgent(a).force);  //p += f(new)*dt/2
            }
        }

        /*
         * Rattle Part II
         */
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            
            IAtomList childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            Boundary boundary = box.getBoundary();
            double[] bondLengths = bondConstraints.bondLengths;

            for (int j=0; j<childList.getAtomCount(); j++) {
                moved[1][j] = true;
            }
            
            for (int iter = 0; iter<maxIterations; iter++) {
                numIterations++;
                boolean success = true;
                for (int j=0; j<childList.getAtomCount(); j++) {
                    moved[0][j] = moved[1][j];
                    moved[1][j] = false;
                }
                for (int j=0; j<bondConstraints.bondedAtoms.length; j++) {
                    int iAtom1 = bondedAtoms[j][0];
                    int iAtom2 = bondedAtoms[j][1];
                    if (!moved[0][iAtom1] && !moved[0][iAtom2]) {
                        continue;
                    }
                    IAtomKinetic atom1 = (IAtomKinetic)childList.getAtom(iAtom1);
                    IAtomKinetic atom2 = (IAtomKinetic)childList.getAtom(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
                    dv.Ev1Mv2(atom2.getVelocity(), atom1.getVelocity());
                    double drdotdv = dr.dot(dv);
                    double mass1 = atom1.getType().getMass();
                    double mass2 = atom2.getType().getMass();
                    double bl2 = bondLengths[j]*bondLengths[j];
                    double g = -drdotdv / ((1.0/mass1+1.0/mass2) * bl2);
                    if (Math.abs(g) > shakeTol) {
                        dr.TE(g);
                        
                        atom2.getVelocity().PEa1Tv1( 1.0/(mass2), dr);
                        atom1.getVelocity().PEa1Tv1(-1.0/(mass1), dr);
                        
                        moved[1][iAtom1] = true;
                        moved[1][iAtom2] = true;
                        success = false;
                    }
                }
                if (success) {
                    break;
                }
                if (iter == maxIterations-1) {
                    System.err.println("failed to converge in rattle for molecule "+i);
                }
            }
        }

        if(isothermal) {
            doThermostatInternal();
        }

        if (printInterval > 0 && stepCount%printInterval == 0) {
            double PE = meterPE.getDataAsScalar();
            double KE = meterKE.getDataAsScalar();
            int moleculeCount = box.getMoleculeList().getMoleculeCount();
            double fac = Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
            System.out.println(currentTime+" "+((double)numIterations)/numBondedMolecules+" "+Kelvin.UNIT.fromSim(2*KE/moleculeCount/6)+" "
                              +fac*KE+" "+fac*PE+" "+fac*(PE+KE));
        }
    }

    public void reset() {
        super.reset();
        /*
         * Rattle Part I
         */
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            bondConstraints.redistributeForces(molecule, agentManager);
            
            IAtomList childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            Boundary boundary = box.getBoundary();
            double[] bondLengths = bondConstraints.bondLengths;

            if (childList.getAtomCount() > moved[0].length) {
                moved = new boolean[2][childList.getAtomCount()];
            }
            for (int j=0; j<childList.getAtomCount(); j++) {
                moved[1][j] = true;
            }
            
            for (int iter = 0; iter<maxIterations; iter++) {
                boolean success = true;
                for (int j=0; j<childList.getAtomCount(); j++) {
                    moved[0][j] = moved[1][j];
                    moved[1][j] = false;
                }
                for (int j=0; j<bondConstraints.bondedAtoms.length; j++) {
                    int iAtom1 = bondedAtoms[j][0];
                    int iAtom2 = bondedAtoms[j][1];
                    if (!moved[0][iAtom1] && !moved[0][iAtom2]) {
                        continue;
                    }
                    IAtomKinetic atom1 = (IAtomKinetic)childList.getAtom(iAtom1);
                    IAtomKinetic atom2 = (IAtomKinetic)childList.getAtom(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
                    dv.Ev1Mv2(atom2.getVelocity(), atom1.getVelocity());
                    double drdotdv = dr.dot(dv);
//                    System.out.println(j+" "+drdotdv);
                    double mass1 = atom1.getType().getMass();
                    double mass2 = atom2.getType().getMass();
                    double bl2 = bondLengths[j]*bondLengths[j];
                    double g = -drdotdv / ((1.0/mass1+1.0/mass2) * bl2);
                    if (Math.abs(g) > shakeTol) {
                        dr.TE(g);
                        
                        atom2.getVelocity().PEa1Tv1( 1.0/(mass2), dr);
                        atom1.getVelocity().PEa1Tv1(-1.0/(mass1), dr);
                        
                        moved[1][iAtom1] = true;
                        moved[1][iAtom2] = true;
                        success = false;
                    }
                }
                if (success) {
                    break;
                }
                if (iter == maxIterations-1) {
//                    System.err.println("failed to converge in shake for molecule "+i);
                }
            }
        }
    }
}
