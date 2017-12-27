/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.*;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Integrator implementing SHAKE algorithm.
 *
 * @author Andrew Schultz
 */
public class IntegratorVelocityVerletShake extends IntegratorMD implements SpeciesAgentManager.AgentSource, AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent> {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationForceSum forceSum;;
    protected final IteratorDirective allAtoms;
    protected final SpeciesAgentManager shakeAgentManager;
    protected AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent> agentManager;
    protected final Vector dr;
    protected double shakeTol;
    protected int maxIterations;
    protected boolean[][] moved;
    protected Vector[] drOld;
    protected final Vector temp;
    public int printInterval = 0;

    public IntegratorVelocityVerletShake(Simulation sim, PotentialMaster potentialMaster, Space _space) {
        this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }
    
    public IntegratorVelocityVerletShake(Simulation sim, PotentialMaster potentialMaster, IRandom random,
                                         double timeStep, double temperature, Space _space) {
        super(potentialMaster,random,timeStep,temperature, _space);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        forceSum = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);

        
        dr = _space.makeVector();
        shakeAgentManager = new SpeciesAgentManager(this, sim);
        setShakeTolerance(1e-14);
        setMaxIterations(20);
        moved = new boolean[2][0];
        drOld = new Vector[0];
        temp = space.makeVector();
    }
    
    public void setForceSum(PotentialCalculationForceSum pc){
        forceSum = pc;
        if(box != null){
            forceSum.setAgentManager(agentManager);
        }
        
    }
    
    public void setBox(Box box) {
        if (this.box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
        }
        super.setBox(box);
        agentManager = new AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent>(this, box);
        forceSum.setAgentManager(agentManager);
    }

    public void setBondConstraints(ISpecies species, int[][] bondedAtoms, double[] bondLengths) {
        shakeAgentManager.setAgent(species, new BondConstraints(bondedAtoms, bondLengths));
    }
    
    public BondConstraints getBondConstratins(ISpecies species) {
        return (BondConstraints)shakeAgentManager.getAgent(species);
    }
    
    public SpeciesAgentManager getShakeAgentManager(){
    	return shakeAgentManager;
    }
    
    public void setShakeTolerance(double newShakeTol) {
        shakeTol = newShakeTol;
    }
    
    public double getShakeTolerance() {
        return shakeTol;
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public void setMaxIterations(int newMaxIterations) {
        maxIterations = newMaxIterations;
    }

    protected void doStepInternal() {
        currentTime += timeStep;

        IMoleculeList molecules = box.getMoleculeList();

        // SHAKE
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints != null) {
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
                if (a.getType().getMass() != 0) {
                    v.PEa1Tv1(0.5*timeStep*a.getType().rm(),agent.force);  // p += f(old)*dt/2
                }
                temp.E(r);
                r.PEa1Tv1(timeStep,v);         // r += p*dt/m
//                System.out.println(iLeaf+" "+r);
                // temporary storage
                v.Ea1Tv1(-1.0, temp);
            }

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
                    IAtom atom1 = childList.getAtom(iAtom1);
                    IAtom atom2 = childList.getAtom(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
//                    if (i==0) System.out.println(iter+" old dr "+Math.sqrt(dr.squared())+" vs "+bondLengths[j]);
                    double dr2 = dr.squared();
                    double bl2 = bondLengths[j]*bondLengths[j];
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

        currentKineticEnergy = 0;
        //Finish integration step
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
//            System.out.println("shook "+iLeaf+" "+a.getPosition());
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            Vector velocity = a.getVelocity();
            // v(t+dt) = (r(t+dt) - r(t))/dt + 0.5 * f(t+dt) / m
            velocity.PE(a.getPosition());
            velocity.TE(1.0/timeStep);
            if (a.getType().getMass() != 0) {
                velocity.PEa1Tv1(0.5*timeStep*a.getType().rm(),agentManager.getAgent(a).force);  //p += f(new)*dt/2
            }
            currentKineticEnergy += a.getType().getMass() * velocity.squared();
        }
        currentKineticEnergy *= 0.5;
        
        if(isothermal) {
            doThermostatInternal();
        }
        if (printInterval > 0 && stepCount%printInterval == 0) {
            double PE = meterPE.getDataAsScalar();
            int moleculeCount = box.getMoleculeList().getMoleculeCount();
            double fac = Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
            System.out.println(currentTime+" "+Kelvin.UNIT.fromSim(2*currentKineticEnergy/moleculeCount/6)+" "
                              +fac*currentKineticEnergy+" "+fac*PE+" "+fac*(PE+currentKineticEnergy));
        }
    }

    public void reset() {
        super.reset();

        forceSum.reset();
        potentialMaster.calculate(box, allAtoms, forceSum);
    }

    public final IntegratorVelocityVerlet.MyAgent makeAgent(IAtom a, Box agentBox) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(IntegratorVelocityVerlet.MyAgent agent, IAtom atom, Box agentBox) {}

    public final Object makeAgent(ISpecies a) {
        return null;
    }
    
    public void releaseAgent(Object agent, ISpecies atom) {}

    public static class BondConstraints {
        public final int[][] bondedAtoms;
        public final double[] bondLengths;
        public BondConstraints(int[][] bondedAtoms, double[] bondLengths) {
            if (bondedAtoms.length != bondLengths.length) {
                throw new IllegalArgumentException("different number of bonded pairs and lengths");
            }
            this.bondedAtoms = bondedAtoms;
            this.bondLengths = bondLengths;
        }

        // redistribute forces to constrained atoms
        // do nothing by default, allow subclasses to override
        public void redistributeForces(IMolecule molecule, AtomLeafAgentManager<MyAgent> agentManager) {}

        // fix atom positions that may have been omitted from constraints
        // do nothing by default, allow subclasses to override
        public void relaxMolecule(IMolecule molecule) {}
    }
}
