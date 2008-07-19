package etomica.integrator;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtomKinetic;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;

/**
 * Integrator implementing SHAKE algorithm.
 *
 * @author Andrew Schultz
 */
public class IntegratorVelocityVerletShake extends IntegratorMD implements AtomTypeAgentManager.AgentSource, AtomLeafAgentManager.AgentSource {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationForceSum forceSum;;
    protected final IteratorDirective allAtoms;
    protected final AtomTypeAgentManager shakeAgentManager;
    protected AtomLeafAgentManager agentManager;
    protected final IVector dr;
    protected double shakeTol;
    protected int maxIterations;
    protected boolean[][] moved;
    protected IVector[] drOld;
    protected final IVector temp;

    public IntegratorVelocityVerletShake(ISimulation sim, IPotentialMaster potentialMaster, ISpace _space) {
        this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }
    
    public IntegratorVelocityVerletShake(ISimulation sim, IPotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature, ISpace _space) {
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
        shakeAgentManager = new AtomTypeAgentManager(this, sim.getSpeciesManager(), sim.getEventManager(), true);
        setShakeTolerance(1e-6);
        setMaxIterations(20);
        moved = new boolean[2][0];
        drOld = new IVector[0];
        temp = space.makeVector();
    }
    
    public void setForceSum(PotentialCalculationForceSum pc){
        forceSum = pc;
        if(box != null){
            forceSum.setAgentManager(agentManager);
        }
        
    }
    
    public void setBox(IBox p) {
        if (box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
        }
        super.setBox(p);
        agentManager = new AtomLeafAgentManager(this,p);
        forceSum.setAgentManager(agentManager);
    }

    public void setBondConstraints(ISpecies species, int[][] bondedAtoms, double[] bondLengths) {
        shakeAgentManager.setAgent(species, new BondConstraints(bondedAtoms, bondLengths));
    }
    
    public BondConstraints getBondConstratins(ISpecies species) {
        return (BondConstraints)shakeAgentManager.getAgent(species);
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

    public void doStepInternal() {
        currentTime += timeStep;

        IAtomSet molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getAtomCount(); i++) {
            IMolecule molecule = (IMolecule)molecules.getAtom(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            
        }

        // SHAKE
        for (int i=0; i<molecules.getAtomCount(); i++) {
            IMolecule molecule = (IMolecule)molecules.getAtom(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints != null) {
                IAtomSet childList = molecule.getChildList();
                IBoundary boundary = box.getBoundary();

                if (drOld.length < bondConstraints.bondedAtoms.length) {
                    IVector[] newDrOld = new IVector[bondConstraints.bondedAtoms.length];
                    System.arraycopy(drOld, 0, newDrOld, 0, drOld.length);
                    for (int j=drOld.length; j<newDrOld.length; j++) {
                        newDrOld[j] = space.makeVector();
                    }
                    drOld = newDrOld;
                }

                for (int j=0; j<bondConstraints.bondedAtoms.length; j++) {
                    IAtomPositioned atom0 = (IAtomPositioned)childList.getAtom(bondConstraints.bondedAtoms[j][0]);
                    IAtomPositioned atom1 = (IAtomPositioned)childList.getAtom(bondConstraints.bondedAtoms[j][1]);
                    drOld[j].Ev1Mv2(atom1.getPosition(), atom0.getPosition());
                    boundary.nearestImage(drOld[j]);
                }
            }
            
            IAtomSet leafList = molecule.getChildList();
            int nLeaf = leafList.getAtomCount();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                MyAgent agent = (MyAgent)agentManager.getAgent(a);
                IVector r = a.getPosition();
                IVector v = a.getVelocity();
                if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                    System.out.println("first "+a+" r="+r+", v="+v+", f="+agent.force);
                }
                v.PEa1Tv1(0.5*timeStep*((IAtomTypeLeaf)a.getType()).rm(),agent.force);  // p += f(old)*dt/2
                temp.E(r);
                r.PEa1Tv1(timeStep,v);         // r += p*dt/m
//                System.out.println(iLeaf+" "+r);
                // temporary storage
                v.Ea1Tv1(-1.0, temp);
            }

            IAtomSet childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            IBoundary boundary = box.getBoundary();
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
                    IAtomPositioned atom1 = (IAtomPositioned)childList.getAtom(iAtom1);
                    IAtomPositioned atom2 = (IAtomPositioned)childList.getAtom(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
//                    if (i==1) System.out.println(iter+" old dr "+Math.sqrt(dr.squared())+" vs "+bondLengths[j]);
                    double dr2 = dr.squared();
                    double bl2 = bondLengths[j]*bondLengths[j];
                    double diffSq = bl2 - dr2;
                    if (Math.abs(diffSq/bl2) > shakeTol) {
                        double mass1 = ((IAtomTypeLeaf)atom1.getType()).getMass();
                        double mass2 = ((IAtomTypeLeaf)atom2.getType()).getMass();
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

        
        
        forceSum.reset();
        //Compute forces on each atom
        potential.calculate(box, allAtoms, forceSum);
        
        currentKineticEnergy = 0;
        //Finish integration step
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
//            System.out.println("shook "+iLeaf+" "+a.getPosition());
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            IVector velocity = a.getVelocity();
            // v(t+dt) = (r(t+dt) - r(t))/dt + 0.5 * f(t+dt) / m
            velocity.PE(a.getPosition());
            velocity.TE(1.0/timeStep);
            velocity.PEa1Tv1(0.5*timeStep*((IAtomTypeLeaf)a.getType()).rm(),((MyAgent)agentManager.getAgent(a)).force);  //p += f(new)*dt/2
            currentKineticEnergy += ((IAtomTypeLeaf)a.getType()).getMass() * velocity.squared();
        }
        currentKineticEnergy *= 0.5;
        
        if(isothermal) {
            doThermostat();
        }
        if (stepCount%100 == 0) {
            double PE = meterPE.getDataAsScalar();
            int moleculeCount = box.getMoleculeList().getAtomCount();
            double fac = Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
            System.out.println(currentTime+" "+Kelvin.UNIT.fromSim(2*currentKineticEnergy/moleculeCount/6)+" "
                              +fac*currentKineticEnergy+" "+fac*PE+" "+fac*(PE+currentKineticEnergy));
        }
    }

    public void reset() throws ConfigurationOverlapException{
        if(!initialized) return;
        
        super.reset();

        forceSum.reset();
        potential.calculate(box, allAtoms, forceSum);
    }

    public Class getAgentClass() {
        return MyAgent.class;
    }

    public final Object makeAgent(IAtomLeaf a) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}

    public Class<BondConstraints> getTypeAgentClass() {
        return BondConstraints.class;
    }

    public final Object makeAgent(IAtomType a) {
        return null;
    }
    
    public void releaseAgent(Object agent, IAtomType atom) {}

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
    }
}
