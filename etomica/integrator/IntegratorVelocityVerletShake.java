package etomica.integrator;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBoundary;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomTypeAgentManager;
import etomica.space.Space;

/**
 * Integrator implementing SHAKE algorithm.  Use adiabatically at your own risk.
 *
 * @author Andrew Schultz
 */
public class IntegratorVelocityVerletShake extends IntegratorVelocityVerlet implements AtomTypeAgentManager.AgentSource {

    private static final long serialVersionUID = 2L;
    protected final AtomTypeAgentManager shakeAgentManager;
    protected final IVector dr;
    protected double shakeTol2;
    protected int maxIterations;
    protected boolean[][] moved;

    public IntegratorVelocityVerletShake(ISimulation sim, IPotentialMaster potentialMaster, Space _space) {
        this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }
    
    public IntegratorVelocityVerletShake(ISimulation sim, IPotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature, Space _space) {
        super(potentialMaster,random,timeStep,temperature, _space);
        dr = _space.makeVector();
        shakeAgentManager = new AtomTypeAgentManager(this, sim.getSpeciesManager(), sim.getEventManager(), true);
        setShakeTolerance(1e-10);
        setMaxIterations(10);
        moved = new boolean[2][0];
    }
    
    public void setBondConstraints(ISpecies species, int[][] bondedAtoms, double[] bondLengths) {
        shakeAgentManager.setAgent(species, new BondConstraints(bondedAtoms, bondLengths));
    }
    
    public BondConstraints getBondConstratins(ISpecies species) {
        return (BondConstraints)shakeAgentManager.getAgent(species);
    }
    
    public void setShakeTolerance(double newShakeTol) {
        shakeTol2 = newShakeTol*newShakeTol;
    }
    
    public double getShakeTolerance() {
        return Math.sqrt(shakeTol2);
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public void setMaxIterations(int newMaxIterations) {
        maxIterations = newMaxIterations;
    }

    public void doStepInternal() {
        super.doStepInternal();
        IAtomSet molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getAtomCount(); i++) {
            IMolecule molecule = (IMolecule)molecules.getAtom(i);
            BondConstraints bondConstraints = (BondConstraints)shakeAgentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            
            IAtomSet childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            IBoundary boundary = box.getBoundary();
            double[] bondLengths = bondConstraints.bondLengths;

            if (childList.getAtomCount() > moved.length) {
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
                    double dr2 = dr.squared();
                    double bl2 = bondLengths[j]*bondLengths[j];
                    double diffSq = dr2 - bl2;
                    if (Math.abs(diffSq/bl2) > shakeTol2) {
                        double mass1 = ((IAtomTypeLeaf)atom1.getType()).getMass();
                        double mass2 = ((IAtomTypeLeaf)atom2.getType()).getMass();
                        double rMass = 1.0/mass1 + 1.0/mass2;
                        double R = Math.sqrt(dr2);
                        double fac = (1.0 - bondLengths[j]/R)/rMass;
                        atom2.getPosition().PEa1Tv1(-fac/mass2, dr);
                        atom1.getPosition().PEa1Tv1(fac/mass1, dr);
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
//        if (stepCount%10 == 0) {
//            double PE = meterPE.getDataAsScalar();
//            double KE = meterKE.getDataAsScalar();
//            int moleculeCount = box.getMoleculeList().getAtomCount();
//            double fac = Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
//            System.out.println(currentTime+" "+Kelvin.UNIT.fromSim(2*KE/moleculeCount/9)+" "
//                              +fac*KE+" "+fac*PE+" "+fac*(PE+KE));
//        }
    }

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
