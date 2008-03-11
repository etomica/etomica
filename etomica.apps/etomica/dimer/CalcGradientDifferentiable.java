package etomica.dimer;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.FunctionMultiDimensionalDifferentiable;
import etomica.util.numerical.FiniteDifferenceDerivative;

/**
 * Uses finite difference methods to determine the second order differential of the potential (i.e. dF/dx).
 * Part of a larger scheme, a user may employ this to fill a two-dimensional array of changes in molecule 
 * A's Force with respect to the movement of other molecules.
 * 
 * @author msellers and ajschultz
 *
 */


public class CalcGradientDifferentiable implements FunctionMultiDimensionalDifferentiable, AgentSource {

    public IBox box;
    public PotentialMaster potentialMaster;
    public double forceConstant;
    int  derivativeOrder;
    FiniteDifferenceDerivative finiteDifferenceDerivative;
    public IteratorDirective allAtoms;
    PotentialCalculationForceSum force;
    AtomAgentManager atomAgent;
    int gradDcomponent, startAtom, stopAtom;
    IAtomSet movableSet;
    private final Space space;
    
    
    public CalcGradientDifferentiable(IBox aBox, PotentialMaster aPotentialMaster, IAtomSet movableSet, Space _space){
        this.box = aBox;
        this.potentialMaster = aPotentialMaster;
        this.movableSet = movableSet;
        this.space = _space;
       
        
        force = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        atomAgent = new AtomAgentManager(this, box);
        force.setAgentManager(atomAgent);
        
        finiteDifferenceDerivative = new FiniteDifferenceDerivative(this);
        finiteDifferenceDerivative.setH(0.00001);
        finiteDifferenceDerivative.setHOptimizer(true);
        finiteDifferenceDerivative.setNtab(10);
    }
    

    public void setComponent(int aGradDcomponent){
        this.gradDcomponent = aGradDcomponent;
    }
    
    
    public double f(double [] position){
        
        for(int i=0; i<position.length/3; i++){
           for(int j=0; j<3; j++){
        	   ((IAtomPositioned)((IMolecule)movableSet.getAtom(i)).getChildList().getAtom(0)).getPosition().setX(j, position[(3*i)+j]);
           }
        }
        force.reset();
        potentialMaster.calculate(box, allAtoms, force);
        
        //not used
        return ((IntegratorVelocityVerlet.MyAgent)atomAgent.getAgent(((IMolecule)movableSet.getAtom(gradDcomponent/3)).getChildList().getAtom(0))).force().x(gradDcomponent%3);

    }
    
    public double df(int [] d, double [] position){
        return finiteDifferenceDerivative.df(d, position);
    }
    
    /**
     * Uses the potential's force calculation at different displacements of a molecule in X, Y and Z
     * to determine the second derivative of the potential.  Here, H (our displacement distance) is
     * equal to 0.0001.
     * 
     * @param d A one-dimensional array describing what column of our larger, global dF/dx array we are working with.
     * @param position A one dimensional array of doubles describing the molecules positions.
     * @return
     */
    public double[] df2(int [] d, double [] position){
        double newH = 0.0001;
        double[] forceRow = new double[d.length];
        int elem = 0;
        
        for(int i=0; i<d.length; i++){
            if(d[i]==1){
                elem = i;
                break;                
            }
        }
        
        position[elem]+=newH;
        f(position);
        
        for(int j=0; j<forceRow.length; j++){
            forceRow[j] = -((IntegratorVelocityVerlet.MyAgent)atomAgent.getAgent(((IMolecule)movableSet.getAtom(j/3)).getChildList().getAtom(0))).force.x(j%3);
        }
        
        position[elem]-=2*newH;
        f(position);
        
        for(int j=0; j<forceRow.length; j++){
            forceRow[j] -= -((IntegratorVelocityVerlet.MyAgent)atomAgent.getAgent(((IMolecule)movableSet.getAtom(j/3)).getChildList().getAtom(0))).force().x(j%3);
            forceRow[j] /= (2.0*newH);
        }
        
        return forceRow;
    }
    
    public int getDimension(){
        return space.D();
    }

    
    public Class getAgentClass() {
        return IntegratorVelocityVerlet.MyAgent.class;
    }

    public Object makeAgent(IAtom a) {
        return new IntegratorVelocityVerlet.MyAgent(space);
    }

    public void releaseAgent(Object agent, IAtom atom) {
        // TODO Auto-generated method stub  
    }
    
}
