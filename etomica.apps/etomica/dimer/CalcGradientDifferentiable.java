package etomica.dimer;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.FunctionMultiDimensionalDifferentiable;
import etomica.util.numerical.FiniteDifferenceDerivative;

/**
 * 
 * 
 * @author msellers
 *
 */


public class CalcGradientDifferentiable implements FunctionMultiDimensionalDifferentiable, AgentSource {

    public Box box;
    public PotentialMaster potentialMaster;
    public double forceConstant;
    int  derivativeOrder;
    FiniteDifferenceDerivative finiteDifferenceDerivative;
    public IteratorDirective allAtoms;
    PotentialCalculationForceSum force;
    AtomAgentManager atomAgent;
    int gradDcomponent, startAtom, stopAtom;
    AtomSet movableSet;
    private final Space space;
    
    
    public CalcGradientDifferentiable(Box aBox, PotentialMaster aPotentialMaster, AtomSet movableSet, Space _space){
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
