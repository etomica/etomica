/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.math.function.FunctionMultiDimensionalDifferentiable;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Uses finite difference methods to determine the second order differential of the potential (i.e. dF/dx).
 * Part of a larger scheme, a user may employ this to fill a two-dimensional array of changes in molecule 
 * A's Force with respect to the movement of other molecules.
 * 
 * @author msellers and ajschultz
 *
 */


public class CalcGradientDifferentiable implements FunctionMultiDimensionalDifferentiable {

    public Box box;
    public PotentialMaster potentialMaster;
    public double forceConstant;
    int  derivativeOrder;
    FiniteDifferenceDerivative finiteDifferenceDerivative;
    public IteratorDirective allAtoms;
    PotentialCalculationForceSum force;
    AtomLeafAgentManager<Vector> atomAgent;
    int gradDcomponent, startAtom, stopAtom;
    IMoleculeList movableSet;
    private final Space space;
    
    
    public CalcGradientDifferentiable(Box aBox, PotentialMaster aPotentialMaster, IMoleculeList movableSet, Space _space){
        this.box = aBox;
        this.potentialMaster = aPotentialMaster;
        this.movableSet = movableSet;
        this.space = _space;
        
        if(potentialMaster instanceof PotentialMasterList){
            ((PotentialMasterList)potentialMaster).getNeighborManager(box).reset();
         }
        
        force = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        atomAgent = new AtomLeafAgentManager<>(a -> space.makeVector(), box);
        force.setAgentManager(atomAgent);
        
        finiteDifferenceDerivative = new FiniteDifferenceDerivative(this);
        finiteDifferenceDerivative.setH(0.00001);
    }
    

    public void setComponent(int aGradDcomponent){
        this.gradDcomponent = aGradDcomponent;
    }
    
    
    public double f(double [] position){
        
        for(int i=0; i<position.length/3; i++){
           for(int j=0; j<3; j++){
        	   movableSet.getMolecule(i).getChildList().get(0).getPosition().setX(j, position[(3*i)+j]);
           }
        }
        force.reset();
        potentialMaster.calculate(box, allAtoms, force);
        
        //not used
        return atomAgent.getAgent(movableSet.getMolecule(gradDcomponent/3).getChildList().get(0)).getX(gradDcomponent%3);

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
            forceRow[j] = -atomAgent.getAgent(movableSet.getMolecule(j/3).getChildList().get(0)).getX(j%3);
        }
        
        position[elem]-=2*newH;
        f(position);
        
        for(int j=0; j<forceRow.length; j++){
            forceRow[j] -= -atomAgent.getAgent(movableSet.getMolecule(j/3).getChildList().get(0)).getX(j%3);
            forceRow[j] /= (2.0*newH);
        }
        
        return forceRow;
    }
    
    public int getDimension(){
        return space.D();
    }
}
