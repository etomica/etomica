package etomica.integrator;

import etomica.action.AtomAction;
import etomica.api.IAtomSet;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.exception.ConfigurationOverlapException;
import etomica.space.Space;

/**
 * Integrator that generates atom trajectories from an analytic formula.
 * Takes an IntegratorAnalytic.AtomAction instance and performs action on all atoms in each
 * time step; intention is for the action to set the position of the atom for the
 * current time (but it could do anything).
 *
 * @author David Kofke
 * @author Nancy Cribbin
 */
public class IntegratorAnalytic extends IntegratorMD {
    
    private static final long serialVersionUID = 1L;
    private AtomTimeAction action;
    
    public IntegratorAnalytic(ISimulation sim, IPotentialMaster potentialMaster, Space _space) {
        this(potentialMaster, sim.getRandom(), 0.05, _space);
    }
    
    public IntegratorAnalytic(IPotentialMaster potentialMaster, IRandom random,
                              double timeStep, Space _space) {
        super(potentialMaster,random,timeStep,0, _space);
    }
    
    public void doStepInternal() {
        super.doStepInternal();
        if(action == null) return;
        elapsedTime += getTimeStep();
        action.setTime(elapsedTime);
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            action.actionPerformed(leafList.getAtom(iLeaf));
        }
    }
    
    public void reset() throws ConfigurationOverlapException {
        elapsedTime = 0.0;
        super.reset();
    }
    
    public void setAction(AtomTimeAction action) {this.action = action;}
    
    public AtomTimeAction getAction() {return action;}
    
    private double elapsedTime = 0.0;
    
    /**
     * Extends AtomAction class to add a method to set the time.
     */
    public static abstract class AtomTimeAction implements AtomAction {
        protected double time;
        public void setTime(double t) {time = t;}
    }
    
 }//end of IntegratorAnalytic
 
