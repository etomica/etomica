/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Apr 12, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.adsorption;

import etomica.atom.IAtomList;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.nbr.PotentialMasterHybrid;



/**
 * @author ecc4
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
public class IntegratorHybrid extends IntegratorBox {
	
	protected final IntegratorMC integratorMC;
	protected final IntegratorMD integratorMD;
	private MyMCMove mcMoveIDA, mcMoveIDB;
    private final PotentialMasterHybrid potentialMasterHybrid;
	private int MDStepCount, MDStepRepetitions;
    
	public IntegratorHybrid(PotentialMasterHybrid potentialMaster, IntegratorMD integratorMD, IntegratorMC integratorMC, double temperature) {
		super(potentialMaster, temperature, integratorMD.getBox());
		this.integratorMD = integratorMD;
		this.integratorMC = integratorMC;
		potentialMasterHybrid = potentialMaster;
        setMDStepRepetitions(50);
    }
	
	public void setMCMoveInsertDelete(MyMCMove mcMoveIDA, MyMCMove mcMoveIDB) {
	    this.mcMoveIDA = mcMoveIDA;
        this.mcMoveIDB = mcMoveIDB;
	}
    
    public void setMDStepRepetitions(int interval) {
        MDStepRepetitions = interval;
        if (MDStepCount > interval || MDStepCount == 0) MDStepCount = interval;
    }
    
    protected void setup() {
        super.setup();
        potentialMasterHybrid.setUseNbrLists(false);
        integratorMC.reset();
        potentialMasterHybrid.setUseNbrLists(true);
        integratorMD.reset();
    }
    
    public void setTemperature(double t) {
        super.setTemperature(t);
        if (integratorMC != null) {
            integratorMC.setTemperature(t);
        }
    }
    
    public void setIsothermal(boolean b){
    	super.setIsothermal(b);
    	integratorMD.setIsothermal(b);
    }

	protected void doStepInternal() {
		if (potentialMasterHybrid != null) {
			potentialMasterHybrid.setUseNbrLists(MDStepCount > 0);
        }
		if(MDStepCount == 0){
		    MDStepCount = MDStepRepetitions;
            mcMoveIDA.setupActiveAtoms();
			mcMoveIDB.setupActiveAtoms();
			for(int i=0; i<1; i++) {
                integratorMC.doStep();
            }
			IAtomList allAtoms = box.getLeafList();
			for (int i = 0; i<allAtoms.size(); i++) {
			    if (allAtoms.get(i).getPosition().getX(2) < -40) {
			        throw new RuntimeException(i+" "+allAtoms.get(i)+" "+allAtoms.get(i).getPosition());
			    }
			}
            potentialMasterHybrid.setUseNbrLists(true);
            potentialMasterHybrid.getNeighborManager(box).reset();
            integratorMD.reset();
	 	} else {
            MDStepCount--;
	 		integratorMD.doStep();
            IAtomList allAtoms = box.getLeafList();
            for (int i = 0; i<allAtoms.size(); i++) {
                if (allAtoms.get(i).getPosition().getX(2) < -40) {
                    throw new RuntimeException(i+" "+allAtoms.get(i)+" "+allAtoms.get(i).getPosition());
                }
            }
		} 
	}

	public double getCurrentTime() {
		return integratorMD.getCurrentTime();
	}
	
	public void reset() {
	    super.reset();
        potentialMasterHybrid.setUseNbrLists(false);
		integratorMC.reset();
        potentialMasterHybrid.setUseNbrLists(true);
		integratorMD.reset();
	}

	public void resetStepCount() {
	    super.resetStepCount();
	    integratorMC.resetStepCount();
	    integratorMD.resetStepCount();
	}
}
