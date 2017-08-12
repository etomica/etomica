/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.IOException;
import java.util.Formatter;

public class IntegratorEnergyMap extends IntegratorBox implements AgentSource{

    IAtom adatom;
    public MeterPotentialEnergy energy;
    String fileTail;
    private final Space space;

    public IntegratorEnergyMap(Simulation aSim, PotentialMaster potentialMaster,
                               IAtom aAdatom, String aFileTail,
                               Space _space) {
        super(potentialMaster, 1.0);
        this.fileTail = aFileTail;
        this.adatom = aAdatom;
        this.space = _space;
    }

    protected void doStepInternal() {
        try {
           
            Formatter formatter = new Formatter("energy-"+fileTail);
            Vector pos = adatom.getPosition();
            // Move atom along Y-axis, steps by 0.1
            for(int i=0; i<292; i++){ //292
                
                // Return atom to original Z position
                adatom.getPosition().setX(2, -1.6);
                
                // Move atom along Z-axis, steps by 0.1
                for(int j=0; j<213; j++){  //213
                    // --PRINT-- 
                    formatter.format("%f %7.2f %7.2f %7.2f \n",new Object[] {energy.getDataAsScalar(),pos.getX(0), pos.getX(1), pos.getX(2)});
                    
                    // Step atom by 0.1 along Z-axis
                    adatom.getPosition().setX(2, adatom.getPosition().getX(2) +0.02);
                }
                // Step atom by 0.1 along Y-axis
                adatom.getPosition().setX(1, adatom.getPosition().getX(1) + 0.02);
     
            }
            formatter.close();
        }
        catch (IOException e) {
            
        }
    }
    
    
    protected void setup() {
        super.setup();
    
        // Set variables for energy
        energy = new MeterPotentialEnergy(potentialMaster);
        energy.setBox(box);
        
        
    }

    
    
    public Class getAgentClass() {
        return IntegratorVelocityVerlet.MyAgent.class;
    }

    public Object makeAgent(IAtom a, Box agentBox) {
        return new IntegratorVelocityVerlet.MyAgent(space);
    }

    public void releaseAgent(Object agent, IAtom atom, Box agentBox) {
        // TODO Auto-generated method stub  
    }
    
}
