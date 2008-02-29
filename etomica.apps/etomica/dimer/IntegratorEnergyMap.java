package etomica.dimer;

import java.io.IOException;
import java.util.Formatter;

import etomica.api.IVector;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.Space;

public class IntegratorEnergyMap extends IntegratorBox implements AgentSource{

    IAtomPositioned adatom;
    public MeterPotentialEnergy energy;
    String fileTail;
    private final Space space;

    public IntegratorEnergyMap(ISimulation aSim, PotentialMaster potentialMaster,
    		                   IAtomPositioned aAdatom, String aFileTail,
    		                   Space _space) {
        super(potentialMaster, 1.0);
        this.fileTail = aFileTail;
        this.adatom = aAdatom;
        this.space = _space;
    }
    
    public void doStepInternal(){
        try {
           
            Formatter formatter = new Formatter("energy-"+fileTail);
            IVector pos = adatom.getPosition();
            // Move atom along Y-axis, steps by 0.1
            for(int i=0; i<292; i++){ //292
                
                // Return atom to original Z position
                adatom.getPosition().setX(2, -1.6);
                
                // Move atom along Z-axis, steps by 0.1
                for(int j=0; j<213; j++){  //213
                    // --PRINT-- 
                    formatter.format("%f %7.2f %7.2f %7.2f \n",new Object[] {energy.getDataAsScalar(),pos.x(0), pos.x(1), pos.x(2)});
                    
                    // Step atom by 0.1 along Z-axis
                    adatom.getPosition().setX(2, adatom.getPosition().x(2) +0.02);
                }
                // Step atom by 0.1 along Y-axis
                adatom.getPosition().setX(1, adatom.getPosition().x(1) + 0.02);
     
            }
            formatter.close();
        }
        catch (IOException e) {
            
        }
    }
    
    
    protected void setup() throws ConfigurationOverlapException{
        super.setup();
    
        // Set variables for energy
        energy = new MeterPotentialEnergy(this.potential);
        energy.setBox(box);
        
        
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
