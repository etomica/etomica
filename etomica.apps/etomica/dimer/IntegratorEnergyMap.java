package etomica.dimer;

import java.io.IOException;
import java.util.Formatter;

import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVector;

public class IntegratorEnergyMap extends IntegratorBox implements AgentSource{

    IAtom adatom;
    public MeterPotentialEnergy energy;
    String fileTail;
    
    public IntegratorEnergyMap(ISimulation aSim, PotentialMaster potentialMaster, IAtom aAdatom, String aFileTail) {
        super(potentialMaster, 1.0);
        this.fileTail = aFileTail;
        this.adatom = aAdatom;
    
        
    
    
    }
    
    public void doStepInternal(){
        try {
           
            Formatter formatter = new Formatter("energy-"+fileTail);
            IVector pos = ((IAtomPositioned)adatom).getPosition();
            // Move atom along Y-axis, steps by 0.1
            for(int i=0; i<292; i++){ //292
                
                // Return atom to original Z position
                ((IAtomPositioned)adatom).getPosition().setX(2, -1.6);
                
                // Move atom along Z-axis, steps by 0.1
                for(int j=0; j<213; j++){  //213
                    // --PRINT-- 
                    formatter.format("%f %7.2f %7.2f %7.2f \n",new Object[] {energy.getDataAsScalar(),pos.x(0), pos.x(1), pos.x(2)});
                    
                    // Step atom by 0.1 along Z-axis
                    ((IAtomPositioned)adatom).getPosition().setX(2, ((IAtomPositioned)adatom).getPosition().x(2) +0.02);
                }
                // Step atom by 0.1 along Y-axis
                ((IAtomPositioned)adatom).getPosition().setX(1, ((IAtomPositioned)adatom).getPosition().x(1) + 0.02);
     
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
        return new IntegratorVelocityVerlet.MyAgent(box.getSpace());
    }

    public void releaseAgent(Object agent, IAtom atom) {
        // TODO Auto-generated method stub  
    }
    
}
