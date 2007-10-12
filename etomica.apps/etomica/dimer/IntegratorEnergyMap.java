package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;

public class IntegratorEnergyMap extends IntegratorBox implements AgentSource{

    IAtom adatom;
    public MeterPotentialEnergy energy;
    
    public IntegratorEnergyMap(ISimulation aSim, PotentialMaster potentialMaster, IAtom aAdatom) {
        super(potentialMaster, 1.0);
    
        adatom = aAdatom;
    
        
    
    
    }
    
    public void doStepInternal(){
        try {
        
        FileWriter writer = new FileWriter("energy-cu-63-45cut");
        
        // Move atom along Y-axis, steps by 0.1
        for(int i=0; i<500; i++){
            
            // Move atom along Z-axis, steps by 0.1
            for(int j=0; j<500; j++){
                
                // Step atom by 0.1 along Z-axis
                ((IAtomPositioned)adatom).getPosition().setX(2, ((IAtomPositioned)adatom).getPosition().x(2) +0.01);
                
                // --PRINT-- 
                writer.write(energy.getDataAsScalar()+"     "+((IAtomPositioned)adatom).getPosition()+"\n");
            }
            
            // Return atom to original Z position
            ((IAtomPositioned)adatom).getPosition().setX(2, -1.2);
            
            // Step atom by 0.1 along Y-axis
            ((IAtomPositioned)adatom).getPosition().setX(1, ((IAtomPositioned)adatom).getPosition().x(1) + 0.01);
            
            // --PRINT--
            writer.write(energy.getDataAsScalar()+"     "+((IAtomPositioned)adatom).getPosition()+"\n");
            
        }
        
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
