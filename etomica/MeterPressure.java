package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;

public class MeterPressure extends simulate.Meter
{
    //private P2 momentumP2;
    //private SpeciesMeterWalls momentumSpecies;
    private final double vScale = Constants.SCALE*Constants.SCALE*Constants.SCALE;
    private double momentumSum = 0.0;
    private double timeSum = 0.0;
    double diameter = 0.15;
    public int meterIndex=0;

    public MeterPressure()
    {
        super();
        setLabel("Pressure (bar)");
    }
    
    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
        timeSum += evt.integrator.drawTimeStep * evt.integrator.integrationInterval;
        updateStatistics(phase);}

    public double currentValue()
    {
        double flux=0.0;
        int count = 0;
        for(Species s=phaseSpace.firstSpecies(); s!=null; s=s.nextSpecies()) {
           if(s.speciesIndex == meterIndex) {
              for(Atom a=s.firstAtom(); a!=terminationAtom(); a=a.nextAtom()) {
                if(Atom.ia instanceof IntegratorHard.Agent) {
                    IntegratorHard.Agent ia = (IntegratorHard.Agent)a.ia;
                    flux = 0.5*ia.pAccumulator*Constants.SCALE/(timeSum * Constants.SCALE * Constants.DEPTH);
                    ia.pAccumulator = 0.0;
                    count++;
                }
              }  //else handle pressure for soft potential
              timeSum = 0.0;
              fluxCalc=flux/(count*Constants.BAR2SIM);   //count should be total area instead of number of atoms
              break;
           }
        }
        
        return flux;
    }
    
    public int getMeterIndex() {return meterIndex;}
    
    public void setMeterIndex(int index) {meterIndex = index;}
}