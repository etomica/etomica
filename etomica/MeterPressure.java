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
    
//    public void initialize() {
//        momentumP2 = new P2();
//        momentumP2.setSpecies1Index(0);
//        momentumP2.setSpecies2Index(1);
//        phase.add(momentumP2);
//        momentumSpecies = new SpeciesMeterWalls();
//        momentumSpecies.setSpeciesIndex(1);
//        phase.add(momentumSpecies);
//        phase.nMoleculeTotal -= 2;
//        phase.nAtomTotal -= 2;
//    }

    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
        timeSum += evt.integrator.drawTimeStep * evt.integrator.integrationInterval;
        updateStatistics(phase);}

    public double currentValue()
    {
        double fluxCalc=0.0;
        for(Species s=phase.firstSpecies(); s!=null; s=s.getNextSpecies()) {
           if(s.speciesIndex == meterIndex) {
 //               s = (SpeciesWalls)s;
              double flux = 0.5*((AtomHardWall)s.firstAtom()).pAccumulator*Constants.SCALE/(timeSum * Constants.SCALE * Constants.DEPTH);
              ((AtomHardWall)s.firstAtom()).pAccumulator = 0.0;
              timeSum = 0.0;
              fluxCalc=flux/Constants.BAR2SIM;
           }
        }
        
        return fluxCalc;
    }
    
    public int getMeterIndex() {return meterIndex;}
    
    public void setMeterIndex(int index) {meterIndex = index;}
}