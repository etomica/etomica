package simulate;

import simulate.*;
import java.beans.*;

//Number of integrationIntervalEvent's fired

public class MeterCycles extends simulate.Meter
{
//    private int count = 0;
    
    public MeterCycles()
    {
        super();
        setLabel("Cycles");
    }

    public double currentValue()
    {
 //       count++;        
        return (double)count;  
        //count is defined in superclass
        //it is updated every time integrationIntervalEvent is received
        //  (assuming updateInterval is kept at 1)
    }
}