package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;

// Writes current values of all meters to console

public class DisplayToConsole extends simulate.Display
{
    Meter[] meter = new Meter[0];
    int nMeters = 0;
    
    public DisplayToConsole()
    {
        super();
    }
    
    public void doUpdate() {
    }
    
        public void addMeter(Meter m) {
            nMeters++;
            Meter[] temp = new Meter[nMeters];
            for(int i=0; i<meter.length; i++) {temp[i] = meter[i];}
            temp[nMeters-1] = m;
            meter = temp;
        }
        
        public void setPhase(Phase p) {
            for(Meter m=p.firstMeter; m!=null; m=m.nextMeter()) {this.addMeter(m);}
        }
    
    public void doPaint(Graphics g) {
        for(int i=0; i<nMeters; i++) {
            System.out.println(meter[i].getLabel() + " " + meter[i].currentValue());
        }
        System.out.println();
    }
    
    public void paint(Graphics g) {
      createOffScreen();
      doPaint(g);
    }
}