package etomica.graphics;
import etomica.*;

// Writes current values of all meters to console

public class DisplayToConsole extends Display implements Meter.MultiUser, EtomicaElement
{
    public String getVersion() {return "DisplayToConsole:01.02.10.0/"+Display.VERSION;}
    Meter[] meter;
    int nMeters = 0;
        
    public DisplayToConsole() {
        this(Simulation.instance);
    }
    public DisplayToConsole(Simulation sim) 
    {
        super(sim);
        if(meter == null) meter = new Meter[0];
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pipes meter data to console");
        return info;
    }
        
    public void doUpdate() {
        for(int i=0; i<meter.length; i++) {
            System.out.println(meter[i].getLabel() + " " + meter[i].currentValue());
        }
        System.out.println();
    }
    
    public void setMeters(Meter[] m) {
        meter = m;
    }
    public Meter[] getMeters() {return meter;}
    
    public void addMeter(Meter m) {
        if(meter == null) meter = new Meter[0];
        nMeters++;
        Meter[] temp = new Meter[nMeters];
        for(int i=0; i<meter.length; i++) {temp[i] = meter[i];}
        temp[nMeters-1] = m;
        meter = temp;
    }
            
    public void setPhase(Phase p) {
        super.setPhase(p);
//        Simulation.meterManager.addObserver(new Observer() {  //need to re-examine this
//            public void update(Observable mgr, Object obj) {
//                if(obj instanceof Meter) DisplayToConsole.this.addMeter((Meter)obj);  //don't want to add instance of MeterFunction
//            }
//        });
    }
        
    public void doPaint(java.awt.Graphics g) {
    }
}