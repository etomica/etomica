package etomica;
import java.awt.Color;

    /**
     * Performs explicit, rectangle rule integration of all configurations of atoms in a phase
     * Assumes a 2-dimensional Space
     */

public class IntegratorRectangleRule extends Integrator {
    
    private double[] xMin, xMax;  //lower and upper limit of integration
    private double[] deltaX;
    private double[][] x;         //x values used in quadrature  [no of points][dimension]
    private int nPointsD;       //number of points for each dimension of integration
    private Space2D.Vector[] positions;
    private int D;              //total dimensions of integration
    private int iieCount;
    private int totalCount = 0;
    private int totalSteps;
    private boolean delay = false;  //slows down integration to better visualize
    public double current = -1.0;

    public IntegratorRectangleRule() {
        this(Simulation.instance);
    }
    public IntegratorRectangleRule(Simulation sim) {
        super(sim);
        nPointsD = 3;
        xMin = new double[2];
        xMax = new double[2];
        deltaX = new double[2];
    }
    
    public void initialize() {
        reset();
    }
    protected void doReset() {
        int nAtoms = firstPhase.atomCount();
        positions = new Space2D.Vector[nAtoms];
        D = 2*nAtoms;
        setN(nPointsD); 
        int delC = 0;
        if(nAtoms > 1) delC = 255/(nAtoms-1);
        int Rvalue = 255;
//        int i=nAtoms-1;
        int i=0;
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
//            a.setColor(Constants.RandomColor());
            a.setColor(new Color(Rvalue,0,255-Rvalue));
            Rvalue -= delC;
            positions[i] = (Space2D.Vector)a.position();
            positions[i].x = x[0][0];
            positions[i].y = x[0][1];
//            i--;
            i++;
        }
        totalCount = 0;
        totalSteps = (int)Math.pow(nPointsD,2*nAtoms);
    }
            
    public void doStep() {}
    
    public void run() {
        iieCount = interval;
        rectangleRecurse(0);
        fireIntervalEvent(intervalEvent);
        parentController.reset();
        runner = new Thread(this);
    }
    
    private void increment() {
        totalCount++;
        
        while(pauseRequested) doWait();
        if(--iieCount == 0) {
            fireIntervalEvent(intervalEvent);
            iieCount = interval;
        }
        if(delay) {
            try { Thread.sleep(300); }
            catch (InterruptedException e) { }
            fireIntervalEvent(intervalEvent);
        }
        Thread.yield();
    }
    
    public void setDelayed(boolean d) {delay = d;}
    public boolean isDelayed() {return delay;}
        

    private double rectangleRecurse(int k) {  //rectangle-rule quadrature formula
        double delX;
        double sum = 0.0;
        boolean doX;
        Space2D.Vector r;
        if(k % 2 == 0) {  //k is even
            r = positions[k >> 1];  //shift bits to divide by 2
            doX = true;
            delX = deltaX[0];
        }
        else {            //k is odd
            r = positions[(k-1) >> 1];
            doX = false;
            delX = deltaX[1];
        }
        for(int i=0; i<nPointsD; i++) {
            if(doX) r.x = x[i][0];
            else    r.y = x[i][1];
            if(k == D-1) { 
                current = Math.exp(-firstPhase.energy.potential()/temperature);
                sum += current;
                increment();  //keep track of number of steps taken, to fire interval event
            }
            else {
                sum += rectangleRecurse(k+1);
            }
        }
        return sum*delX;
        
    }
    
    public int getN() {return nPointsD;}
    public void setN(int n) {
        boolean reInitialize = (nPointsD != n);
        nPointsD = n;
        if(reInitialize) this.initialize();
        xMax[0] = firstPhase.dimensions().component(0);
        xMax[1] = firstPhase.dimensions().component(1);
        deltaX[0] = xMax[0]/(double)n;
        deltaX[1] = xMax[1]/(double)n;
        x = new double[2*firstPhase.atomCount()][2];
        for(int i=0; i<n; i++) {
            x[i][0] = (i+0.5)*deltaX[0];
            x[i][1] = (i+0.5)*deltaX[1];
        }
    }
    
    public int percentComplete() {
        return (int)(100*(double)totalCount/(double)totalSteps);
    }
    
    public Integrator.Agent makeAgent(Atom a) {return null;}

}