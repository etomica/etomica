package simulate.electrostatics;

import simulate.*;
import utility.Functions;

// Ewald sum; contains some elements particular to a 2D system

public class EwaldSum extends MeterMultiFunction {
    
    AtomPair.Iterator.All iteratorAll;
    Atom.Iterator atomIterator;
    double[] xDirect, xRecip, yDirect, yRecip;
    int nDirect, nRecip;
    private double z2Sum;
    private static final double selfC = -1.0/Math.sqrt(Math.PI);
    private static final double twoPi = 2.0*Math.PI;
    private double alpha;
    private double m2Max;
    public static final int TOTAL = 0;
    public static final int DIRECT = 1;
    public static final int RECIP = 2;
    public static final int SELF = 3;
    
    public EwaldSum() {
        super(4);
        setNDirect(5);
        setNRecip(5);
//        setAlpha(5.0);
        meters[TOTAL] = new TotalMeter();
        meters[DIRECT] = new DirectMeter();
        meters[RECIP] = new RecipMeter();
        meters[SELF] = new SelfMeter();
    }
    
    public double[] currentValue() {return meters[TOTAL].currentValue();}
        
      
    public void setNDirect(int n) {
        nDirect = n;
        xDirect = new double[n];
        yDirect = new double[n];
        for(int i=0; i<nDirect; i++) {xDirect[i] = i;}
        meters[DIRECT] = new DirectMeter();
        meters[TOTAL] = new TotalMeter();
    }
    public void setNRecip(int n) {
        nRecip = n;
        xRecip = new double[n];
        yRecip = new double[n];
        for(int i=0; i<nRecip; i++) {xRecip[i] = i;}
        meters[RECIP] = new RecipMeter();
        meters[TOTAL] = new TotalMeter();
    }
    
    public void setAlpha(double a) {if(phase!=null) alpha = a/phase.boundary().dimensions().component(0);}
    public double getAlpha() {return alpha;}
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        iteratorAll = p.iterator.makeAtomPairIteratorAll();
        atomIterator = p.iterator.makeAtomIteratorUp();
        setAlpha(5.);
    }
    
    public double self() {
        atomIterator.reset();
        z2Sum = 0.0;
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            double z = ((Monopole)a.type.electroType()).z();
            z2Sum += z*z;
        }
        return selfC*alpha*z2Sum;
    }
    
    private double recip(Space2D.Vector M) {
        double m2 = M.squared();
        if(m2 > m2Max) return 0.0;
        atomIterator.reset();
        double cSum = 0.0;
        double sSum = 0.0;
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            double z = ((Monopole)a.type.electroType()).z();
            double arg = M.dot(a.r);
            cSum += z*Math.cos(arg);
            sSum += z*Math.sin(arg);
        }
        return twoPi/phase.volume()*Math.exp(-0.25*m2/(alpha*alpha))/m2 * (cSum*cSum + sSum*sSum);
    }
    
    double[] recipSum() { 
        Space2D.Vector M = (Space2D.Vector)phase.parentSimulation().space.makeVector();
        double L = phase.boundary().dimensions().component(0);
        double c = twoPi/L;
        m2Max = 1.001*((nRecip-1)*(nRecip-1)*c*c);
        double peSum = 0.0;
        yRecip[0] = peSum;
        for(int m=1; m<nRecip; m++) {
            M.x = c*m;
            M.y = c*m;
            peSum += recip(M);
            M.x = -c*m;
            peSum += recip(M);
            M.y = -c*m;
            peSum += recip(M);
            M.x = +c*m;
            peSum += recip(M);
            for(int i=-m+1; i<m; i++) {
                M.x = +c*m;
                M.y = +c*i;
                peSum += recip(M);
                M.x = -c*m;
                peSum += recip(M);
                M.x = +c*i;
                M.y = +c*m;
                peSum += recip(M);
                M.y = -c*m;
                peSum += recip(M);
            }
            yRecip[m] = peSum;
        }
        return yRecip;
     }
     
    /**
     * Performs a sum of the energy over periodic images such that each additional
     * shell yields a cubic shaped region
     * Assumes a 2-D cubic central shell (same length and width)
     */
    double[] directSum() { 
        Space2D.Vector M = (Space2D.Vector)phase.parentSimulation().space.makeVector();
        double L = phase.boundary().dimensions().component(0);
        m2Max = 1.001*((nDirect-1)*(nDirect-1)*L*L);
        double peSum = direct(M);
        yDirect[0] = peSum;
        for(int m=1; m<nDirect; m++) {
            M.x = m*L;
            M.y = m*L;
            peSum += direct(M);
            M.x = -m*L;
            peSum += direct(M);
            M.y = -m*L;
            peSum += direct(M);
            M.x = +m*L;
            peSum += direct(M);
            for(int i=-m+1; i<m; i++) {
                M.x = +m*L;
                M.y = +i*L;
                peSum += direct(M);
                M.x = -m*L;
                peSum += direct(M);
                M.x = +i*L;
                M.y = +m*L;
                peSum += direct(M);
                M.y = -m*L;
                peSum += direct(M);
            }
            yDirect[m] = peSum;
        }
        return yDirect;
     }
     
     private double direct(Space2D.Vector M) {
        double m2 = M.squared();
        if(m2 > m2Max) return 0.0;
        double pe = 0.0;
        iteratorAll.reset();
        while(iteratorAll.hasNext()) {
            AtomPair pair = iteratorAll.next();
            if(m2 > -10.0) ((Space2D.CoordinatePair)pair.cPair).reset(M);
            double r = Math.sqrt(pair.r2());
            double z1 = ((Monopole)pair.atom1.type.electroType()).z();
            double z2 = ((Monopole)pair.atom2.type.electroType()).z();
            pe += z1*z2*Functions.erfc(alpha*r)/r;  
        }
        double m = Math.sqrt(m2);
        double self = (m2 > 0.0) ? z2Sum*Functions.erfc(alpha*m)/m : 0.0;  //self-image interaction  //erfc goes here too
        return 2.0*pe + self;
     }
     
     private class DirectMeter extends MeterFunction {
        DirectMeter() {
            setNPoints(nDirect);
            setLabel("Direct");
        }
        public double[] currentValue() {return yDirect;}
        public void setX(double xmin, double xmax, int n) {
            super.setX(xmin,xmax,n);
            for(int i=0; i<n; i++) {x[i] = i;}
        }
     }
     
     private class RecipMeter extends MeterFunction {
        RecipMeter() {
            setNPoints(nRecip);
            setLabel("Reciprocal");
        }
        public double[] currentValue() {return yRecip;}
        public void setX(double xmin, double xmax, int n) {
            super.setX(xmin,xmax,n);
            for(int i=0; i<n; i++) {x[i] = i;}
        }
     }
     
     private class SelfMeter extends MeterFunction {
        SelfMeter() {
            setNPoints(2);
            setLabel("Self");
        }
        public double[] currentValue() {
           y[0] = y[1] = self();
           return y;
        }
        public void setX(double xmin, double xmax, int n) {
            super.setX(xmin,xmax,n);
            for(int i=0; i<n; i++) {x[i] = i;}
        }
     }
     
     private class TotalMeter extends MeterFunction {
        TotalMeter() {
            setNPoints(Math.max(nRecip, nDirect));
            setLabel("Total");
        }
        public double[] currentValue() {
            double s = self();
            double[] r = recipSum();
            double[] d = directSum();
            for(int i=0; i<nPoints; i++) {
                y[i] = s + r[Math.min(i,nRecip-1)] + d[Math.min(i,nDirect-1)];
            }
            return y;
        }
        public void setX(double xmin, double xmax, int n) {
            super.setX(xmin,xmax,n);
            for(int i=0; i<n; i++) {x[i] = i;}
        }
     }
        
}