package etomica.veos;

public class VanderWaals implements BetaSource {

    public VanderWaals(double temperature) {
        this(temperature, Integer.MAX_VALUE);
    }
    
    public VanderWaals(double temperature, int m) {
        maxIndex = m;
        T = temperature;
    }
    
    public void setTemperature(double t) {
        T = t;
    }
    
    public double beta(int k) {
        if(k <= 0) throw new IllegalArgumentException("invalid k in VanderWaals.beta");
        if(k == 1) {
            return -2*(1 - 1/T);
        } else if(k <= maxIndex) {
            return -(double)(k+1)/k;
        } else {
            return 0.0;
        }
    }
    
    public double kbeta(int k) {
        if(k <= 0) throw new IllegalArgumentException("invalid k in VanderWaals.beta");
        if(k == 1) {
            return -2*(1 - 1/T);
        } else if(k <= maxIndex) {
            return -(double)(k+1);
        } else {
            return 0.0;
        }
    }


    public void setMaxIndex(int m) {
        maxIndex = m;
    }

    public int maxIndex() {
        return maxIndex;
    }

    private int maxIndex;
    private double T;
}
