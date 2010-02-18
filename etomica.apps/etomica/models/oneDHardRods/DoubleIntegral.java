package etomica.models.oneDHardRods;

public class DoubleIntegral {

    double xStart, yStart, xEnd, yEnd; 
    int xN, yN;
    
    public DoubleIntegral(double xstart, double xend, double ystart, double yend,
            int xn, int yn){
        xStart = xstart;
        yStart = ystart;
        xEnd = xend;
        yEnd = yend;
        xN = xn;
        yN = yn;
        
        double xPrefix = (xEnd - xStart) / xN;
        double yPrefix = (yEnd - yStart) / yN;
        
        double total = 0;
        
        
        
        
        
        
    }
    
        
    private double function(double x, double y){
        double value = x + y ;
        return value;
    }
    
    public static void main() {
        
    }
}
