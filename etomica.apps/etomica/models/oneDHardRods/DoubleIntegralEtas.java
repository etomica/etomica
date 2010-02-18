package etomica.models.oneDHardRods;

public class DoubleIntegralEtas {

    public DoubleIntegralEtas(double xStart, double xEnd, double yStart, double yEnd,
            int xN, int yN){
        double total = 0.0;
        double xValue = 0.0;
        double yValue = 0.0;
        
        //we skip i = 0 because the endpoints do not double.
        for(int i = 1; i < xN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            for (int j = 1; j < yN; j++) {
                yValue = yStart + j *(yEnd - yStart) / yN;
                total += 4 * function(xValue, yValue);
            }
        }
        
        //Now do the "edge" stuff, where either the x or y value, but not both,
        //  are the start or end value
        xValue = xStart;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            total += 2 * function(xValue, yValue);
        }
        xValue = xEnd;
        for (int j = 1; j < yN; j++) {
            yValue = yStart + j *(yEnd - yStart) / yN;
            total += 2 * function(xValue, yValue);
        }
        yValue = yStart;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            total += 2 * function(xValue, yValue);
        }
        yValue = yEnd;
        for (int i = 1; i < yN; i++) {
            xValue = xStart + i *(xEnd - xStart) / xN;
            total += 2 * function(xValue, yValue);
        }
        
        //Now the four values for which both the x and y values are
        //  the start or end value
        xValue = xStart;
        yValue = yStart;
        total += function(xValue, yValue);
        
        yValue = yEnd;
        total += function(xValue, yValue);
        
        xValue = xEnd;
        yValue = yStart;
        total += function(xValue, yValue);
        
        yValue = yEnd;
        total += function(xValue, yValue);
        
        //Now we do the prefix thing
        double prefix = (xEnd - xStart)*(yEnd - yStart) / (4* xN * yN);
        
        total *= prefix;
        
        System.out.println("Integral = " + total);
    }
    
        
    private double function(double x, double y){
        double value = x * y ;
        return value;
    }
    
    public static void main(String[] args) {
      double xStart = -1.0;
      double yStart = -1.0;
      double xEnd = 1.0;
      double yEnd = 1.0;
      int xN = 100000;
      int yN = 100000;
        

      DoubleIntegralEtas di = new DoubleIntegralEtas(xStart, xEnd, yStart, yEnd, xN,yN);
        
    }
}
