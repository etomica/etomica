package etomica;

public interface PotentialCalculation {
    
    public interface Sum extends PotentialCalculation {
        public double sum();
    }
}//end of PotentialCalculation