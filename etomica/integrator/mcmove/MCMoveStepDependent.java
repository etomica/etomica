package etomica.integrator.mcmove;

public interface MCMoveStepDependent {

    /**
     * Sets the step size of the move.
     * 
     * @throws IllegalArgumentException if the given step size is less than the
     * minimum step size or greater than the maximum step size.
     */
    public void setStepSize(double newStepSize);

    /**
     * returns the current step size.
     */
    public double getStepSize();

    /**
     * Sets the maximum allowable step size.
     */
    public void setStepSizeMax(double step);

    /**
     * Returns the maximum allowable steps size.
     */
    public double getStepSizeMax();

    /**
     * Returns the maximum allowable steps size.
     */
    public void setStepSizeMin(double step);

    /**
     * Returns the minimum allowable steps size.
     */
    public double getStepSizeMin();

}