package etomica.lattice.crystal;

public interface Primitive3D {
    
    //primitive-vector lengths
    public void setA(double a);
    public double getA();
    
    public void setB(double b);
    public double getB();
    
    public void setC(double c);
    public double getC();
    
    //flags indicating which fields can be edited independently
    public boolean isEditableA();
    public boolean isEditableB();
    public boolean isEditableC();
    public boolean isEditableAlpha();
    public boolean isEditableBeta();
    public boolean isEditableGamma();
    //angles
    /**
     * alpha is the angle between vectors 1 and 2 (b and c)
     */
    public void setAlpha(double t);
    public double getAlpha();
    
    /**
     * beta is the angle between vectors 0 and 2 (a and c)
     */
    public void setBeta(double t);
    public double getBeta();
    
    /**
     * gamma is the angle between vectors 0 and 1 (a and b)
     */
    public void setGamma(double t);
    public double getGamma();
    
}