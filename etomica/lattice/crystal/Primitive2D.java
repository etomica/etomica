package etomica.lattice.crystal;

public interface Primitive2D {
    
    public void setA(double a);
    public double getA();
    
    public void setB(double b);
    public double getB();
    
    //angles
    public void setAlpha(double t);
    public double getAlpha();
}