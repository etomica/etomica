package simulate;

public interface AtomWall {

    public void setAngle(int t);
    public int getAngle();
        
    public boolean isVertical();
    public boolean isHorizontal();
    
    public int getThickness();
    public void setThickness(int thickness);
    
    public double getLength();
    public void setLength(double length);
    
    public double getTemperature();
    public void setTemperature(double t);
    
}
