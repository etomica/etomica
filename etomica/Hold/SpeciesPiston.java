package simulate; 
    
public class SpeciesPiston extends SpeciesWall {
    
    public double density, pressure;
    
    public SpeciesPiston() {
        super();
    }
    
    public double getDensity() {return density;}
    public void setDensity(double d) {
        density = d / 1e30 * 1e3 * Constants.AVOGADRO;
        setMass(density*firstAtom.diameter*firstAtom.diameter*1e10);
        computeForce();
    }
    
    public double getPressure() {return pressure/Constants.BAR2SIM;}
    public void setPressure(double p) {
        pressure = p*Constants.BAR2SIM;
        computeForce();
    }
    public void setPressureI(int p) {setPressure((double)p);}
    
    public double getVelocity() {return firstAtom.p[coordinateIndex]*firstAtom.rm;}
    public void setVelocity(double v) {
        firstAtom.p[coordinateIndex]=firstAtom.mass*v;}
    
    public final double getMass() {return firstAtom.getMass();}
    public void setMass(double m) {firstAtom.setMass(m);}
    
    public void setHorizontal(boolean b) {
        super.setHorizontal(b);
        computeForce();
    }
    public void setVertical (boolean b) {
        super.setVertical(b);
        computeForce();
    }
    
    void computeForce() {
        if(((AtomWall)firstAtom).isHorizontal()) {
            firstAtom.f[0] = 0.0;
            firstAtom.f[1] = getMass()*Constants.G + pressure*firstAtom.diameter*firstAtom.diameter;
            firstAtom.setForceFree(firstAtom.f[1] == 0.0);
        }
        else {
            firstAtom.f[0] = pressure*firstAtom.diameter*firstAtom.diameter;
            firstAtom.f[1] = 0.0;
            firstAtom.setForceFree(firstAtom.f[0] == 0.0);
        }
    }
}