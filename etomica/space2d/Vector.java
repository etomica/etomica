package simulate.space2D;
import simulate.Space;
import java.util.Random;
    
public final class Vector extends Space.Vector {  //declared final for efficient method calls
    public static final Random random = new Random();
    public static final Vector ORIGIN = new Vector(0.0,0.0);  //anything using WORK is not thread-safe
    public static final Vector WORK = new Vector();
    public double x, y;
    public Vector () {x = 0.0; y = 0.0;}
    public Vector (double x, double y) {this.x = x; this.y = y;}
    public Vector (double[] a) {x = a[0]; y = a[1];}//should check length of a for exception
    public int length() {return Space2D.D;}
    public int D() {return Space2D.D;}
    public double component(int i) {return (i==0) ? x : y;}
    public void setComponent(int i, double d) {if(i==0) x=d; else y=d;}
    public void E(Vector u) {x = u.x; y = u.y;}
    public void E(double a) {x = a; y = a;}
//        public void E(int i, double a) {if(i==0) x = a; else y = a;}  //assumes i = 0 or 1
    public void Ea1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x = a1*u1.x; y = a1*u1.y;}
    public void PEa1Tv1(double a1, Space.Vector u) {Vector u1=(Vector)u; x += a1*u1.x; y += a1*u1.y;}
    public void PE(Vector u) {x += u.x; y += u.y;}
    public void ME(Vector u) {x -= u.x; y -= u.y;}
    public void PE(int i, double a) {if(i==0) x += a; else y += a;}
    public void TE(double a) {x *= a; y *= a;}
    public void TE(int i, double a) {if(i==0) x *= a; else y *= a;}
    public void DE(double a) {x /= a; y /= a;}
    public void DE(Vector u) {x /= u.x; y /= u.y;}
    public double squared() {return x*x + y*y;}
    public double dot(Vector u) {return x*u.x + y*u.y;}
    public void normalize() {
        double norm = Math.sqrt(1/(x*x + y*y));
        x *= norm;
        y *= norm;
    }
    public void randomStep(double d) {x += (2.*random.nextDouble()-1.0)*d; y+= (2.*random.nextDouble()-1.0)*d;} //uniformly distributed random step in x and y, within +/- d
    public void setRandom(double d) {x = random.nextDouble()*d; y = random.nextDouble()*d;}
    public void setRandom(double dx, double dy) {x = random.nextDouble()*dx; y = random.nextDouble()*dy;}
    public void setRandom(Vector u) {setRandom(u.x,u.y);}
    public void setRandomCube() {
        x = random.nextDouble() - 0.5; 
        y = random.nextDouble() - 0.5;
    }
    public void setRandomSphere() {
        x = Math.cos(2*Math.PI*random.nextDouble()); 
        y = Math.sqrt(1.0 - x*x);
        if(random.nextDouble() < 0.5) y = -y;
    }
    public void randomDirection() {
        x = Math.cos(Math.PI*random.nextDouble()); 
        y = Math.sqrt(1.0 - x*x);
        if(random.nextDouble() < 0.5) y = -y;
    }
    public void E(Space.Vector u) {E((Vector)u);}
    public void PE(Space.Vector u) {PE((Vector)u);}
    public void ME(Space.Vector u) {ME((Vector)u);}
    public void TE(Space.Vector u) {TE((Vector)u);}
    public void DE(Space.Vector u) {DE((Vector)u);}
    public double dot(Space.Vector u) {return dot((Vector)u);}
}
    
