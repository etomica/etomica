package simulate;
import java.util.*;
import java.awt.Graphics;

public class MoleculeDiatomic extends Molecule {

    private final double[] e = new double[Space.D];  //normalized orientation vector
    final static int nAtoms = 2;
    final Atom[] atom = new Atom[nAtoms];
    
    public MoleculeDiatomic(int index, double rm, double[] diameter,
                            double L, double[] e) {
        super();
        this.speciesIndex = index;
        double mass = 0.0;
        for(int i=0; i<nAtoms; i++) {
            atom[i] = new Atom();
            atom[i].setRm(rm);
            mass += 1.0/rm;
            atom[i].setDiameter(diameter[i]);
        }
        this.rm = 1.0/mass;
        Space.uEv1(this.e,e);
        setAtomPositions(L,e);
    }
    
    public void setDiameter(int i, double d) {atom[i].setDiameter(d);}

/*    public void zeroForce(int i) {atom[i].zeroForce();}
    public void addForce(int i, double[] force) {atom[i].addForce(force);}
    public void subtractForce(int i, double[] force) {atom[i].subtractForce(force);}

    public void zeroForce() {   //something here causes problem
        super.zeroForce();      //with bean instantiation
        for(int i=nAtoms; --i>=0; ) {atom[i].zeroForce();}
    }
*/    public void translate(int i, double[] dr) {
        atom[i].translate(dr);
        Space.uPEa1Tv1(displacement,this.rm/atom[i].rm,dr);
    }
    
    public void setPosition(double[] rNew) {
        Space.uEv1(this.r,rNew);
        updateOrientation();
        double L = getBondLength();   //inefficient recalc of L
        setAtomPositions(L, e);
    }     
    
    public void setOrientation(double ex) {  //assumes Space.D = 2
        double ey = Math.sqrt(1.0-ex*ex);
        this.setOrientation(ex, ey);
    }
    public void setOrientation(double ex, double ey) { //assumes Space.D=2;
        double L = getBondLength();                    //assumes ex^2+ey^2=1
        e[0] = ex; e[1] = ey;
        setAtomPositions(L, e);
    }
    public void updateOrientation() {
        double rL = 1.0/getBondLength();
        Space.uEa1T_v1Mv2_(this.e,rL,atom[0].r,atom[1].r);
    }
    public double[] getOrientation() {
        updateOrientation();
        return this.e;
    }
      
    public void setBondLength(double L) {
        updateOrientation();
        setAtomPositions(L, e);
    }
    public double getBondLength() {
        double w1 = Math.sqrt(Space.v1Mv2_S(atom[0].r,atom[1].r));
        return w1;
    }
    
    public void setAtomPositions(double L, double[] e) {
        Space.uEv1Pa1Tv2(atom[0].r,this.r,L*this.rm/atom[0].rm,e);
        Space.uEv1Ma1Tv2(atom[1].r,this.r,L*this.rm/atom[1].rm,e);
    }

    public void accelerate(int i, double[] dp) {atom[i].accelerate(dp);}
    public void decelerate(int i, double[] dp) {atom[i].decelerate(dp);}

    public double getKineticEnergy() {
        return atom[0].getKineticEnergy() + atom[1].getKineticEnergy();
    }
}