package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.exception.MethodNotImplementedException;
import etomica.space.IVector;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Two-centered Lennard Jones molecule with a quadrupole.
 *
 * @author Jayant K. Singh
 */
public class P2LJQ extends Potential2 implements Potential2Soft {

    public P2LJQ(Space space) {
        this(space, 1, 1, 1);
    }

    public P2LJQ(Space space, double sigma, double epsilon,  double moment) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setQuadrupolarMomentSquare(moment);
        dr = space.makeVector();
    }

    public void setHardCoreDiamterSq(double val){
        hsdiasq=val;
    }

    public void setBox(Box box) {
        nearestImageTransformer = box.getBoundary();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(AtomSet pair){
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)pair.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        if(r2<hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        double s2 = sigma2/(r2);
        double s6 = s2*s2*s2;
        double ener = epsilon4*s6*(s6 - 1.0);

        if(Q2!=0.0){
            // normalize dr, the vector between the molecules
            dr.normalize();

            // v1 is the orientation of molecule 1
            IVector v1 = atom1.getOrientation().getDirection();

            // v2 is the orientation of molecule 2
            IVector v2 = atom2.getOrientation().getDirection();

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and dr
            double cos1 = v1.dot(dr);
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and dr
            double cos2 = v2.dot(dr);
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to dr
            // between the molecules multiplied by sin1 and sin2
            
            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double temp = v1.dot(v2) - 5*cos1*cos2;
            
            double uqq = (3.0*Q2)/(4.0*r2*r2*Math.sqrt(r2));
            double uQuad = (uqq)*(1.0-5.0*(cos1*cos1+cos2*cos2+3*cos1*cos1*cos2*cos2) +2*(temp*temp));
            ener += uQuad;
        }
        return ener;
    }

    protected double calculateEnergy(double r2){
        if(r2<hsdiasq) {
            return Double.POSITIVE_INFINITY;
        }
        double s2 = sigma2/(r2);
        double s6 = s2*s2*s2;
        return epsilon4*s6*(s6 - 1.0);
    }

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }

    public double getEpsilon() {return epsilon;}

    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4*epsilon;
    }

    public final void setQuadrupolarMomentSquare(double moment){
        Q2=moment;	
    }
    
    public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
        throw new MethodNotImplementedException();
    }

    public IVector[] gradient(AtomSet atoms) {
        throw new MethodNotImplementedException();
    }

    public double virial(AtomSet atoms) {
        throw new MethodNotImplementedException();
    }

    public double hyperVirial(AtomSet pair) {
        throw new MethodNotImplementedException();
    }

    public double integral(double rc) {
        // TODO Auto-generated method stub
        return 0;
    }

    private static final long serialVersionUID = 1L;
    private double sigma , sigma2;
    private double epsilon, epsilon4;
    private double hsdiasq=1.0/Math.sqrt(2);
    private double Q2;
    private NearestImageTransformer nearestImageTransformer;
    private final IVector dr;
}
