package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.exception.MethodNotImplementedException;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.IVectorRandom;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Coulomb;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.util.RandomNumberGenerator;

/**
 * Two-centered Lennard Jones molecule with a quadrupole.
 *
 * @author Jayant K. Singh
 */
public class P2LJQ extends Potential2 implements Potential2Soft {

    public P2LJQ(Space space) {
        this(space, 1, 1, 1);
    }

    public P2LJQ(Space space, double sigma, double epsilon,  double momentSquared) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setQuadrupolarMomentSquare(momentSquared);
        gradient = new IVector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        drunit = space.makeVector();
        dcos1dr = space.makeVector();
        dcos2dr = space.makeVector();
    }

    public void setHardCoreDiamterSq(double val){
        hsdiasq=val;
    }

    public void setBox(IBox box) {
        nearestImageTransformer = box.getBoundary();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IAtomSet pair){
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

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }

    public double getEpsilon() {return epsilon;}

    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4*epsilon;
        epsilon48 = 48*epsilon;
    }

    public void setQuadrupolarMomentSquare(double moment){
        Q2=moment;
    }
    
    public double getQuadrupolarMomentSquare() {
        return Q2;
    }
    
    public IVector[] gradient(IAtomSet pair, Tensor pressureTensor) {
        return gradient(pair);
    }

    public IVector[] gradient(IAtomSet pair) {
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)pair.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        double rdudr = epsilon48*s6*(s6 - 0.5);
        gradient[0].Ea1Tv1(rdudr/r2,dr);

        if(Q2!=0.0){
            // normalize dr, the vector between the molecules
            drunit.E(dr);
            double r12 = Math.sqrt(r2);
            drunit.TE(1/r12);

            // v1 is the orientation of molecule 1
            IVector v1 = atom1.getOrientation().getDirection();

            // v2 is the orientation of molecule 2
            IVector v2 = atom2.getOrientation().getDirection();

            // cos1 and sin1 are the cosine and sine of the angle (theta1)
            // between v1 and dr
            double cos1 = v1.dot(drunit);
            // cos2 and sin2 are the cosine and sine of the angle (theta2)
            // between v2 and dr
            double cos2 = v2.dot(drunit);
            // cos12sin1sin2 is the cosine of phi12, the angle between the
            // projections of v1 and v2 onto the plane perpendicular to dr
            // between the molecules multiplied by sin1 and sin2
            
            // temp = sin1*sin2*cos(phi12) - 4*cos1*cos2
            // cos(phi12) = v1 dot v2 - (v1 dot dr) * (v1 dot dr) / (sqrt(1-(v1 dot dr)^2)sqrt(1-(v2 dot dr)^2)
            // [ do some algebra!]
            // temp = v1 dot v2 - 5*cos1*cos2
            double temp = v1.dot(v2) - 5*cos1*cos2;
            
            double uqq = (3.0*Q2)/(4.0*r2*r2*r12);
            double uQuad = 5/r12*(uqq)*(1.0-5.0*(cos1*cos1+cos2*cos2+3*cos1*cos1*cos2*cos2) +2*(temp*temp));
            gradient[0].PEa1Tv1(uQuad, drunit);

            // calculate dcos1 / dr1  (gradient of cos1)
            dcos1dr.Ea1Tv1(-r2, v1);
            dcos1dr.PEa1Tv1(v1.dot(dr), dr);
            dcos1dr.TE(-1/(r12*r2));

            // calculate dcos2 / dr1  (gradient of cos2)
            dcos2dr.Ea1Tv1(-r2, v2);
            dcos2dr.PEa1Tv1(v2.dot(dr), dr);
            dcos2dr.TE(-1/(r12*r2));

            // now actually add terms to the gradients
            gradient[0].PEa1Tv1(uqq*10*cos1, dcos1dr);
            gradient[0].PEa1Tv1(uqq*10*cos2, dcos2dr);
            gradient[0].PEa1Tv1(uqq*20*v1.dot(v2)*cos2, dcos1dr);
            gradient[0].PEa1Tv1(uqq*20*v1.dot(v2)*cos1, dcos2dr);
            gradient[0].PEa1Tv1(-uqq*70*cos1*cos2*cos2, dcos1dr);
            gradient[0].PEa1Tv1(-uqq*70*cos1*cos1*cos2, dcos2dr);
        }

        // pairwise additive, so
        gradient[1].Ea1Tv1(-1,gradient[0]);
        
        return gradient;
    }

    public double virial(IAtomSet atoms) {
        gradient(atoms);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)atoms.getAtom(1);

        // LJ contributation

        dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        nearestImageTransformer.nearestImage(dr);
        double v = gradient[1].dot(dr);
        if (Double.isInfinite(v)) {
            throw new RuntimeException("oops "+v);
        }
        return gradient[1].dot(dr);
    }

    public double hyperVirial(IAtomSet pair) {
        throw new MethodNotImplementedException();
    }

    /**
     * Sets the temperature used for Boltzmann-weighting of the orientational
     * average energy used in u(double) and integral(double)
     */
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    // return a Boltzmann-weighted orientational average of the energy
    // for the given distance
    public double u(double r2) {
        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        double r4 = r2*r2;
        return epsilon4*s6*(s6 - 1.0) - (14.0/(5.0*temperature))*Q2*Q2/(r4*r4*r2);
    }
    
    public double integral(double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        double rc = sigma/rC;
        double sigmaD = space.powerD(sigma);
        double rc3 = rc*rc*rc;
        double rc6 = rc3*rc3;
        double rc12 = rc6*rc6;
        // LJ part
        double uInt = 4.0*epsilon*sigmaD*A*(rc12/9 - rc6/3)/rc3;  //complete LRC is obtained by multiplying by N1*N2/V
        double r2 = rC;
        double r4 = r2*r2;
        return uInt - (7.0/(25.0*temperature))*Q2*Q2/(r4*r4*r2*rC);
    }

    private static final long serialVersionUID = 1L;
    private double sigma , sigma2;
    private double epsilon, epsilon4, epsilon48;
    private double hsdiasq=1.0/Math.sqrt(2);
    private double Q2;
    private NearestImageTransformer nearestImageTransformer;
    private final IVector dr, drunit, dcos1dr, dcos2dr;
    private final IVector[] gradient;
    protected double temperature;

    public static void main(String[] args) {
        //-13.4x10^(-40) C m2
        double Q = Coulomb.UNIT.toSim(13.4E-40)*1E20;
        System.out.println(Q+" "+(Q*Q/Math.pow(3.0354,5)/Kelvin.UNIT.toSim(125.317)));
        System.out.println(Debye.UNIT.fromSim(Q));
        System.out.println(Q*Q/(Math.pow(3.91,5)*Kelvin.UNIT.toSim(148.2)));
        System.out.println("1D "+Debye.UNIT.toSim(1));
        System.out.println("3.33564E-30Cm "+Coulomb.UNIT.toSim(3.33564E-30)*1E10);
        System.exit(1);
        RandomNumberGenerator random = new RandomNumberGenerator();
        Space3D space = Space3D.getInstance();
        Simulation sim = new Simulation(space, false);
        IBox box = new Box(new BoundaryRectangularNonperiodic(space, random), space);
        sim.addBox(box);
        SpeciesSpheresRotating species = new SpeciesSpheresRotating(sim, space);
        sim.getSpeciesManager().addSpecies(species);
        box.setNMolecules(species, 2);
        
        P2LJQ potential = new P2LJQ(sim.getSpace());
        
        IAtomSet leafAtoms = box.getLeafList();
        IAtomOriented atom1 = (IAtomOriented)leafAtoms.getAtom(1);
        potential.setBox(box);
        
        IVector grad1 = space.makeVector();
        IVector oldPosition = space.makeVector();
        IVectorRandom ran = (IVectorRandom)space.makeVector();

        atom1.getPosition().setX(0, 2);
        for (int i=0; i<100; i++) {
            // calculate the gradient
            System.out.println("do gradient");
            IVector[] gradients = potential.gradient(leafAtoms);
            grad1.E(gradients[1]);
            if (grad1.isNaN()) {
                throw new RuntimeException("oops "+grad1);
            }
            oldPosition.E(atom1.getPosition());
            // calculate the gradient numerically for 1 direction
            int d = random.nextInt(3);
            double h = 0.001;
            System.out.println("d="+d);
            atom1.getPosition().setX(d, oldPosition.x(d)+h);
            double uplus = potential.energy(leafAtoms);
            System.out.println("U plus "+uplus);
            atom1.getPosition().setX(d, oldPosition.x(d)-h);
            double uminus = potential.energy(leafAtoms);
            System.out.println("U minus "+uminus);
            double du = (uplus-uminus)/(2*h);
            if (Double.isNaN(du)) {
                throw new RuntimeException("oops "+du+" "+uminus+" "+uplus);
            }
            System.out.println(du+" "+grad1.x(d));
            // check that the analytical and numerical gradients are equal
            if (Math.abs(du-grad1.x(d))/(Math.abs(du)+Math.abs(grad1.x(d))) > 1E-5) {
                System.err.println("hmm");
                if (Math.abs(du-grad1.x(d))/(Math.abs(du)+Math.abs(grad1.x(d))) > 1E-3) {
                    throw new RuntimeException("oops "+du+" "+grad1.x(d));
                }
            }
            else {
                System.out.println("success");
            }
            do {
                // move the atom1 to an entirely random position within 5sigma of atom0
                ((IVectorRandom)atom1.getPosition()).setRandomInSphere(random);
                atom1.getPosition().TE(5);
                ran.setRandomSphere(random);
                atom1.getOrientation().setDirection(ran);
                System.out.println("direction = "+ran);
            }
            while (Double.isInfinite(potential.energy(leafAtoms)));
        }
    }
}
