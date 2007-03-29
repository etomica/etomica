package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.integrator.IntegratorMC;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.space.Tensor;
import etomica.species.Species;
import etomica.util.IRandom;

/**
 * Contains specifics for the solid hexane model:
 * <ul>
 * <li>bond length = 0.4
 * <li>rejects moving the whole molecule as a CBMC move
 * </ul>
 * @author cribbin
 *
 */
public class CBMCGrowSolidHexane extends CBMCGrowStraightAlkane {

    public CBMCGrowSolidHexane(PotentialMaster p, IRandom random, IntegratorMC integrator, Phase phs, Species species, int NTrials){
        super(p, random, integrator, phs, species, 6, NTrials);
        
        if(p.getSpace().D() != 3){
            throw new IllegalArgumentException("Torsional bond is only used in 3D simulations");
        }
        
        setBondLength(0.4);
//        setPrefactor(1.0);
        
        phi = (180 - 109.47) / 360.0 * 2 * Math.PI; //makes sure the vector is pointing in the right direction on the cosine section
        
//        pots = p.getPotentials();
//        potMast = p;
    
        rotor = p.getSpace().makeTensor();
        temp2 = p.getSpace().makeVector();
    }
    
    protected int calcStartAtom(){
        return random.nextInt(chainlength-2) + 1 ;//function will not pick chainlength;
    }
    
    //Different because we know the bond angle
    //All moves are accepted
    protected IVector calcRandomBondWithAngle(AtomLeaf a, AtomLeaf b){
        //temp will be the radial vector
        //vex will be the axial vector
        
        temp.E(b.getCoord().getPosition());
        temp.ME(a.getCoord().getPosition());
        vex.E(temp);     //Store it because we need it, and we are going to change temp in the next line of code.
       
        temp.E(getNormal(temp));

//      Subtract the projection of b-c(temp) onto a-b(vex) from temp,
        // which leaves the radial part of bc as temp
        vex.normalize();    
        //This is the projection bit. (equivalent to multiplying by the cosine of the angle between the radial and axial vectors)
        vex.TE(vex.dot(temp));
        vex.TE(1/getBondlength()/getBondlength());
        //This is the subtraction bit.
        temp.ME(vex);
        
        //Create the rotation matrix for an arbitrary unit vector and an 
        // arbitrary angle, and apply it
        double randomAngle = random.nextDouble() * 2 * Math.PI;
        double cosRA = Math.cos(randomAngle);
        double sinRA = Math.sin(randomAngle);
        
        rotor.E(new double[]{
                cosRA+(1-cosRA)*vex.x(0)*vex.x(0),                  //xx
                (1-cosRA)*vex.x(0)*vex.x(1)-sinRA*vex.x(2),         //xy
                (1-cosRA)*vex.x(0)*vex.x(2)+sinRA*vex.x(1),         //xz
                (1-cosRA)*vex.x(1)*vex.x(0)+sinRA*vex.x(2),         //yx
                cosRA+(1-cosRA)*vex.x(1)*vex.x(1),                  //yy
                (1-cosRA)*vex.x(1)*vex.x(2)-sinRA*vex.x(0),         //yz
                (1-cosRA)*vex.x(2)*vex.x(0)-sinRA*vex.x(1),         //zx
                (1-cosRA)*vex.x(2)*vex.x(1)+sinRA*vex.x(0),         //zy
                cosRA+(1-cosRA)*vex.x(2)*vex.x(2)                   //zz
        });
                
        //Mulitply the rotation tensor by temp to get the rotated vector.
        //We then have a vector perpendicular to the a-b axis, rotated a random
        //angle.  This is the radial vector.
        rotor.transform(temp);
        
        //we normalize our radial vector (temp). Our axial vector (vex) is 
        //  already normalized.
        temp.normalize();
        
        //we multiply the axial vector by the axial length and the radial
        // vector by the radial length (worked out by hand, based on a bond
        // length of 0.4 and a bond angle of 109.47 
        temp.TE(getBondlength() * Math.cos(phi));
        vex.TE(getBondlength() * Math.sin(phi));
        
        //Add these together, and return that value:
        temp.PE(vex);
        return temp;

    }
    
    protected IVector calcRandomBondWithAngleAndTorsion(AtomLeaf a, AtomLeaf b, 
            AtomLeaf c){
        //Get a random number, and place it between the limits on the new atom's
        // placement.
        //The angle must be between 108.6919204 degrees, and 251.23080979 deg.
        double randomAngle = (251.23080979 - 108.6919204) * random.nextDouble() + 108.6919204;
//        System.out.println(randomAngle);
        //we convert from degrees to radians:
        randomAngle = randomAngle * 2 * Math.PI / 360.0;
        double cosRA = Math.cos(randomAngle);
        double sinRA = Math.sin(randomAngle);
        
        //get a normal vector to the a-b vector and the b-c vector
        //This vector is, by definition, perpendicular to the a-b vector, which
        // makes it a radius of a circle centered on that axis.
        vex.E(b.getCoord().getPosition());
        vex.ME(a.getCoord().getPosition());
        temp.E(c.getCoord().getPosition());
        temp.ME(b.getCoord().getPosition());
        
        //get the cosine of the angle between a-to-b (vex, axial vector)
        // and b-c (temp, radial vector)
//        double cosAng = vex.dot(temp);
//        cosAng /= (getBondlength() * getBondlength());
        
        //Subtract the projection of b-c(temp) onto a-b(vex) from temp,
        // which leaves the radial part of bc as temp
        vex.normalize();
//        vex.TE(getBondlength()*cosAng);     
        //This is the projection bit. (equivalent to multiplying by the cosine of the angle between the radial and axial vectors)
        vex.TE(vex.dot(temp));
        vex.TE(1/getBondlength()/getBondlength());
        temp.ME(vex);           //This is the subtraction bit.
        
        //Create the rotation matrix for an arbitrary unit vector
        vex.normalize();
        rotor.E(new double[]{
                cosRA+(1-cosRA)*vex.x(0)*vex.x(0),                  //xx
                (1-cosRA)*vex.x(0)*vex.x(1)-sinRA*vex.x(2),         //xy
                (1-cosRA)*vex.x(0)*vex.x(2)+sinRA*vex.x(1),         //xz
                (1-cosRA)*vex.x(1)*vex.x(0)+sinRA*vex.x(2),         //yx
                cosRA+(1-cosRA)*vex.x(1)*vex.x(1),                  //yy
                (1-cosRA)*vex.x(1)*vex.x(2)-sinRA*vex.x(0),         //yz
                (1-cosRA)*vex.x(2)*vex.x(0)-sinRA*vex.x(1),         //zx
                (1-cosRA)*vex.x(2)*vex.x(1)+sinRA*vex.x(0),         //zy
                cosRA+(1-cosRA)*vex.x(2)*vex.x(2)                   //zz
        });
                
        //Mulitply the rotation tensor by temp to get the rotated vector.
        //We then have a vector perpendicular to the a-b axis, rotated a random
        //angle.  This is the radial vector.
        rotor.transform(temp);
        
        //we normalize our radial vector (temp). Our axial vector (vex) is 
        //  already normalized.
        temp.normalize();
        
        //we multiply the axial vector by the axial length and the radial
        // vector by the radial length (worked out by hand, based on a bond
        // length of 0.4 and a bond angle of 109.47
        temp.TE(getBondlength() * Math.cos(phi));
        vex.TE(getBondlength() * Math.sin(phi));
        
        //Add these together, and return that value:
        temp.PE(vex);
        return temp;
    }
    
    protected double calcBondAngleEnergy(double d){
        throw new RuntimeException("calcBondAngleEnergy should not be called in CBMCGrowSolidHexane");
        /*
         * We are using a set bond length, and a rigid bond angle,
         * we don't need to calculate an energy because we wrote this program
         * to only use the proper bond length and angle
         */
    }
    
    protected double calcBondTorsionalEnergy(IVector v){
        throw new RuntimeException("calcBondTorsionalEnergy should not be called in CBMCGrowSolidHexane");
        /*
         * We are using a set bond length, and a rigid bond angle,
         * we don't need to calculate an energy because we wrote this program
         * to only use the proper bond length and angle
         */
    }
    
    protected double calcExternalEnergy(Atom a){
        //make sure that the proper potential is enabled.  Really,
        //we only have the one potential, so this line is unnecessary
        //but I want it in here for reference when I am extending
        //this code.
//        potMast.setEnabled(pots[0], true);
        externalMeter.setPhase(phase);
        externalMeter.setTarget(a);
        
        double blind = externalMeter.getDataAsScalar();
//        potMast.setEnabled(pots[0], false);
        
        return blind;
    }
    //Since the bondlength is constant, this returns an ordinary getter
    protected double calcBondL(){
        return getBondlength();
    }
    /**
     * Returns a unit normal to the argument vector
     * @param vector
     * @return
     */
    protected IVector getNormal(IVector vect){
        //Determine the smallest component
        int min = 0;
        if(vect.x(1) < vect.x(0)) {min = 1;}
        if(vect.x(2) < vect.x(min)) {min = 2;}

        //create the unit vector in that direction
        temp2.E(0.0);
        temp2.setX(min, 1.0);
        return temp2;
    }
    
    public void setBondLength(double d){
        bondlength = d;
    }
    public double getBondlength(){
        return bondlength;
    }

    public double energyChange(){
        return 0.0;
    }
    
    private static final long serialVersionUID = 1L;
//    private Potential[] pots;
//    private PotentialMaster potMast;
    private double phi;
    Tensor rotor;
    IVector temp2;
}