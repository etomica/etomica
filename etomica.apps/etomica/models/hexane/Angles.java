package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomTreeNodeGroup;
import etomica.phase.Phase;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * ASSUMES 3D!!!
 * 
 * @author cribbin
 * 
 */
public class Angles {

    public Angles(Phase phase) {
        dim = phase.space().D();
        if (dim != 3) {
            throw new RuntimeException("Angles class requires a 3D space!");
        }
        ang = phase.space().makeVector();
        vex = phase.space().makeVector();

        vex0 = phase.space().makeVector();
        vex1 = phase.space().makeVector();
        vex2 = phase.space().makeVector();
        vex3 = phase.space().makeVector();

        vex4 = phase.space().makeVector();
        vex5 = phase.space().makeVector();
        vex6 = phase.space().makeVector();
        vex7 = phase.space().makeVector();
   
        axis0 = phase.space().makeVector();
        axis1 = phase.space().makeVector();
        axis2 = phase.space().makeVector();
        
        answer = new double[dim];
        
//        vex0 = new Vector[dim];
//        vex1 = new Vector[dim];
//        for(int i = 0; i < dim; i++){
//            vex0[i] = phase.space().makeVector(); 
//            vex1[i] = phase.space().makeVector();
//        }
        
        alpha = 0.0;
        beta = 0.0;
        gamma = 0.0;
        alphaRef = 0.0;
        betaRef = 0.0;
        gammaRef = 0.0;
        alphaMeas = 0.0;
        betaMeas = 0.0;
        gammaMeas = 0.0;

        temp = phase.space().makeVector();
        rot = phase.space().makeTensor();
        
    }

    
    public double[] getAngles(AtomPair atompair){
        return getAngles(atompair.getAtom(0), atompair.getAtom(1));
    }
    
    public double[] getAngles(Atom atom0, Atom atom1){
        //make sure they're zero at the start.
        answer[0] = 0.0;
        answer[1] = 0.0;
        answer[2] = 0.0;
        
        //Set up all the axes based on the molecule atom0, the reference molecule
        //Long rotational axis of atom 0
        vex0.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(0)).coord.position() );
        vex1.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(1)).coord.position() );
        vex2.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(2)).coord.position() );
        vex3.E(vex2);
        vex3.ME(vex0);
        //vex3 should now go from the 0th atom on the molecule to the 2nd atom on the molecule
        axis0.E(vex3);
        //Now we take the midpoint of vex3.
        vex3.E(vex2);
        vex3.PE(vex0);
        vex3.TE(-0.5);
        //Then we subtract the location of the 1st atom on the molecule to get our final vector.
        vex3.PE(vex1);
        
        axis1.E(vex3);  //AKA isosceleshappyvector
        
        axis2.E(axis0.cross(axis1));
        
        //Normalize our axes
        double l0, l1, l2;
        l0 = getVectorLength(axis0);
        l1 = getVectorLength(axis1);
        l2 = getVectorLength(axis2);
        
        axis0.TE(1/l0);
        axis1.TE(1/l1);
        axis2.TE(1/l2);
        
        //Now we play with the molecule we are measuring.
        
        //Long rotational axis of atom 1
        vex0.E( ((AtomLeaf)((AtomTreeNodeGroup)atom1.node).childList.get(0)).coord.position() );
        vex1.E( ((AtomLeaf)((AtomTreeNodeGroup)atom1.node).childList.get(1)).coord.position() );
        vex2.E( ((AtomLeaf)((AtomTreeNodeGroup)atom1.node).childList.get(2)).coord.position() );
        vex3.E(vex2);
        vex3.ME(vex0);
        //vex3 should now go from the 0th atom on the molecule to the 2nd atom on the molecule
        
        //now we do the whole projection thing!!  notes 11/3/06 will tell you why this is possible
        // mostly because the axes are  unit vectors.
        //Project vex3 onto axis1 and axis2
        answer[1] = vex3.dot(axis1);
        answer[2] = vex3.dot(axis2);
        
        //do the match and rotate thing
        //pull the reference molecule's tail to the origin
        vex0.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(0)).coord.position() );
        vex1.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(1)).coord.position() );
        vex2.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(2)).coord.position() );
        vex2.ME(vex0);
        vex1.ME(vex0);
        vex0.ME(vex0);
        
        //pull the measured molecule's tail to the origin
        vex4.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(0)).coord.position() );
        vex5.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(1)).coord.position() );
        vex6.E( ((AtomLeaf)((AtomTreeNodeGroup)atom0.node).childList.get(2)).coord.position() );
        vex6.ME(vex4);
        vex5.ME(vex4);
        vex4.ME(vex4);
        
        //ROTATION of the measured molecule to align with the reference molecule
        
        //First, find the angles of each backbone from each axis.
        //we need to rotate the measured atom to match the reference 
        
        //alpha is the rotation around the x axis
        //beta is the rotation around the y axis
        //gamma is the rotation around the z axis
        
        ((Vector3D)temp).E(1.0, 0.0, 0.0);
        alphaRef = dotForAngle(temp, vex3);
        
        ((Vector3D)temp).E(1.0, 0.0, 0.0);
        alphaMeas = dotForAngle(temp, vex7);
        
        alpha = alphaMeas - alphaRef;
        
        ((Vector3D)temp).E(0.0, 1.0, 0.0);
        betaRef = dotForAngle(temp, vex3);
        
        ((Vector3D)temp).E(0.0, 1.0, 0.0);
        betaMeas = dotForAngle(temp, vex7);
        
        beta = betaMeas - betaRef;
        
        ((Vector3D)temp).E(0.0, 0.0, 1.0);
        gammaRef = dotForAngle(temp, vex3);
        
        ((Vector3D)temp).E(0.0, 0.0, 1.0);
        gammaMeas = dotForAngle(temp, vex7);
        
        gamma = gammaMeas - gammaRef;
        
        //Rotate the measured vectors  //nan we could collapse these into a single tensor
        //Rotate the vector around the x axis
        double[] dork0 =  {  1.0, 0.0, 0.0,
                            0.0, Math.cos(alpha), -1.0*Math.sin(alpha), 
                            0.0, Math.sin(alpha), Math.cos(alpha)   
                            };
        rot.E(dork0);
        vex4.transform(rot);
        vex5.transform(rot);
        vex6.transform(rot);
        
        //rotate the vector around the y axis
        double[] dork1 = {    Math.cos(beta), 0.0, Math.sin(beta),
                            0.0, 1.0, 0.0,
                            -1.0*Math.sin(beta), 0.0, Math.cos(beta)    
                };
        rot.E(dork1);
        vex4.transform(rot);
        vex5.transform(rot);
        vex6.transform(rot);
        
        //rotate the vector around the z axis
        double[] dork2 = {    Math.cos(gamma), -1.0*Math.sin(gamma), 0.0,
                        Math.sin(gamma), Math.cos(gamma), 0.0,
                        0.0, 0.0, 1.0
                        };
        rot.E(dork2);
        vex4.transform(rot);
        vex5.transform(rot);
        vex6.transform(rot);
        
        //Calculate the happyisocelesvector of the reference molecule
        vex3.E(vex2);
        vex3.PE(vex0);
        vex3.TE(-0.5);
        vex3.PE(vex1);
        
        //Calculate the happyisocelesevector of the measured molecule
        vex7.E(vex6);
        vex7.PE(vex4);
        vex7.TE(-0.5);
        vex7.PE(vex5);
        
        //Calculate the angle between the happyisocelesvectors
        answer[0] = dotForAngle(vex7, vex3);
        
        //Throw a party!
     
        return answer;
    }

    private double getVectorLength(Vector v) {
        double l = 0.0;
        for (int i = 0; i < dim; i++) {
            l += v.x(i) * v.x(i);
        }
        return Math.sqrt(l);
    }

    private double dotForAngle(Vector v, Vector u) {
        double l = 0.0;
        for (int i = 0; i < dim; i++) {
            l += u.x(i) * v.x(i);
        }
        l = l / getVectorLength(u) / getVectorLength(v);
        return Math.acos(l);
    }


    private Vector vex;

    private int dim;
    private Vector ang; // holds the delta angles
    private Vector vex0;
    private Vector vex1;
    private Vector vex2;
    private Vector vex3;
    private Vector vex4;
    private Vector vex5;
    private Vector vex6;
    private Vector vex7;
    double[] answer;
    
    private Vector axis0;
    private Vector axis1;
    private Vector axis2;
    
    private double alpha, beta, gamma;
    private double alphaRef, betaRef, gammaRef;
    private double alphaMeas, betaMeas, gammaMeas;
    
    Vector temp;
    Tensor rot;
}

