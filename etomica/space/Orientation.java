package etomica.space;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class Orientation {
    public abstract void E(Orientation o); //copies the given orientation to this
    public abstract Vector[] bodyFrame();//body-frame axes in the space-fixed frame
    public abstract double[] angle();//set of angles describing the orientation
    public abstract void rotateBy(double[] t); //rotate all angles by amounts in t array
    public abstract void rotateBy(int i, double dt); //rotate angle i by given amount
    public abstract void randomRotation(double t); //rotate by random amount in solid angle theta on present position
    public abstract void convertToBodyFrame(Vector v); //changes the components of v from space frame to body frame
    public abstract void convertToSpaceFrame(Vector v);//changes the components of v from body frame to space frame
}