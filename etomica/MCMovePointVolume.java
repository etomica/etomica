package etomica;

import java.util.Random;
import etomica.units.Dimension;
import etomica.lattice.*;
import etomica.utility.OdeSolver;
import java.awt.*;

/** 
 * Monte Carlo move that changes the system volume by distorting the space, causing
 * more volume to be place at one point selected at random.  Space is distorted from
 * this point until reaching the boundaries, where it becomes a uniform translation.
 *
 * @author Nandou Lu
 * @author David Kofke
 */

public class MCMovePointVolume extends MCMove implements DisplayPhase.DrawingObject {
    
    private final Random rand = new Random();
    protected double pressure;
    protected PhaseAction.Inflate inflate;
    MyLattice lattice;

    public MCMovePointVolume() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(Default.PRESSURE);
        lattice = new MyLattice(6,new VelocityField(2,20), 0.015);
    }
        
    public void thisTrial() {
        double hOld, hNew, vOld, vNew;
        vOld = phase.volume();
        hOld = phase.energy.potential() + pressure*vOld;
        
        if(Math.random() < 0.5) { //transform square to distorted
        }
        else { //transform distorted back to square
        }
        
        double vScale = (2.*Math.random()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/(double)phase.parentSimulation().space().D());
        inflate.actionPerformed(phase,rScale);
        hNew = phase.energy.potential() + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+(phase.moleculeCount+1)*vScale)
                < Math.random()) 
            {  //reject
              inflate.retractAction();
            }
        nAccept++;   //accept
    }
    
    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
    
    /**
     * Method to draw the grid of original/deformation points.
     * To use, add the instance of this class to the displayPhase via addDrawingObject.
     */
    public void draw(Graphics g, int[] origin, double scale) {
        int toPixels = (int)(scale * etomica.units.BaseUnit.Length.Sim.TO_PIXELS);
        double diameter = 1.0;
        int sigmaP = (int)(toPixels*diameter);
        double xL = phase.boundary().dimensions().component(0)/lattice.scale;
        double yL = phase.boundary().dimensions().component(1);
        for(int i=0; i<lattice.size; i++) {
            for(int j=0; j<lattice.size; j++) {
                MySite site = lattice.sites[i][j];
                Space.Vector r = site.originalPosition;
                int xP = origin[0] + (int)(toPixels*(xL*r.component(0)-0.5*diameter));
                int yP = origin[1] + (int)(toPixels*(yL*r.component(1)-0.5*diameter));
                g.setColor(java.awt.Color.red);
                g.fillOval(xP,yP,sigmaP,sigmaP);
                
                r = site.deformedPosition;
                xP = origin[0] + (int)(toPixels*(xL*r.component(0)-0.5*diameter));
                yP = origin[1] + (int)(toPixels*(yL*r.component(1)-0.5*diameter));
                g.setColor(java.awt.Color.black);
                g.fillOval(xP,yP,sigmaP,sigmaP);
            }
        }
    }
    
    public static void main(String args[]) {
        
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        
        Default.ATOM_SIZE = 1.0;
        
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        Species species = new SpeciesDisks(sim);
        Potential2 potential = new P2SimpleWrapper(sim,new PotentialHardDisk(sim));
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom();
        MCMovePointVolume mcMovePointVolume = new MCMovePointVolume();
        Controller controller = new Controller(sim);
        Phase phase = new Phase(sim);
        DisplayPhase displayPhase = new DisplayPhase(sim);
        displayPhase.setScale(0.8);
        DeviceSlider slider = new DeviceSlider(mcMovePointVolume, "pressure");
        slider.setMinimum(0);
        slider.setMaximum(100);
        
        mcMovePointVolume.setPhase(phase);
        integrator.add(mcMoveAtom);
        integrator.add(mcMovePointVolume);

		MyLattice lattice = mcMovePointVolume.lattice;
        species.setNMolecules(0);
        
        displayPhase.addDrawingObject(mcMovePointVolume);
        
		Simulation.instance.elementCoordinator.go();
		
		phase.boundary().dimensions().TE(0,lattice.scale);
		
/*		ColorSchemeBySpecies colorScheme = new ColorSchemeBySpecies();
		colorScheme.addSpecies(species, new ColorScheme.Simple(java.awt.Color.red));
		colorScheme.addSpecies(deformed, new ColorScheme.Simple(java.awt.Color.black));
		displayPhase.setColorScheme(colorScheme);
		
		Atom atom = phase.firstAtom();
		for(int i=0; i<lattice.size; i++) {
		    for(int j=0; j<lattice.size; j++) {
		        atom.coordinate.position().E(lattice.sites[i][j].originalPosition);
		        atom.coordinate.position().TE(phase.boundary().dimensions());
		        atom = atom.nextAtom();
		    }
		}
		for(int i=0; i<lattice.size; i++) {
		    for(int j=0; j<lattice.size; j++) {
		        atom.coordinate.position().E(lattice.sites[i][j].deformedPosition);
		        atom.coordinate.position().TE(phase.boundary().dimensions());
		        atom = atom.nextAtom();
		    }
		}
*/		        
        f.add(Simulation.instance.panel()); //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
    
    /**
     * Array of 2D triangular cells formed by bisecting the cells of a square lattice.
     */
    private class MyLattice implements AbstractLattice {
        
        SquareLattice squareLattice;
        SiteIterator.List iterator = new SiteIterator.List();
        MySite[][] sites;
        int size;
        double scale;
        /**
         * Size is the number of sites in each dimension; considered non-periodically
         */
        MyLattice(int size, final VelocityField vField, double deltaT) {
            this.size = size;
            squareLattice = new SquareLattice(size, new MySiteFactory(), 1.0/(double)(size-1));
            sites = new MySite[size][size];
            int[] index = new int[2];
            for(index[0]=0; index[0]<size; index[0]++) { //make array of sites for convenience of access
                for(index[1]=0; index[1]<size; index[1]++) {
                    sites[index[0]][index[1]] = (MySite)squareLattice.site(index);
                }
            }
            //assign site positions and make vectors for deformed positions
            for(int i=0; i<size; i++) {    
                for(int j=0; j<size; j++) {
                    MySite site = sites[i][j];
                    site.originalPosition = (Space2D.Vector)((BravaisLattice.Coordinate)site.coordinate()).position();
         //           site.originalPosition.PE(-0.5);
                    site.deformedPosition = new Space2D.Vector();
                    site.deformedPosition.E(site.originalPosition);
                } 
            }
            for(int i=0; i<size-1; i++) {    //loops don't do last row/column of lattice
                for(int j=0; j<size-1; j++) {
                    MySite site = sites[i][j];
                    site.c1 = new MyCell(site, sites[i][j+1], sites[i+1][j+1]);
                    site.c2 = new MyCell(site, sites[i+1][j], sites[i+1][j+1]);
                } 
            }
            
            //Deform lattice
            //make a local rhs class from the velocity field suitable for input to the ode solver
            OdeSolver.Rhs rhs = new OdeSolver.Rhs() {
                public double[] dydx(OdeSolver.Variables xy) {
                    return vField.v(xy.y);
                }
            };
            for(int i=0; i<size; i++) {
                for(int j=0; j<size; j++) {
                    MySite site = sites[i][j];
                    OdeSolver.Variables xy0 = new OdeSolver.Variables(0.0, site.originalPosition.toArray());
         //           System.out.println(i + "  " + j + " " + xy0.y[0] + " " + xy0.y[1]);
                    OdeSolver.Variables[] xy = OdeSolver.rungeKuttaAdaptive(xy0,deltaT,1.e-7,rhs);
                    site.deformedPosition.E(xy[xy.length-1].y);
                } 
            }
            //reshape to square
            double avg = 0.0;
            for(int i=0; i<size; i++) {
                double x0 = sites[0][i].deformedPosition.component(0);
                for(int j=0; j<size; j++) {
                    sites[j][i].deformedPosition.PE(0,-x0);
                }
                avg += sites[size-1][i].deformedPosition.component(0);
            }
            avg /= size;
            for(int i=0; i<size; i++) {
                double factorX = avg/sites[size-1][i].deformedPosition.component(0);
                for(int j=0; j<size; j++) {
                    Space.Vector position = sites[j][i].deformedPosition;
                    position.setComponent(0,factorX*position.component(0));//scale each row to uniform width
                }
            }
            for(int i=0; i<size; i++) {
                double y0 = sites[i][0].deformedPosition.component(1);
                double factorY = 1.0/(sites[i][size-1].deformedPosition.component(1)-y0);
                for(int j=0; j<size; j++) {
                    Space.Vector position = sites[i][j].deformedPosition;
                    position.PE(1,-y0);
                    position.setComponent(1,factorY*position.component(1));
                }
            }
         //   scale = sites[size-1][0].deformedPosition.component(0) - sites[0][0].deformedPosition.component(0);
            scale = avg;
        }//end of MyLattice constructor
        
        public MyCell getDeformedCell(Space.Vector r) {
            //find and return cell containing r in deformed lattice
            return null;
        }
        public MyCell getOriginalCell(Space2D.Vector r) {
            //find and return cell containing r in undeformed lattice
            int ix = (int)r.x * size;
            int iy = (int)r.y * size;
            MySite site = sites[ix][iy];
            double dx = r.x - site.originalPosition.x;
            double dy = r.x - site.originalPosition.y;
            return (dy > dx) ? site.c1 : site.c2;
        }
        
        private void deform(MySite site, OdeSolver.Rhs rhs) {
            //deformation algorithm
        }
        
        //AbstractLattice methods
        public int D() {return 2;}
        public int siteCount() {return 2*squareLattice.siteCount();}
        public Site site(AbstractLattice.Coordinate coord) {return null;} //not used
        public Site randomSite() {return iterator.randomSite();}
        public SiteIterator iterator() {return iterator;}     //iterator for all sites in lattice
    }
    
    private class MySite extends Site {
        MyCell c1, c2;
        Space2D.Vector originalPosition, deformedPosition;
        MySite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
            super(parent, iterator, coord);
        }
    }
    private class MySiteFactory implements SiteFactory {
       public Site makeSite(AbstractLattice parent, 
                            SiteIterator.Neighbor iterator, 
                            AbstractLattice.Coordinate coord) {
           return new MySite(parent, iterator, coord);
       }
    }
    private class MyCell {
        public MySite origin;
        public double jacobian;
        public double a11, a12, a21, a22;
        public double deltaX, deltaY;
        MyCell(MySite s0, MySite s1, MySite s2) {
            origin = s0;
            deltaX = origin.deformedPosition.x - origin.originalPosition.x;
            deltaY = origin.deformedPosition.y - origin.deformedPosition.y;
            double xbn = s1.deformedPosition.x - s0.deformedPosition.x;
            double xan = s2.deformedPosition.x - s0.deformedPosition.x;
            double ybn = s1.deformedPosition.y - s0.deformedPosition.y;
            double yan = s2.deformedPosition.y - s0.deformedPosition.y;
            double xb  = s1.originalPosition.x - s0.originalPosition.x;
            double xa  = s2.originalPosition.x - s0.originalPosition.x;
            double yb  = s1.originalPosition.y - s0.originalPosition.y;
            double ya  = s2.originalPosition.y - s0.originalPosition.y;
            double denominator = ya*xb - yb*xa;
            a11 = (xbn*ya - xan*yb)/denominator;
            a12 = (xan*xb - xa*xbn)/denominator;
            a21 = (ya*ybn - yan*yb)/denominator;
            a22 = (yan*xb - xa*ybn)/denominator;
            jacobian = a11*a22 - a12*a21;
        }
        
        public void transform(Space2D.Vector r, double scale) {
            double dx = r.x - origin.originalPosition.x;
            double dy = r.y - origin.originalPosition.y;
            double dxn = a11*dx + a12*dy;
            double dyn = a21*dx + a22*dy;
            r.x = origin.originalPosition.x + scale*(deltaX + dxn);
            r.y = origin.originalPosition.y + scale*(deltaY + dyn);
        }
    }
}

class VelocityField {
    
    double[][] source;
    double[] velocity;
    int D;
    int nSource;
    double m = 1;
    VelocityField(int dim, int nSourceHalf) {
        D = dim;
        velocity = new double[D];
        nSource = 2*nSourceHalf + 1;
        source = new double[nSource][D];
        int k=0;
        for(int i=-nSourceHalf; i<=nSourceHalf; i++) {
            for(int j=0; j<D-1; j++) {
                source[k][j] = 0.5;
            }
            source[k][D-1] = 0.5+(double)i;
            k++;
        }
    }
        //written for 2D
    public double[] v(double[] r) {
        velocity[0] = 0.0;
        velocity[1] = 0.0;
        for(int k=0; k<nSource; k++) {
            double dx = r[0] - source[k][0];
            double dy = r[1] - source[k][1];
            double r2 = dx*dx + dy*dy;
            velocity[0] += dx/r2;
            velocity[1] += dy/r2;
        }
        velocity[0] *= m;
        velocity[1] *= m;
        return velocity;
    }
}//end of VelocityField