package etomica.modules.crystalviewer;

import javax.swing.JTabbedPane;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.api.IVector;

import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicDiamond;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.LatticeHcp;
import etomica.lattice.LatticePlane;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;


public class CrystalViewer extends SimulationPanel {
    
	final static String APP_NAME = "Crystal Viewer";
    protected final ISimulation sim;

    protected SpeciesSpheresMono species;
    protected IBox box;
    protected IVector center;
    protected LatticePlane latticePlane;
    protected ClipPlaneEditor clipPlaneEditor;
    protected DisplayBox displayBox;
    protected LatticeEditor latticeEditor;


    public CrystalViewer() {
        super(APP_NAME);
        Space space = Space3D.getInstance();
        sim = new Simulation(space);
        center = space.makeVector();

        species = new SpeciesSpheresMono(sim, space);
        sim.getSpeciesManager().addSpecies(species);

        BasisMonatomic basisMonatomic = new BasisMonatomic(space);
        
        BravaisLattice[] lattices = new BravaisLattice[] {
                new LatticeCubicSimple(),
                new BravaisLatticeCrystal(new PrimitiveTetragonal(space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveHexagonal(space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveOrthorhombic(space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveMonoclinic(space), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveTriclinic(space), basisMonatomic),
                new LatticeCubicFcc(),
                new LatticeCubicBcc(),
                new LatticeHcp(),
                new LatticeCubicDiamond()
            };

        double[]  boxSize = new double[] { 10.0, 10.0, 10.0 };
        
        box  = new Box(new BoundaryDeformableLattice(lattices[0].getPrimitive(),
        		                                         (etomica.api.IRandom)null,
        		                                         boxSize), space);
        sim.addBox(box);

        String[] latticeNames = new String[]{
                "Simple Cubic", "Tetragonal", "Hexagonal", "Orthorhombic", "Monoclinic", "Triclinic", "FCC", "BCC", "HCP", "Diamond"};

        displayBox = new DisplayBox(box, space);
        displayBox.setPixelUnit(new Pixel(20));
        displayBox.setResizeOnNewBox(false);

        // we pass these to make LatticePlane happy.  they'll get whacked by update() later
        latticePlane = new LatticePlane(lattices[0].getPrimitive(), new int[] {1,0,0});
        
        clipPlaneEditor = new ClipPlaneEditor(latticePlane, displayBox);
        
        latticeEditor = new LatticeEditor(this, lattices, latticeNames, space);
        
        JTabbedPane controlTabs = new JTabbedPane();
        controlTabs.add("Crystal", latticeEditor.getPanel());
        controlTabs.add("Plane", clipPlaneEditor.getPanel());

        controlPanel.add(controlTabs);
        graphicsPanel.add(displayBox.graphic());
        toolbar.addContributor("Colin Tedlock");
    }

    public void update(BravaisLattice currentLattice) {
        latticePlane.setPrimitive(currentLattice.getPrimitive());
        displayBox.setLabel(currentLattice.toString());
        clipPlaneEditor.update();
        displayBox.repaint();
    }    
    
    public static class Applet extends javax.swing.JApplet {

	    public void init() {
            CrystalViewer viewer = new CrystalViewer();
		    getContentPane().add(viewer);
	    }
    }
    
 
    public static void main(String[] args) {
        
        SimulationPanel viewer = new CrystalViewer();
        SimulationGraphic.makeAndDisplayFrame(viewer, APP_NAME);
    }
    
}