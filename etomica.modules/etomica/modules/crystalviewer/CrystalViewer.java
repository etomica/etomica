package etomica.modules.crystalviewer;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

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
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;


public class CrystalViewer extends SimulationPanel {
    
	final static String APP_NAME = "Crystal Viewer";
    protected final ISimulation sim;
    private JPanel mainPanel;

    protected SpeciesSpheresMono species;
    protected Box box;
    protected IVector center;
    protected LatticePlane latticePlane;
    protected ClipPlaneEditor clipPlaneEditor;
    protected DisplayBox displayBox;
    protected LatticeEditor latticeEditor;


    public CrystalViewer() {
        super(APP_NAME);
        sim = new Simulation(Space3D.getInstance());
        center = sim.getSpace().makeVector();

        species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);

        BasisMonatomic basisMonatomic = new BasisMonatomic(sim.getSpace());
        
        BravaisLattice[] lattices = new BravaisLattice[] {
                new LatticeCubicSimple(),
                new BravaisLatticeCrystal(new PrimitiveTetragonal(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveHexagonal(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveOrthorhombic(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveMonoclinic(sim.getSpace()), basisMonatomic),
                new BravaisLatticeCrystal(new PrimitiveTriclinic(sim.getSpace()), basisMonatomic),
                new LatticeCubicFcc(),
                new LatticeCubicBcc(),
                new LatticeHcp(),
                new LatticeCubicDiamond()
            };

        double[]  boxSize = new double[] { 10.0, 10.0, 10.0 };
        
        box  = new Box(new BoundaryDeformableLattice(lattices[0].getPrimitive(),
        		                                         (etomica.util.IRandom)null,
        		                                         boxSize));
        sim.addBox(box);

        String[] latticeNames = new String[]{
                "Simple Cubic", "Tetragonal", "Hexagonal", "Orthorhombic", "Monoclinic", "Triclinic", "FCC", "BCC", "HCP", "Diamond"};

        displayBox = new DisplayBox(box);
        displayBox.setPixelUnit(new Pixel(20));
        displayBox.setResizeOnNewBox(false);

        // we pass these to make LatticePlane happy.  they'll get whacked by update() later
        latticePlane = new LatticePlane(lattices[0].getPrimitive(), new int[] {1,0,0});
        
        clipPlaneEditor = new ClipPlaneEditor(latticePlane, displayBox);
        
        latticeEditor = new LatticeEditor(this, lattices, latticeNames);
        
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