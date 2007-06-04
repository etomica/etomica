package etomica.modules.junit.crystalviewer;

import junit.framework.TestCase;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtomPositioned;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.space.IVector;


public class BLCPrimitiveOrthorhombicLatticePlaneTest extends TestCase {

	private final int DEFAULT_SIZE = 5;
	private final int DEFAULT_MILLER[] = {0,1,0};
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};
	
	private String funcName = "";

	private double epsilon = 1.0E-5;;

	private LatticePlaneTestUtility lptu = null;

	public BLCPrimitiveOrthorhombicLatticePlaneTest(String name) {
		super(name);
		funcName = name;
	}

	protected void setUp() throws Exception {
		super.setUp();
		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, DEFAULT_MILLER, DEFAULT_BOX);
	        lptu.setDimensions(DEFAULT_SIZE);
		}
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	private double[] makeArray(IVector v) {
	    return new double[] {v.x(0), v.x(1), v.x(2)};
	}

    /*
     * Miller indices = 0, 1, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 5
     * plane = 1.0
     */
    public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 1.0;
    	AtomArrayList leafList = null;

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
                if(a.getPosition().x(1) >= spacePos-epsilon &&
                   a.getPosition().x(1) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(a.getPosition().x(1) >= spacePos-epsilon &&
               a.getPosition().x(1) <= spacePos+epsilon) {
            	System.out.println(funcName + " -> Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }

    } // End testStandard()

    /*
     * Miller indices = 0, 1, 0
     * size of cell (A) = 1.5
     * size of cell (B) = 1.75
     * size of cell (C) = 2.0
     * cells per side = 5
     * plane = -2.0
     */
    public void testCellSizeAllDifferentIncrease() {

    	int idx = 0;
    	double cubicSizeA = 1.5;
    	double cubicSizeB = 1.75;
    	double cubicSizeC = 2.0;
    	double plane = -2.0;
    	AtomArrayList leafList = null;

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);

        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
                if(a.getPosition().x(1) >= spacePos-epsilon &&
                   a.getPosition().x(1) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(a.getPosition().x(1) >= spacePos-epsilon &&
               a.getPosition().x(1) <= spacePos+epsilon) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }
    } // End testABCellSizeIncrease()

    /*
     * Miller indices = 1, 3, 2
     * size of cell (A, B, C) = 1.0
     * cells per side = 9
     * plane = 8.0
     */
    public void testOddMillerIndicesDistantPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	AtomArrayList leafList = null;
    	int size = 9;
    	double plane = 8.0;
    	int itemsFound = 0;
    	int[] millerIndices = new int[] { 1, 3, 2 };
        double actualPlane[][] =
                { { -4.0, 2.0, 3.0 }, { -4.0, 4.0, 0.0 }, { -3.0, 1.0, 4.0},
                  { -3.0, 3.0, 1.0 }, { -2.0, 2.0, 2.0 }, { -2.0, 4.0, -1.0},
                  { -1.0, 1.0, 3.0 }, { -1.0, 3.0, 0.0 }, { 0.0, 0.0, 4.0},
                  { 0.0, 2.0, 1.0 }, { 0.0, 4.0, -2.0 }, { 1.0, 1.0, 2.0},
                  { 1.0, 3.0, -1.0 }, { 2.0, 0.0, 3.0 }, { 2.0, 2.0, 0.0},
                  { 2.0, 4.0, -3.0 }, { 3.0, -1.0, 4.0 }, { 3.0, 1.0, 1.0},
                  { 3.0, 3.0, -2.0 }, { 4.0, 0.0, 2.0 }, { 4.0, 2.0, -1.0},
                  { 4.0, 4.0, -4.0} };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
	            	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(dd.contains(makeArray(a.getPosition()))) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }

        assertEquals(actualPlane.length, itemsFound);

    } // End testOddMillerIndicesDistantPlane()

    /*
     * Miller indices = 0, 1, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 4
     * plane = 0.0
     */
    public void testEvenAtomsPerSideZeroPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	AtomArrayList leafList = null;
    	int dimensionSize = 4;
    	double plane = 0.0;

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            	assertFalse(lptu.getLatticePlane().inPlane(
            	    	(etomica.space3d.Vector3D)(a.getPosition())));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
             fail();
        }
    	
    } // End testEvenAtomsPerSideZeroPlane

    /*
     * Miller indices = 0, 1, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 4
     * plane = 0.5
     */
    public void testEvenAtomsPerSideZeroPt5Plane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	AtomArrayList leafList = null;
    	int dimensionSize = 4;
    	double plane = 0.5;

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

                if(a.getPosition().x(1) >= spacePos-epsilon &&
                   a.getPosition().x(1) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(a.getPosition().x(1) >= spacePos-epsilon &&
               a.getPosition().x(1) <= spacePos+epsilon) {
                System.out.println(funcName + " -> Atom position : " + a.getPosition() +
                 			" should be in plane but is not.");
             }
             else {
                 System.out.println(funcName + " ->Atom position : " + a.getPosition() +
                 			" should not be in plane but is.");       	
             }
         	fail();
        }
    	
    } // End testEvenAtomsPerSideZeroPt5Plane

    /*
     * Miller indices = 2, 2, 1
     * size of cell (A) = 0.7
     * size of cell (B) = 0.8
     * size of cell (C) = 0.9
     * cells per side = 7
     * plane = 2.95
     */
    public void testPlaneMinusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 0.7;
    	double cubicSizeB = 0.8;
    	double cubicSizeC = 0.9;
    	AtomArrayList leafList = null;
    	int size = 7;
    	double plane = 2.95;
    	int[] millerIndices = new int[] { 2, 2, 1 };

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
         	fail();
        }


    } // End testPlaneMinusFiveHundreths()

    /*
     * Miller indices = 2, 2, 1
     * size of cell (A) = 0.7
     * size of cell (B) = 0.8
     * size of cell (C) = 0.9
     * cells per side = 7
     * plane = 3.0
     */
    public void testPlane() {

    	int idx = 0;
    	double cubicSizeA = 0.7;
    	double cubicSizeB = 0.8;
    	double cubicSizeC = 0.9;
    	AtomArrayList leafList = null;
    	int size = 7;
    	double plane = 3.0;
    	int itemsFound = 0;
    	int[] millerIndices = new int[] { 2, 2, 1 };
        double actualPlane[][] =
                  { { -2.1, 2.4, 2.7 }, { -1.4, 1.6, 2.7 }, { -1.4, 2.4, 0.9 },
                    { -0.7, 0.8, 2.7 }, { -0.7, 1.6, 0.9 }, { -0.7, 2.4, -0.9 },
                    { 0.0, 0.0, 2.7}, { 0.0, 0.8, 0.9 }, { 0.0, 1.6, -0.9 },
                    { 0.0, 2.4, -2.7 }, { 0.7, -0.8, 2.7 }, { 0.7, 0.0, 0.9 },
                    { 0.7, 0.8, -0.9 }, { 0.7, 1.6, -2.7 }, { 1.4, -1.6, 2.7 },
                    { 1.4, -0.8, 0.9 }, { 1.4, 0.0, -0.9 }, { 1.4, 0.8, -2.7 },
                    { 2.1, -2.4, 2.7 }, { 2.1, -1.6, 0.9 }, { 2.1, -0.8, -0.9 },
                    { 2.1, 0.0, -2.7} };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
	            	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(dd.contains(makeArray(a.getPosition()))) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }

        assertEquals(actualPlane.length, itemsFound);

    } // End testPlane()

    /*
     * Miller indices = 2, 2, 1
     * size of cell (A) = 0.7
     * size of cell (B) = 0.8
     * size of cell (C) = 0.9
     * cells per side = 7
     * plane = 3.05
     */
    public void testPlanePlusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 0.7;
    	double cubicSizeB = 0.8;
    	double cubicSizeC = 0.9;
    	AtomArrayList leafList = null;
    	int size = 7;
    	double plane = 3.05;
    	int[] millerIndices = new int[] { 2, 2, 1 };

        lptu.createLatticeAndPhase(lptu.ORTHORHOMBIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveOrthorhombic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getPhase().getSpeciesMaster().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
         	fail();
        }


    } // End testPlanePlusFiveHundreths()

    public class DoubleTwoDArray {
    	private double[][] array;
    	private double epsilon = 1.0E-5;

    	public DoubleTwoDArray(double[][] in) {
    		array = in;
    	}
    	
    	public boolean contains(double[] val) {
    		boolean b = false;

    		if(array[0].length == val.length) {
    			for(int i = 0; i < array.length; i++) {
    				boolean matching = true;
    				for(int j = 0; j < array[0].length; j++) {
    					if(val[j] < array[i][j] - epsilon ||
    					   val[j] > array[i][j] + epsilon) {
    						matching = false;
    						break;
    					}
    				}
    				if(matching == true) {
    					b = true;
    					break;
    				}
    			}
    		}
            return b;
    	}
    	
    	public void setEpsilon(double e) {
    		this.epsilon = e;
    	}

    }  // End public class DoubleTwoDArray

}
