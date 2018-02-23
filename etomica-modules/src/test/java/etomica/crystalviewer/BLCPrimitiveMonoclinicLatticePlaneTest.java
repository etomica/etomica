/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.crystalviewer;

import etomica.space.Vector;
import junit.framework.TestCase;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.lattice.crystal.PrimitiveMonoclinic;

public class BLCPrimitiveMonoclinicLatticePlaneTest extends TestCase {

	private final int DEFAULT_SIZE = 7;
	private final int DEFAULT_MILLER[] = {1,0,0};
	private final double DEFAULT_BETA = Math.PI * 2 * (90.0 / 360.0);
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};
	
	private String funcName = "";

	private double epsilon = 1.0E-5;;

	private LatticePlaneTestUtility lptu = null;

	public BLCPrimitiveMonoclinicLatticePlaneTest(String name) {
		super(name);
		funcName = name;
	}

	protected void setUp() throws Exception {
		super.setUp();
		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndBox(lptu.MONOCLINIC, DEFAULT_MILLER, DEFAULT_BOX);
	        lptu.setDimensions(DEFAULT_SIZE);
		}
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	private double[] makeArray(Vector v) {
	    return new double[] {v.getX(0), v.getX(1), v.getX(2)};
	}

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A,B,C) = 1.0
     * cells per side = 7
     * plane = 1.0
     * beta = 90 degrees
     */
    public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 1.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.MONOCLINIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(DEFAULT_BETA);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(0) >= spacePos-epsilon &&
                   a.getPosition().getX(0) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(0) >= spacePos-epsilon &&
               a.getPosition().getX(0) <= spacePos+epsilon) {
                System.out.println(funcName + " -> Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            	fail();
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is."); 
            	fail();
            }
        }

    } // End testStandard()

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A) = 1.5
     * size of cell (B) = 1.75
     * size of cell (C) = 2.0
     * cells per side = 7
     * plane = -2.0
     * beta = 90 degrees
     */
    public void testCellSizeAllDifferent() {

    	int idx = 0;
    	double cubicSizeA = 1.5;
    	double cubicSizeB = 1.75;
    	double cubicSizeC = 2.0;
    	double plane = -2.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.MONOCLINIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(DEFAULT_BETA);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);

                if(a.getPosition().getX(0) >= spacePos-epsilon &&
                   a.getPosition().getX(0) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(0) >= spacePos-epsilon &&
               a.getPosition().getX(0) <= spacePos+epsilon) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            	fail();
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            	fail();
            }
        }
    } // End testCellSizeAllDifferent()

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A) = 1.25
     * size of cell (B) = 1.5
     * size of cell (C) = 1.75
     * cells per side = 7
     * plane = 3.0
     * beta = 125 degrees
     */
    public void testCellSizeAllDifferentBeta125() {

    	int idx = 0;
    	double cubicSizeA = 1.25;
    	double cubicSizeB = 1.5;
    	double cubicSizeC = 1.75;
    	double plane = 3.0;
    	int itemsFound = 0;
    	double beta = Math.PI * 2.0 * (125.0 / 360.0);
    	IAtomList leafList = null;
        double actualPlane[][] =
    	{ { 6.761276290842992, -4.5, -4.300548232517207 }, { 5.757517527228662, -4.5, -2.8670321550114717},
      	  { 4.75375876361433, -4.5, -1.433516077505736 }, { 3.749999999999999, -4.5, 0.0},
    	  { 2.7462412363856687, -4.5, 1.4335160775057352 }, { 1.7424824727713384, -4.5, 2.8670321550114712},
    	  { 0.7387237091570071, -4.5, 4.300548232517207 }, { 6.761276290842992, -3.0, -4.300548232517207},
    	  { 5.757517527228662, -3.0, -2.8670321550114717 }, { 4.75375876361433, -3.0, -1.433516077505736},
    	  { 3.749999999999999, -3.0, 0.0 }, { 2.7462412363856687, -3.0, 1.4335160775057352},
    	  { 1.7424824727713384, -3.0, 2.8670321550114712 }, { 0.7387237091570071, -3.0, 4.300548232517207},
    	  { 6.761276290842992, -1.5, -4.300548232517207 }, { 5.757517527228662, -1.5, -2.8670321550114717},
    	  { 4.75375876361433, -1.5, -1.433516077505736 }, { 3.749999999999999, -1.5, 0.0},
    	  { 2.7462412363856687, -1.5, 1.4335160775057352 }, { 1.7424824727713384, -1.5, 2.8670321550114712},
    	  { 0.7387237091570071, -1.5, 4.300548232517207 }, { 6.761276290842992, 0.0, -4.300548232517207},
    	  { 5.757517527228662, 0.0, -2.8670321550114717 }, { 4.75375876361433, 0.0, -1.433516077505736},
    	  { 3.749999999999999, 0.0, 0.0 }, { 2.7462412363856687, 0.0, 1.4335160775057352},
    	  { 1.7424824727713384, 0.0, 2.8670321550114712 }, { 0.7387237091570071, 0.0, 4.300548232517207},
    	  { 6.761276290842992, 1.5, -4.300548232517207 }, { 5.757517527228662, 1.5, -2.8670321550114717},
    	  { 4.75375876361433, 1.5, -1.433516077505736 }, { 3.749999999999999, 1.5, 0.0},
    	  { 2.7462412363856687, 1.5, 1.4335160775057352 }, { 1.7424824727713384, 1.5, 2.8670321550114712},
    	  { 0.7387237091570071, 1.5, 4.300548232517207 }, { 6.761276290842992, 3.0, -4.300548232517207},
    	  { 5.757517527228662, 3.0, -2.8670321550114717 }, { 4.75375876361433, 3.0, -1.433516077505736},
    	  { 3.749999999999999, 3.0, 0.0 }, { 2.7462412363856687, 3.0, 1.4335160775057352},
    	  { 1.7424824727713384, 3.0, 2.8670321550114712 }, { 0.7387237091570071, 3.0, 4.300548232517207},
    	  { 6.761276290842992, 4.5, -4.300548232517207 }, { 5.757517527228662, 4.5, -2.8670321550114717},
    	  { 4.75375876361433, 4.5, -1.433516077505736 }, { 3.749999999999999, 4.5, 0.0},
    	  { 2.7462412363856687, 4.5, 1.4335160775057352 }, { 1.7424824727713384, 4.5, 2.8670321550114712},
    	  { 0.7387237091570071, 4.5, 4.300548232517207} };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.MONOCLINIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
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

    } // End testCellSizeAllDifferentBeta125()

    /*
     * Miller indices = 3, 2, 1
     * size of cell (A,B,C) = 1.0
     * cells per side = 5
     * plane = 1.0
     * beta = 147 degrees
     */
    public void testOddMillerIndicesBeta147() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	IAtomList leafList = null;
    	int size = 5;
    	double plane = 1.0;
    	int itemsFound = 0;
    	double beta = Math.PI * 2.0 * (147.0 / 360.0);
    	int[] millerIndices = new int[] { 3, 2, 1 };
        double actualPlane[][] =
           { { -2.677341135890848, 1.0, 1.0892780700300546 },
             { -1.0, 2.0, 0.0 },
             { -0.838670567945424, 0.0, 0.5446390350150272 },
             { 0.8386705679454236, 1.0, -0.5446390350150273 },
             { -0.6773411358908481, -2.0, 1.0892780700300544 },
             { 0.9999999999999996, -1.0, 0.0 },
             { 2.6773411358908477, 0.0, -1.0892780700300546 },
             { 2.8386705679454236, -2.0, -0.5446390350150272 } };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.MONOCLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
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

    } // End testOddMillerIndicesBeta147()

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A,B,C) = 1.0
     * cells per side = 4
     * plane = 0
     * beta = 90 degrees
     */
    public void testEvenAtomsPerSideZeroPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	IAtomList leafList = null;
    	int dimensionSize = 4;
    	double plane = 0.0;

        lptu.createLatticeAndBox(lptu.MONOCLINIC, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(DEFAULT_BETA);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
            	assertFalse(lptu.getLatticePlane().inPlane(
            	    	a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            fail();
        }
    	
    } // End testEvenAtomsPerSideZeroPlane


   /*
    * Miller indices = 1, 0, 0
    * size of cell (A,B,C) = 1.0
    * cells per side = 4
    * plane = 1.5
    * beta = 115 degrees
    */
    public void testEvenAtomsPerSidePt5Plane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	IAtomList leafList = null;
    	int dimensionSize = 4;
    	int itemsFound = 0;
    	double beta = Math.PI * 2.0 * (115.0 / 360.0);
    	double plane = 1.5;
        double actualPlane[][] =
          { { 2.1339273926110494, -1.5, -1.359461680554975 }, { 1.71130913087035, -1.5, -0.453153893518325},
            { 1.2886908691296506, -1.5, 0.453153893518325 }, { 0.8660726073889511, -1.5, 1.359461680554975},
            { 2.1339273926110494, -0.5, -1.359461680554975 }, { 1.71130913087035, -0.5, -0.453153893518325},
            { 1.2886908691296508, -0.5, 0.453153893518325 }, { 0.8660726073889511, -0.5, 1.359461680554975},
            { 2.1339273926110494, 0.5, -1.359461680554975 }, { 1.7113091308703499, 0.5, -0.453153893518325},
            { 1.2886908691296506, 0.5, 0.453153893518325 }, { 0.8660726073889511, 0.5, 1.359461680554975},
            { 2.1339273926110494, 1.5, -1.359461680554975 }, { 1.7113091308703499, 1.5, -0.453153893518325},
            { 1.2886908691296506, 1.5, 0.453153893518325 }, { 0.8660726073889511, 1.5, 1.359461680554975} };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.MONOCLINIC, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
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
    	
    } // End testEvenAtomsPerSideZeroPt5Plane

    /*
     * Miller indices = 1, 1, 1
     * size of cell (A) = 1.1
     * size of cell (B) = 1.2
     * size of cell (C) = 1.3
     * cells per side = 7
     * plane = 2.95
     * beta = 155 degrees
     */
    public void testPlaneMinusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 2.95;
    	double beta = Math.PI * 2.0 * (155.0 / 360.0);
    	int[] millerIndices = new int[] { 1, 1, 1 };

        lptu.createLatticeAndBox(lptu.MONOCLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            fail();
        }

    } // End testPlaneMinusFiveHundreths()

    /*
     * Miller indices = 1, 1, 1
     * size of cell (A) = 1.1
     * size of cell (B) = 1.2
     * size of cell (C) = 1.3
     * cells per side = 7
     * plane = 3.0
     * beta = 155 degrees
     */
    public void testPlane() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 3.0;
    	int itemsFound = 0;
    	double beta = Math.PI * 2.0 * (155.0 / 360.0);
    	int[] millerIndices = new int[] { 1, 1, 1 };
        double actualPlane[][] =
          { { -6.834600369442935, 3.6, 1.648211220788728 }, { -5.734600369442934, 2.4, 1.648211220788728 },
        	{ -4.556400246295291, 3.6, 1.0988074805258188 }, { -4.634600369442935, 1.2, 1.648211220788728},
            { -3.4564002462952903, 2.4, 1.0988074805258188 }, { -2.278200123147645, 3.6, 0.5494037402629093},
            { -3.5346003694429347, 0.0, 1.648211220788728 }, { -2.35640024629529, 1.2, 1.0988074805258188},
            { -1.1782001231476449, 2.4, 0.5494037402629093 }, { 0.0, 3.6000000000000005, 0.0},
            { -2.4346003694429346, -1.2, 1.648211220788728 }, { -1.2564002462952901, 0.0, 1.0988074805258188},
            { -0.0782001231476448, 1.2, 0.5494037402629093 }, { 1.1, 2.4000000000000012, 0.0},
            { 2.278200123147645, 3.6, -0.5494037402629093 }, { -1.334600369442935, -2.4, 1.648211220788728},
            { -0.15640024629529048, -1.2, 1.0988074805258188 }, { 1.0217998768523548, 0.0, 0.5494037402629093},
            { 2.1999999999999997, 1.2, 0.0 }, { 3.3782001231476446, 2.4, -0.5494037402629093},
            { 4.556400246295289, 3.6, -1.0988074805258186 }, { -0.2346003694429344, -3.6, 1.648211220788728},
            { 0.94359975370471, -2.4, 1.0988074805258188 }, { 2.1217998768523554, -1.2, 0.5494037402629093},
            { 3.3000000000000003, 0.0, 0.0 }, { 4.478200123147645, 1.2, -0.5494037402629093},
            { 5.6564002462952905, 2.4, -1.0988074805258186 }, { 6.834600369442935, 3.6, -1.648211220788728} };
        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.MONOCLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
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
     * Miller indices = 1, 1, 1
     * size of cell (A) = 1.1
     * size of cell (B) = 1.2
     * size of cell (C) = 1.3
     * cells per side = 7
     * plane = 3.05
     * beta = 155 degrees
     */
    public void testPlanePlusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 3.05;
    	double beta = Math.PI * 2.0 * (155.0 / 360.0);
    	int[] millerIndices = new int[] { 1, 1, 1 };

        lptu.createLatticeAndBox(lptu.MONOCLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveMonoclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
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
