/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;

/**
 * A polyhedron with six faces, all of which are four-sided polygons.
 * Other six-face polyhedra (e.g., trigonal bipyrimid) cannot
 * be derived from this class.  This polytope has 12 edges and 8 vertices.
 *  
 * @author kofke
 *  
 */
public abstract class Hexahedron extends Polyhedron {

    /**
     * Constructs a hexahedron to reside in the given space, which
     * determines the type of vertices it will have (e.g. based on Vector3D).
     * Faces are formed from the following sets of vertices (referring to
     * the index of the vertices array): (0,1,2,3), (1,3,7,5), (4,5,7,6),
     * (0,4,6,2), (2,3,7,6), and (0,1,5,4) 
     */
    public Hexahedron(Space embeddedSpace) {
        super(makeFaces(embeddedSpace));
    }
    
    private static Polygon[] makeFaces(Space embeddedSpace) {
        Vector[] vertices = embeddedSpace.makeVectorArray(8);
        LineSegment[] edges = new LineSegment[12];
        Polygon[] faces = new Polygon[6];
        //note that it is important that the vertices get used in order, 0, 1, 2, etc.
        for(int i=0; i<vertices.length; i++) {vertices[i].setX(0,i);}
        edges[0] = new LineSegment(embeddedSpace, vertices[0],vertices[1]);
        edges[1] = new LineSegment(embeddedSpace, vertices[2],vertices[3]);
        edges[2] = new LineSegment(embeddedSpace, vertices[1],vertices[3]);
        edges[3] = new LineSegment(embeddedSpace, vertices[0],vertices[2]);
        edges[4] = new LineSegment(embeddedSpace, vertices[0],vertices[4]);
        edges[5] = new LineSegment(embeddedSpace, vertices[4],vertices[5]);
        edges[6] = new LineSegment(embeddedSpace, vertices[1],vertices[5]);
        edges[7] = new LineSegment(embeddedSpace, vertices[4],vertices[6]);
        edges[8] = new LineSegment(embeddedSpace, vertices[6],vertices[7]);
        edges[9] = new LineSegment(embeddedSpace, vertices[5],vertices[7]);
        edges[10] = new LineSegment(embeddedSpace, vertices[2],vertices[6]);
        edges[11] = new LineSegment(embeddedSpace, vertices[3],vertices[7]);
        //note that it is important that the edges get used in order, 0, 1, 2, etc.
        faces[0] = new PolygonGeneral(new LineSegment[] {edges[0],edges[1],edges[2],edges[3]});
        faces[1] = new PolygonGeneral(new LineSegment[] {edges[0],edges[4],edges[5],edges[6]});
        faces[2] = new PolygonGeneral(new LineSegment[] {edges[5],edges[7],edges[8],edges[9]});
        faces[3] = new PolygonGeneral(new LineSegment[] {edges[1],edges[8],edges[10],edges[11]});
        faces[4] = new PolygonGeneral(new LineSegment[] {edges[3],edges[4],edges[7],edges[10]});
        faces[5] = new PolygonGeneral(new LineSegment[] {edges[2],edges[6],edges[9],edges[11]});
        return faces;
    }

}
