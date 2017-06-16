/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on June 2, 2005
 */
package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * 
 * A truncated octahedron is an Archimedean, space-filling polyhedron with 24
 * vertices, 36 edges, and 14 faces. All edges are of the same length. Eight of
 * the faces are regular hexagons, and six are squares.
 * 
 * @author Katherine Schadel
 */
public class TruncatedOctahedron extends Polyhedron {

    /**
     * Constructs a truncated octahedron with vertices defined in the embedded
     * space and with unit edge lengths.
     */
    public TruncatedOctahedron(Space embeddedSpace) {
        this(embeddedSpace, 1.0);
    }

    /**
     * Constructs a truncated octahedron with vertices defined in the embedded
     * space and with edge lengths equal to the given value.
     */

    public TruncatedOctahedron(Space embeddedSpace, double edgeLength) {
        super(makeFaces(embeddedSpace));
        setEdgeLength(edgeLength);
    }

    /**
     * @return the common length of all edges of the truncated octahedron
     */

    public double getEdgeLength() {
        return edgeLength;
    }

    /**
     * Specifies the common length of all edges of the truncated octahedron.
     * 
     * @param edgeLength the new length of all edges
     */

    public void setEdgeLength(double edgeLength) {
        this.edgeLength = edgeLength;
        updateVertices();
    }

    /**
     * Returns the length of an edge of a cube in which the truncated octahedron
     * is inscribed.  This is the distance from one of the square faces of the truncated
     * octahedron to the square face opposite it.  It is equal to 2*sqrt(2)*edgeLength.
     */
    public double getContainingCubeEdgeLength() {
        return 2 * Math.sqrt(2) * edgeLength;
    }

    /**
     * Sets the size of the truncated octahedron in terms of the edge length of
     * a cube that would contain it.
     */
    public void setContainingCubeEdgeLength(double length) {
        edgeLength = length / (2 * Math.sqrt(2));
    }

    /**
     * The following defines the vertices that the truncated octahedron will
     * have, based on Vector3D. Edges are formed from the set of vertices by
     * referring to the index of vertices[], and faces are formed by referring
     * to the index of the edges[]. Thus it is essential that when constructing
     * edges[], the first reference to each vertex must be sequential. For
     * example, the first two edges might look something like this:
     * vertices[0],vertices[1] and vertices[0], vertices[2]. Once the point
     * vertices[0] has appeared for the first time, it may appear at any
     * location in edges[] thereafter. However, the first two edges could not be
     * vertices[0],vertices[2] and vertices[0],vertices[1], because vertices[2]
     * was called before vertices[1] was called for the first time. This same
     * logic applies to the ordering of the edges when defining the faces.
     * 
     * The variable "n" has been incorporated into the definitions of the
     * vertices so that the size of the truncated octahedron can be adjusted.
     * 
     * The vertices have been numbered such that the faces are defined in terms
     * of vertices as follows: (0,1,2,3), (0,4,5,6,7,8), (1,8,9,10,11,12),
     * (2,12,13,14,15,16), (3,4,16,17,18,19), (5,19,20,21), (7,9,22,23),
     * (11,13,24,25), (15,17,26,27), (6,21,23,28,29,30), (10,22,25,28,31,32),
     * (14,24,27,32,33,34), (18,20,26,30,34,35), and (29,31,33,35).
     * 
     * The physical order of the edges on each face does not need to be
     * sequential. For example, when naming the edges of a square (let's call it
     * faces[4]), the clockwise order of the edges could be
     * edges[17],edges[5],edges[8],edges[3]. However, when listing the edges in
     * faces[4], the order could be edges[3], edges[5], edges[8], edges[17].
     * This greatly simplifies the task of numbering vertices.
     */

    public void updateVertices() {
        final double nSqrt2 = edgeLength * Math.sqrt(2);
        final double n2Sqrt2 = 0.5*nSqrt2;
        ((Vector3D) vertices[0]).E(n2Sqrt2, 0, nSqrt2);
        ((Vector3D) vertices[1]).E(0, -n2Sqrt2, nSqrt2);
        ((Vector3D) vertices[2]).E(-n2Sqrt2, 0, nSqrt2);
        ((Vector3D) vertices[3]).E(0, n2Sqrt2, nSqrt2);
        ((Vector3D) vertices[4]).E(nSqrt2, 0, n2Sqrt2);
        ((Vector3D) vertices[5]).E(nSqrt2, -n2Sqrt2, 0);
        ((Vector3D) vertices[6]).E(n2Sqrt2, -nSqrt2, 0);
        ((Vector3D) vertices[7]).E(0, -nSqrt2, n2Sqrt2);
        ((Vector3D) vertices[8]).E(-n2Sqrt2, -nSqrt2, 0);
        ((Vector3D) vertices[9]).E(-nSqrt2, -n2Sqrt2, 0);
        ((Vector3D) vertices[10]).E(-nSqrt2, 0, n2Sqrt2);
        ((Vector3D) vertices[11]).E(-nSqrt2, n2Sqrt2, 0);
        ((Vector3D) vertices[12]).E(-n2Sqrt2, nSqrt2, 0);
        ((Vector3D) vertices[13]).E(0, nSqrt2, n2Sqrt2);
        ((Vector3D) vertices[14]).E(n2Sqrt2, nSqrt2, 0);
        ((Vector3D) vertices[15]).E(nSqrt2, n2Sqrt2, 0);
        ((Vector3D) vertices[16]).E(nSqrt2, 0, -n2Sqrt2);
        ((Vector3D) vertices[17]).E(0, -nSqrt2, -n2Sqrt2);
        ((Vector3D) vertices[18]).E(-nSqrt2, 0, -n2Sqrt2);
        ((Vector3D) vertices[19]).E(0, nSqrt2, -n2Sqrt2);
        ((Vector3D) vertices[20]).E(0, -n2Sqrt2, -nSqrt2);
        ((Vector3D) vertices[21]).E(n2Sqrt2, 0, -nSqrt2);
        ((Vector3D) vertices[22]).E(-n2Sqrt2, 0, -nSqrt2);
        ((Vector3D) vertices[23]).E(0, n2Sqrt2, -nSqrt2);
        applyTranslationRotation();
    }

    /**
     * The for loop below creates the edges (and thus the faces) listed in the
     * arrays below in 3D space.
     * 
     * @param embeddedSpace
     * @return
     */

    private static Polygon[] makeFaces(Space embeddedSpace) {
        Vector[] vertices = embeddedSpace.makeVectorArray(24);
        LineSegment[] edges = new LineSegment[36];
        Polygon[] faces = new Polygon[14];
        for (int i = 0; i < vertices.length; i++) {
            vertices[i].setX(0, i);
        }
        //Remember, the first appearances of the vertices must be sequential.
        edges[0] = new LineSegment(embeddedSpace, vertices[0], vertices[1]);
        edges[1] = new LineSegment(embeddedSpace, vertices[1], vertices[2]);
        edges[2] = new LineSegment(embeddedSpace, vertices[2], vertices[3]);
        edges[3] = new LineSegment(embeddedSpace, vertices[0], vertices[3]);
        edges[4] = new LineSegment(embeddedSpace, vertices[0], vertices[4]);
        edges[5] = new LineSegment(embeddedSpace, vertices[4], vertices[5]);
        edges[6] = new LineSegment(embeddedSpace, vertices[5], vertices[6]);
        edges[7] = new LineSegment(embeddedSpace, vertices[6], vertices[7]);
        edges[8] = new LineSegment(embeddedSpace, vertices[1], vertices[7]);
        edges[9] = new LineSegment(embeddedSpace, vertices[7], vertices[8]);
        edges[10] = new LineSegment(embeddedSpace, vertices[8], vertices[9]);
        edges[11] = new LineSegment(embeddedSpace, vertices[9], vertices[10]);
        edges[12] = new LineSegment(embeddedSpace, vertices[2], vertices[10]);
        edges[13] = new LineSegment(embeddedSpace, vertices[10], vertices[11]);
        edges[14] = new LineSegment(embeddedSpace, vertices[11], vertices[12]);
        edges[15] = new LineSegment(embeddedSpace, vertices[12], vertices[13]);
        edges[16] = new LineSegment(embeddedSpace, vertices[3], vertices[13]);
        edges[17] = new LineSegment(embeddedSpace, vertices[13], vertices[14]);
        edges[18] = new LineSegment(embeddedSpace, vertices[14], vertices[15]);
        edges[19] = new LineSegment(embeddedSpace, vertices[4], vertices[15]);
        edges[20] = new LineSegment(embeddedSpace, vertices[15], vertices[16]);
        edges[21] = new LineSegment(embeddedSpace, vertices[5], vertices[16]);
        edges[22] = new LineSegment(embeddedSpace, vertices[8], vertices[17]);
        edges[23] = new LineSegment(embeddedSpace, vertices[6], vertices[17]);
        edges[24] = new LineSegment(embeddedSpace, vertices[11], vertices[18]);
        edges[25] = new LineSegment(embeddedSpace, vertices[9], vertices[18]);
        edges[26] = new LineSegment(embeddedSpace, vertices[14], vertices[19]);
        edges[27] = new LineSegment(embeddedSpace, vertices[12], vertices[19]);
        edges[28] = new LineSegment(embeddedSpace, vertices[17], vertices[20]);
        edges[29] = new LineSegment(embeddedSpace, vertices[20], vertices[21]);
        edges[30] = new LineSegment(embeddedSpace, vertices[16], vertices[21]);
        edges[31] = new LineSegment(embeddedSpace, vertices[20], vertices[22]);
        edges[32] = new LineSegment(embeddedSpace, vertices[18], vertices[22]);
        edges[33] = new LineSegment(embeddedSpace, vertices[22], vertices[23]);
        edges[34] = new LineSegment(embeddedSpace, vertices[19], vertices[23]);
        edges[35] = new LineSegment(embeddedSpace, vertices[21], vertices[23]);
        //Remember, the first appearances of the edges must be sequential.
        faces[0] = new PolygonGeneral(new LineSegment[] { edges[0], edges[1],
                edges[2], edges[3] });
        faces[1] = new PolygonGeneral(new LineSegment[] { edges[0], edges[4],
                edges[5], edges[6], edges[7], edges[8] });
        faces[2] = new PolygonGeneral(new LineSegment[] { edges[1], edges[8],
                edges[9], edges[10], edges[11], edges[12] });
        faces[3] = new PolygonGeneral(new LineSegment[] { edges[2], edges[12],
                edges[13], edges[14], edges[15], edges[16] });
        faces[4] = new PolygonGeneral(new LineSegment[] { edges[3], edges[4],
                edges[16], edges[17], edges[18], edges[19] });
        faces[5] = new PolygonGeneral(new LineSegment[] { edges[5], edges[19],
                edges[20], edges[21] });
        faces[6] = new PolygonGeneral(new LineSegment[] { edges[7], edges[9],
                edges[22], edges[23] });
        faces[7] = new PolygonGeneral(new LineSegment[] { edges[11], edges[13],
                edges[24], edges[25] });
        faces[8] = new PolygonGeneral(new LineSegment[] { edges[15], edges[17],
                edges[26], edges[27] });
        faces[9] = new PolygonGeneral(new LineSegment[] { edges[6], edges[21],
                edges[23], edges[28], edges[29], edges[30] });
        faces[10] = new PolygonGeneral(new LineSegment[] { edges[10],
                edges[22], edges[25], edges[28], edges[31], edges[32] });
        faces[11] = new PolygonGeneral(new LineSegment[] { edges[14],
                edges[24], edges[27], edges[32], edges[33], edges[34] });
        faces[12] = new PolygonGeneral(new LineSegment[] { edges[18],
                edges[20], edges[26], edges[30], edges[34], edges[35] });
        faces[13] = new PolygonGeneral(new LineSegment[] { edges[29],
                edges[31], edges[33], edges[35] });
        return faces;
    }

    public boolean contains(Vector v) {
        double length = getContainingCubeEdgeLength();
        double x = Math.abs(v.getX(0) - position.getX(0)) / length;
        double y = Math.abs(v.getX(1) - position.getX(1)) / length;
        double z = Math.abs(v.getX(2) - position.getX(2)) / length;
        if (x > 0.5 || y > 0.5 || z > 0.5)
            return false;
        if (x + y + z > 0.75)
            return false;
        return true;
    }

    public double getVolume() {
        updateVertices();
        double volume = 8 * Math.sqrt(2)
                * (edgeLength * edgeLength * edgeLength);
        return volume;
    }

    private static final long serialVersionUID = 1L;
    private double edgeLength;
}
