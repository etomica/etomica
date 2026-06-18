/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

/** \brief Pure virtual class from which wall objects are derived.
 *
 * This is a pure virtual class for a generic wall object. A wall object
 * can be specified by deriving a new class from this and specifying the
 * functions.*/
public interface Wall {

     /** A pure virtual function for testing whether a point is
      * inside the wall object. */
     boolean point_inside(double x,double y,double z);
     /** A pure virtual function for cutting a cell without
      * neighbor-tracking with a wall. */
     boolean cut_cell(VoronoiCell c,double x,double y,double z);
     /** A pure virtual function for cutting a cell with
      * neighbor-tracking enabled with a wall. */
     boolean cut_cell(VoronoiCellNeighbor c,double x,double y,double z);

}
