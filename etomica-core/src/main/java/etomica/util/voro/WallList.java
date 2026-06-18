/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

import java.util.Arrays;

/** \brief A class for storing a list of pointers to walls.
 *
 * This class stores a list of pointers to wall classes. It contains several
 * simple routines that make use of the wall classes (such as telling whether a
 * given position is inside all of the walls or not). It can be used by itself,
 * but also forms part of container_base, for associating walls with this
 * class. */
public class WallList {

    /** An array holding pointers to wall objects. */
    Wall[] walls;
    /** A pointer to the next free position to add a wall pointer.
     */
    int wep;

    /** The wall_list constructor sets up an array of pointers to wall classes. */
    public WallList() {
        walls = new Wall[Config.init_wall_size];
        wep = 0;
    }

    /** Adds a wall to the list.
     * \param[in] w the wall to add. */
    public void add_wall(Wall w) {
        if(wep==walls.length) increase_wall_memory();
        walls[wep] = w;
        wep++;
    }
    void add_wall(WallList wl) {
        for (int i=0; i<wl.wep; i++) add_wall(wl.walls[i]);
    }
    /** Determines whether a given position is inside all of the
     * walls on the list.
     * \param[in] (x,y,z) the position to test.
     * \return True if it is inside, false if it is outside. */
    public boolean point_inside_walls(double x,double y,double z) {
        for (int i=0; i<wep; i++) if (!walls[i].point_inside(x,y,z)) return false;
        return true;
    }
    /** Cuts a Voronoi cell by all of the walls currently on
     * the list.
     * \param[in] c a reference to the Voronoi cell class.
     * \param[in] (x,y,z) the position of the cell.
     * \return True if the cell still exists, false if the cell is
     * deleted. */
    public boolean apply_walls(VoronoiCellBase c, double x, double y, double z) {
        if (c instanceof VoronoiCell) for (int i=0; i<wep; i++) if(!(walls[i].cut_cell((VoronoiCell) c,x,y,z))) return false;
        if (c instanceof VoronoiCellNeighbor) for (int i=0; i<wep; i++) if(!(walls[i].cut_cell((VoronoiCellNeighbor) c,x,y,z))) return false;
        return true;
    }
    public void deallocate() {
        Arrays.fill(walls, null);
        wep = 0;
    }

    /** Increases the memory allocation for the walls array. */
    protected void increase_wall_memory() {
        int new_size = walls.length<<1;
        if (new_size>Config.max_wall_size) {
            Common.voro_fatal_error("Wall memory allocation exceeded absolute maximum",Config.Voropp.MEMORY_ERROR);
        }
        walls = Arrays.copyOf(walls, new_size);
    }

}
