/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

public interface NeighborIteratorHard extends NeighborIterator {

    void iterUpNeighbors(int iAtom, NeighborConsumerHard consumer, double falseTime);

    void iterDownNeighbors(int iAtom, NeighborConsumerHard consumer, double falseTime);

    void iterAllNeighbors(int iAtom, NeighborConsumerHard consumer, double falseTime);
}
