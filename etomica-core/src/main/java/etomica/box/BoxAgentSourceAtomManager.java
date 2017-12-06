/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.BoxAgentManager.BoxAgentSource;

/**
 * BoxAgentSource that creates an AtomLeafAgentManager.
 *
 * @param <E> atom agent class
 */
public class BoxAgentSourceAtomManager<E> implements BoxAgentSource<AtomLeafAgentManager<E>> {

    protected final AgentSource<E> atomAgentSource;
    protected final Class atomAgentClass;

    /**
     * @param atomAgentSource object that makes the atom agents
     * @param atomAgentClass Class specified by &lt;E&gt;
     */
    public BoxAgentSourceAtomManager(AgentSource<E> atomAgentSource, Class atomAgentClass) {
        super();
        this.atomAgentSource = atomAgentSource;
        this.atomAgentClass = atomAgentClass;
    }

    public AtomLeafAgentManager<E> makeAgent(Box box) {
        return new AtomLeafAgentManager<E>(atomAgentSource, box, atomAgentClass);
    }

    public void releaseAgent(AtomLeafAgentManager<E> agent) {
        agent.dispose();
    }
}
