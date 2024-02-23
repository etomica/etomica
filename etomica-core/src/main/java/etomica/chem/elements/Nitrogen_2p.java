package etomica.chem.elements;

    public class Nitrogen_2p extends ElementChemical {

        protected Nitrogen_2p(String symbol) {
            this(symbol, 14.01);
        }

        protected Nitrogen_2p(String symbol, double mass) {
            super(symbol, mass, 7);
        }

        public static final Nitrogen_2p INSTANCE = new Nitrogen_2p("N_2p");

        private static final long serialVersionUID = 1L;
}
