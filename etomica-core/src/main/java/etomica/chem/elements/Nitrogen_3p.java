package etomica.chem.elements;

    public class Nitrogen_3p extends ElementChemical {

        protected Nitrogen_3p(String symbol) {
            this(symbol, 14.01);
        }

        protected Nitrogen_3p(String symbol, double mass) {
            super(symbol, mass, 7);
        }

        public static final Nitrogen_3p INSTANCE = new Nitrogen_3p("N_3p");

        private static final long serialVersionUID = 1L;
}
