package resub.substitutionmodel.epochs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import resub.math.ResubMathUtil;


/**
 * Collates information about a single epoch and updates automatically from the state
 */
public class AlphabetEpoch implements Comparable<AlphabetEpoch> {
	
	// The symbol for a state that has been removed from the index list
	public static final String BLANK_STATE = ".";
	
	int index;
	int nStatesFullAlphabet;
	RealParameter transitionAge;
	IntegerParameter alpha;
	IntegerParameter beta;
	IntegerParameter indicator;
	
	int rank;
	List<String> statesAtThisEpoch;
	List<String> statesAtPrevEpoch;
	int nStatesAtEpoch, nStatesAtPrevEpoch;
	
	// Subst model terms
	double[] frequencies;
	private double[] relativeRates;
	double[] transportMatrix;
	
	
	// Inputs to the subst model
	RealParameter substModelFrequencies;
	RealParameter substModelRates;
	GeneralSubstitutionModel substModel;
	
	
	public AlphabetEpoch(int index, RealParameter transitionAge, IntegerParameter alpha, IntegerParameter beta, IntegerParameter indicator, int nStatesFullAlphabet, boolean complexEigen) {
		this.index = index;
		this.transitionAge = transitionAge;
		this.alpha = alpha;
		this.beta = beta;
		this.indicator = indicator;
		this.rank = -1;
		this.frequencies = null;
		this.nStatesFullAlphabet = nStatesFullAlphabet;
		
		initSubstModel(complexEigen);
		
	}
	
	
	private void initSubstModel(boolean complexEigen) {
		
        // Initialise rates for smaller matrix
        List<Double> rates = new ArrayList<>();
        for (int i = 0; i < this.nStatesFullAlphabet*(this.nStatesFullAlphabet-1);  i++) {
        	rates.add(1.0);
        }
        substModelRates = new RealParameter();
        substModelRates.initByName("value", rates);
        
        
        // Initialise freqs
        List<Double> f = new ArrayList<>();
        for (int i = 0; i < this.nStatesFullAlphabet;  i++) {
        	f.add(1.0 / this.nStatesFullAlphabet);
        }
        substModelFrequencies = new RealParameter();
        substModelFrequencies.initByName("value", f);
        Frequencies freqsTopInput = new Frequencies();
        freqsTopInput.initByName("frequencies", substModelFrequencies);
        
        
        // Init reduced subst model
        if (complexEigen) {
        	substModel = new ComplexSubstitutionModel();
        }else {
        	substModel = new GeneralSubstitutionModel();
        }
        
        substModel.initByName("rates", substModelRates, "frequencies", freqsTopInput);
        
        
        // Prepare the transformer matrix
        this.transportMatrix = new double[(this.nStatesFullAlphabet+1)*(this.nStatesFullAlphabet+1)];
        

		
	}


	// Full list of epochs is required to know this
	protected void setStatesAtEpoch(List<String> statesAtEpoch, List<String> statesAtPrevEpoch) {
		this.statesAtThisEpoch = statesAtEpoch;
		this.statesAtPrevEpoch = statesAtPrevEpoch;
		
		this.nStatesAtEpoch = 0;
		for (int i = 0; i < statesAtThisEpoch.size(); i ++) {
			if (!statesAtThisEpoch.get(i).equals(AlphabetEpoch.BLANK_STATE)) nStatesAtEpoch++;
		}
		
		this.nStatesAtPrevEpoch = 0;
		for (int i = 0; i < statesAtPrevEpoch.size(); i ++) {
			if (!statesAtPrevEpoch.get(i).equals(AlphabetEpoch.BLANK_STATE)) nStatesAtPrevEpoch++;
		}
		
		
	}
	
	public int getNStatesAtEpoch() {
		return nStatesAtEpoch;
	}
	
	public int getNStatesAtPrevEpoch() {
		return nStatesAtPrevEpoch;
	}
	
	
	public void setRank(int rank) {
		this.rank = rank;
	}
	
	/**
	 * Epoch rank (ordered by age, with the youngest being rank 0)
	 * @return
	 */
	public int getRank() {
		return this.rank;
	}
	
	
	
	/**
	 * Epoch index (in the state space)
	 * @return
	 */
	public int getIndex() {
		return this.index;
	}
	
	public double getAge() {
		return this.transitionAge.getValue(this.index);
	}
	
	public int getAlpha() {
		return this.alpha.getValue(this.index);
	}
	
	public void setAlpha(int val) {
		this.alpha.setValue(this.index, val);
	}

	public int getBeta() {
		return this.beta.getValue(this.index);
	}
	
	public void setBeta(int val) {
		this.beta.setValue(this.index, val);
	}
	
	public String getAlphaState() {
		int alpha = this.getAlpha();
		return this.statesAtPrevEpoch.get(alpha);
	}
	
	
	public String getAlphaStateSorted() {
		int alpha = this.getAlpha();
		String s = this.statesAtPrevEpoch.get(alpha);
		String[] parts = s.split("/");
		Arrays.sort(parts);
		return String.join("/", parts);
		
	}
	
	public String getBetaState() {
		int beta = this.getBeta();
		return this.statesAtPrevEpoch.get(beta);
	}
	
	public String getBetaStateSorted() {
		int beta = this.getBeta();
		String s = this.statesAtPrevEpoch.get(beta);
		String[] parts = s.split("/");
		Arrays.sort(parts);
		return String.join("/", parts);
	}
	
	public String getTransitionString() {
		
		if (this.getIndicator() == 0) {
			return "null";
		}
		
		
		// Sort the strings so that there is just one representation of a given state
		String a = getAlphaStateSorted();
		String b = getBetaStateSorted();
		if (a.compareTo(b) > 0) {
			String c = a;
			a = b;
			b = c;
		}
		
		if (this.getIndicator() == 1) {
			return a + "/" + b + " -> " + a + "+" + b;
		}
		
		if (this.getIndicator() == 2) {
			return a + " -> " + a + "+" + b;
		}
		
		
		return "";
	}
	
	public String getParentState() {
		
		if (this.getIndicator() == 0) {
			//return getAlphaState() + "+" + getBetaState();
		}
		
		if (this.getIndicator() == 1) {
			return getAlphaState() + "/" + getBetaState();
		}
		
		if (this.getIndicator() == 2) {
			return getAlphaState();
		}
		
		return "null";
	}
	
	public int getIndicator() {
		return this.indicator.getValue(this.index);
	}
	
	
	public double getEffectiveAge() {
		if (!this.isActive()) return Double.POSITIVE_INFINITY;
		return this.getAge();
	}

	@Override
	public int compareTo(AlphabetEpoch o) {
		return Double.compare(this.getEffectiveAge(), o.getEffectiveAge());
	}


	public boolean isActive() {
		return this.getIndicator() != 0;
	}
	
	public boolean isRefinement() {
		return this.getIndicator() == 1;
	}
	
	public boolean isExpansion() {
		return this.getIndicator() == 2;
	}
	
	
	

	
	@Override
	public String toString() {
		String str = "Epoch r=" + this.getRank() + " i=" + this.getIndex() + " ";
		
		if (!this.isActive()) {
			str += "inactive.";
		}
		else {
			if (this.getIndicator() == 1) {
				str += "refinement. ";
			}
			
			if (this.getIndicator() == 2) {
				str += "expansion. ";
			}
			
			str +=  getParentState() + " -> " + getAlphaState() + "+" + getBetaState();
			str += " te=" + this.getAge() + ".";
			
		}
		
		
		str += " N=" + this.getNStatesAtEpoch() + ".";
		
		return str;
	}


	
	
	public String getFrequencyString() {
		String str = "epoch " + this.rank + "\n";
		for (int i = 0; i < nStatesFullAlphabet; i++) {
			if (statesAtThisEpoch.get(i).equals(BLANK_STATE)) continue;
			str += this.statesAtThisEpoch.get(i) + "=" + this.frequencies[i] + "\n";
		}
		return str;
	}
	
	public String getRateString() {
		String str = "epoch " + this.rank + "\n";
		int index = 0;
		for (int i = 0; i < nStatesFullAlphabet; i++) {
			for (int j = 0; j < nStatesFullAlphabet; j++) {
				if (i == j) continue;
				if (j > i && !statesAtThisEpoch.get(i).equals(BLANK_STATE) && !statesAtThisEpoch.get(j).equals(BLANK_STATE)) {
					str += this.statesAtThisEpoch.get(i) + "<=>" + this.statesAtThisEpoch.get(j) + ":" + this.relativeRates[index] + "\n";
				}
				index++;
			}
		}	
		return str;
	}
	
	public List<String> getStateNames(){
		return this.statesAtThisEpoch;
	}
	

	
	public double[] getFrequencies() {
		return this.frequencies;
	}



	public double[] getRelativeRates() {
		return this.relativeRates;
	}

	/**
	 * Set the rates and frequencies, and prepare a substituion model from that
	 * Assuming reversibility
	 * This will be a full 20x20 matrix to assist with multiplication but with some rows and columns blanked out
	 * @param relativeRates
	 * @param frequencies
	 */
	public void prepareSubstitionModel(double[] relativeRates, double[] frequencies, double pStay) {
		
		
//		Frequencies f= substModel.frequenciesInput.get();
//		System.out.print("f before: ");
//		for (int i = 0; i < 20; i ++) {
//			System.out.print(f.getFreqs()[i] + ",");
//		}
//		System.out.println();
		
		
		this.relativeRates = relativeRates;
		this.frequencies = frequencies;
		if (this.frequencies == null || this.relativeRates == null) return;
		
		
		// Ensure the rates have the right dimension
		//this.relativeRates = ResubMathUtil.convertRatesToSymmetric(this.relativeRates, nStatesFullAlphabet);
		
		

		
		// Update freqs in the subst model
		for (int i = 0; i < frequencies.length; i ++) {
			substModelFrequencies.setValue(i, frequencies[i]);
		}
		substModel.frequenciesInput.get().doUpdate();

		
		// Update rates in the subst model
		for (int i = 0; i < relativeRates.length; i ++) {
			substModelRates.setValue(i, relativeRates[i]);
		}
				
		// Update subst model
        this.substModel.doUpdate();
        
        
        
        int alphaIndex = this.getAlpha();
        int betaIndex = this.getBeta();
        
		// Prepare the transformer matrix between this epoch and the next
		// p1 stays in index alpha and p2 changes to beta
		double p1 = pStay;
		double p2 = 1-pStay;
		for (int from = 0; from < this.nStatesFullAlphabet; from++) {
			for (int to = 0; to < this.nStatesFullAlphabet; to++) {
				
				// Diagonal matrix, except for the row being expanded/refined, which has p1 and p2
				int index = ResubMathUtil.getIndex(from, to, nStatesFullAlphabet);
				
				
				if (from != alphaIndex && from == to) {
					this.transportMatrix[index] = 1;
				}
				
				else if (from == alphaIndex && to == alphaIndex) {
					this.transportMatrix[index] = p1;
				}
				
				else if (from == alphaIndex && to == betaIndex) {
					this.transportMatrix[index] = p2;
				}
				
				else {
					this.transportMatrix[index] = 0;
				}
				
			}
		}
		
		
//		System.out.print("f after: ");
//		for (int i = 0; i < 20; i ++) {
//			System.out.print(f.getFreqs()[i] + ",");
//		}
//		System.out.println();
	
		
	}
	
	public GeneralSubstitutionModel getSubstitionModel() {
		if (this.frequencies == null || this.relativeRates == null) return null;
		return this.substModel;
	}
	
	
	public double[] getTransportMatrix() {
		return this.transportMatrix;
	}


	
	
}



