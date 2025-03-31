package resub.substitutionmodel;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import resub.math.ResubMathUtil;
import resub.substitutionmodel.epochs.AlphabetEpoch;
import resub.substitutionmodel.epochs.AlphabetEpochs;


@Description("Assumes multiple coding epochs, where 2 states coalesce at each transition")
public class MultiTransitionResub extends SubstitutionModel.Base implements Loggable {
	
	final public Input<AlphabetEpochs> epochsInput = new Input<>("epochs", "vector of alphabet epochs", Input.Validate.REQUIRED);

	
	AlphabetEpochs alphabetEpochs;
	GeneralSubstitutionModel substModel;
	
	// Temporary matrices that are initialised only once to save time
	double[] tmpMatrixChaining;
	double[] tmpMatrixEpoch;
	double[] tmpMatrixTransported;
	
	
	// Changes when update() is called
	boolean needsUpdate;
	int nOldActiveEpochs;
	double[] epochTransitionAges;
	
	
	@Override
	public void initAndValidate() {
		
		this.needsUpdate = true;
		this.frequencies = frequenciesInput.get();
        this.nrOfStates = this.frequencies.getFreqs().length;
	        
		this.alphabetEpochs = epochsInput.get();
		this.substModel = alphabetEpochs.getSubstModel();
		if (this.substModel == null) {
			throw new IllegalArgumentException("'epochs' must have a 'substModel'");
		}
		
		this.tmpMatrixChaining = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
		this.tmpMatrixEpoch = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
		this.tmpMatrixTransported = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
		this.epochTransitionAges = new double[this.alphabetEpochs.getNOldEpochs()];
		
		update();
		
	}
	
	@Override
    public double[] getFrequencies() {
		
		
		// Get the frequencies of the oldest epoch
		return alphabetEpochs.getRootFrequencies();
        
    }
	
	
	public void update() {
		
		if (!needsUpdate) return;
		
		// Get the age of each epoch (sorted in increasing order of age; where the age is infinity if the epoch is inactive)
		this.nOldActiveEpochs = alphabetEpochs.getNActiveEpochs();
		for (int i = 0; i < alphabetEpochs.getNOldEpochs(); i ++) {
			AlphabetEpoch epoch = alphabetEpochs.getEpochByRank(i);
			this.epochTransitionAges[i] = epoch.getEffectiveAge(); 
		}
		
		needsUpdate = false;
		
	}
	
	
	@Override
	protected boolean requiresRecalculation() {
		if (InputUtil.isDirty(epochsInput)) needsUpdate = true;
		if (InputUtil.isDirty(frequenciesInput)) needsUpdate = true;
		return needsUpdate;
	}
	
	
	@Override
    public void store() {
		super.store();
	}
	
	
	
	@Override
	public void restore() {
		this.needsUpdate = true;
		super.restore();
	}

	

	@Override
	public void getTransitionProbabilities(Node node, double startAge, double endAge, double rate, double[] outMatrix) {
		
		

		double distance = (startAge - endAge) * rate;
		if (distance < 0) {
			throw new IllegalArgumentException("Error 1: negative branch length detected! " + distance);
		}
		
		
		update();
		double bottomEpochStartAge = this.epochTransitionAges[0];
		
		
		// Use the full-alphabet vanilla subst model
		if (nOldActiveEpochs == 0 || startAge <= bottomEpochStartAge) {
			this.substModel.getTransitionProbabilities(node, startAge, endAge, rate, outMatrix);
			ResubMathUtil.tidyAndValidateProbs(outMatrix, "full matrix 1", this.nrOfStates, -1);
			return;
		}
		
		
		//System.out.println(this.getClass().getSimpleName() + " " + bottomEpochStartAge);
		
		
		// For each epoch, ranked in decreasing order of age. rank=-1 is the present day full alphabet
		// We will start from the oldest epoch and work our way down
		// This means that the probability matrices are chained together in the forward time direction
		boolean outMatrixInitialised = false;
		for (int epochRank = nOldActiveEpochs-1; epochRank >= -1; epochRank--) {
			
			//System.out.println(this.getClass().getSimpleName() + " " + epochRank);
			
			
			int indexBeingDropped = (epochRank == -1) ? -1 : this.alphabetEpochs.getEpochByRank(epochRank).getBeta();
			double epochStartAge = (epochRank == nOldActiveEpochs-1) ? Double.POSITIVE_INFINITY : this.epochTransitionAges[epochRank+1];
			double epochEndAge = (epochRank == -1) ? 0 : this.epochTransitionAges[epochRank];
			
			
			// This should not have happened
			if (epochEndAge == Double.POSITIVE_INFINITY) {
				throw new IllegalArgumentException("The active epoch age is inf " + epochRank);
			}
			
			// Does this branch cross this epoch?
			boolean lineageCrossesEpoch = (startAge >= epochEndAge && endAge < epochStartAge);
			if (!lineageCrossesEpoch) continue;
			
			
			// Find the interval of overlap
			double h0 = Math.min(startAge, epochStartAge);
			double h1 = Math.max(endAge, epochEndAge);
			
			
			// Get the substitution model and transformer matrix for this epoch
			GeneralSubstitutionModel substModelEpoch = epochRank == -1 ? this.substModel : alphabetEpochs.getEpochSubstitionModelByRank(epochRank);
			
			//System.out.println("reduced rate matrix:");
			//substModelEpoch.doUpdate();
			
			// Get the probability matrix for this epoch
			substModelEpoch.getTransitionProbabilities(node, h0, h1, rate, tmpMatrixEpoch);
			ResubMathUtil.tidyAndValidateProbs(tmpMatrixEpoch, "reduced matrix at rank " + epochRank, this.nrOfStates, indexBeingDropped);
			
			
			//ResubMathUtil.printProbabilityMatrix(substModelEpoch.getRateMatrix());
			
			
			//System.out.println("reduced matrix:");
			//ResubMathUtil.printProbabilityMatrix(tmpMatrixEpoch, this.nrOfStates);
			
			// Does this branch have a second epoch?
			boolean thereIsAnotherEpoch = epochRank >= 0 && endAge <= epochEndAge + 1e-16;
			if (thereIsAnotherEpoch) {
				
				// Multiply by transporter matrix
				double[] transportMatrix = alphabetEpochs.getEpochTransportMatrixByRank(epochRank);
				
//				ResubMathUtil.tidyAndValidateProbs(transportMatrix, "transport matrix at rank " + epochRank, this.nrOfStates, -1);
//				System.out.println("transportMatrix:");
//				ResubMathUtil.printProbabilityMatrix(transportMatrix, this.nrOfStates);
//				
				
				ResubMathUtil.multiplyMatrices(tmpMatrixEpoch, transportMatrix, tmpMatrixTransported, this.nrOfStates);
				ResubMathUtil.tidyAndValidateProbs(tmpMatrixTransported, "transported matrix at rank " + epochRank, this.nrOfStates, -1);
				ResubMathUtil.copy(tmpMatrixTransported, tmpMatrixEpoch);
				
//				System.out.println("tmpMatrixEpoch:");
//				ResubMathUtil.printProbabilityMatrix(tmpMatrixEpoch, this.nrOfStates);
//			
			}
			
			
			if (outMatrixInitialised) {
				
				// Multiply A*B where A is the pmatrix from the previous (older) epoch and B is this epoch
				ResubMathUtil.multiplyMatrices(outMatrix, tmpMatrixEpoch, tmpMatrixChaining, this.nrOfStates);
				ResubMathUtil.tidyAndValidateProbs(tmpMatrixChaining, "chained matrix at rank " + epochRank, this.nrOfStates, -1);
				ResubMathUtil.copy(tmpMatrixChaining, outMatrix);
				
			}
			
			else {
				
				ResubMathUtil.copy(tmpMatrixEpoch, outMatrix);
				outMatrixInitialised = true;
			}
			
			
			// All done
			if (!thereIsAnotherEpoch) break;
			
			
		}
		
		
				
		
	}

	


	@Override
	public void init(PrintStream out) {
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		
	}

	@Override
	public void close(PrintStream out) {
		
	}
	
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType.getClass().equals(epochsInput.get().getDataType().getClass());
	}
	
	
	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		return null;
	}
	
	@Override
    public boolean canReturnComplexDiagonalization() {
		
		// Set to true so that beagle will let us do the matrix exponentiation work
        return true;
    }

}


