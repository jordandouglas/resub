package resub.math;

import beast.base.core.Log;

public class ResubMathUtil {
	
	/**
	 * Compute result = left * right
	 **/
	public static void multiplyMatrices(double[] left, double[] right, double[] result, int dim) {
		
		
		int index3 = 0;
		for (int rowNum = 0; rowNum < dim; rowNum ++) {
			
			for (int colNum = 0; colNum < dim; colNum ++) {
				
				double sum = 0;
				
				// Move along the columns of left and rows of right 
				int index1 = rowNum * dim;
				int index2 = colNum;
				for (int pos = 0; pos < dim; pos++) {
					double term1 = left[index1];
					double term2 = right[index2]; 
					sum += term1*term2;

					index1 ++;
					index2 += dim;
					
				}

				result[index3] = sum;
				index3 ++;
				
			}
		}
		
		
	}
	
	
	/**
	 * Get the index in a 1D matrix representation of a 2D matrix with dim*dim dimensions
	 * @param from
	 * @param to
	 * @param dim
	 * @return
	 */
	public static int getIndex(int from, int to, int dim) {
		return from*dim + to;
	}
	
	
	
	/**
	 * Copy a double array over (assuming that the latter is already initialised at the right length)
	 * @param copyFrom
	 * @param copyTo
	 */
	public static void copy(double[] copyFrom, double[] copyTo) {
		System.arraycopy(copyFrom, 0, copyTo, 0, Math.min(copyFrom.length, copyTo.length));
	}
	
	
	public static void printProbabilityMatrix(double[] matrix, int dim) {
		
		
		int index = 0;
		for (int i = 0; i < dim; i ++) {
			for (int j = 0; j < dim; j ++) {
				System.out.printf("%.3f\t", matrix[index]);
				index++;
			}
			System.out.println("");
		}
		
		
	}
	

	public static void setRate(double[] rates, int to, int from, double value, int ndim) {
		
		if (from == to) return;
		
		
		int k = 0;
		for (int i = 0; i < ndim; i ++) {
			for (int j = 0; j < ndim; j ++) {
				if (i == j) continue;
				if (i == to && j == from) {
					rates[k] = value;
					return;
				}
				k ++;
			}
		}
	}
	
	
	
	public static double[] convertRatesToSymmetric(double[] ratesIn, int ndim) {
		
//		if (ratesIn.length == ndim * (ndim-1)) {
//			return ratesIn;
//		}
		
		double[] ratesOut = new double[ndim * (ndim - 1)];
        int k = 0;
        for (int i = 0; i < ndim; i++) {
            for (int j = i+1; j < ndim; j++) {
            	double rateVal = ratesIn[k];
            	setRate(ratesOut, i, j, rateVal, ndim);
            	setRate(ratesOut, j, i, rateVal, ndim);
            	k++;
            }
        }
        
        return ratesOut;
		
		
	}
	
	
	/**
	 * For debugging: ensure that each row of the transition probability matrix sums to 1
	 * @param matrix
	 * @param name
	 */
	public static void tidyAndValidateProbs(double[] matrix, String name, int dim, int beta) {
		
		
		// Correct the probs in the row being dropped, just in case of numerical instabilities
		if (beta > -1) {
			int k = 0;
			for (int i = 0; i < dim; i ++) {
				for (int j = 0; j < dim; j ++) {
					if (i == beta || j == beta) {
						if (i == j) {
							matrix[k] = 1;
						}else {
							matrix[k] = 0;
						}
					}
					k++;
				}
				
			}
		}
		
		double maximumError = 0;
		int k = 0;
		for (int i = 0; i < dim; i ++) {
			double rowsum = 0;
			
			int krow = k;
			for (int j = 0; j < dim; j ++) {
				double prob = matrix[k];
				rowsum += prob;
				k++;
			}
			
			maximumError = Math.max(maximumError, Math.abs(1.0-rowsum));
			
			
			if (Math.abs(1.0-rowsum) > 1e-4) {
				//Log.warning(name + " sums to " + rowsum + " on row " + i + " for beta " + beta);
				//throw new IllegalArgumentException();
			}
			
			
			// Normalise it to correct for minor instabilities else they might propagate
			// We really don't want the likelihoods to be too high
			for (int j = 0; j < dim; j ++) {
				matrix[krow] = matrix[krow] / rowsum;
				krow++;
			}
			
		}
		
		// This is way too much error
		if (maximumError > 1e-2) {
			Log.warning("Warning: the probability row sum of " + name + " has an error of " + maximumError + ". If this message persists then the model is unstable.");
		}
		
	}


	
	public static void printProbabilityMatrix(double[][] rateMatrix) {
		for (int i = 0; i < rateMatrix.length; i ++) {
			for (int j = 0; j < rateMatrix.length; j ++) {
				System.out.printf("%.3f\t", rateMatrix[i][j]);
			}
			System.out.println("");
		}
		
	}
	

}
