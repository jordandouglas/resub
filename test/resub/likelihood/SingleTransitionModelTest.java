package resub.likelihood;


import org.junit.Assert;
import org.junit.Test;

import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import resub.substitutionmodel.SVSGeneralSubstitutionModel;
import resub.substitutionmodel.SingleTransitionResub;


public class SingleTransitionModelTest {
	
	public static final double PRECISION = 1e-6;
	public static final String RATES_STR = "4.432 0.3791 0.5617 0.3362 2.015 0.4722 0.2673 0.4467 0.321 1.228 0.3354 0.8666 0.6356 0.3921 3.992 1.748 2.319 0.1529 0.5584 0.2803 0.1374 0.9365 0.72 0.5864 0.5093 0.1167 1.795 1.164 2.392 0.5706 0.3118 0.5543 3.314 2.908 4.416 0.5549 0.7716 3.573 0.01 0.9316 0.5865 0.01 0.8569 0.01 0.01 4.696 0.439 0.9639 0.2175 1.326 0.569 0.1409 0.2392 0.1486 0.1159 0.1641 0.4819 0.01 1.411 0.1136 0.08969 1.256 0.327 2.988 0.4909 0.6071 0.5711 0.2406 0.01 0.1393 0.01 0.9002 1.036 0.01 1.698 1.859 0.01 0.07376 0.01 0.07911 0.1126 0.1634 0.7343 2.039 5.874 0.2773 0.01 0.4422 0.03281 0.1001 1.566 0.2708 0.1631 0.07232 1.353 0.3508 0.01 0.09905 0.06602 0.2757 0.7037 0.2954 0.4792 4.619 0.2693 2.937 0.8973 1.689 1.271 0.3938 0.844 3.762 0.1303 3.875 3.457 0.1901 0.1478 0.3567 0.1924 0.01 0.7648 8.736 0.09717 0.2367 0.2786 0.4948 2.082 0.7338 3.489 4.43 1.33 0.8052 0.2321 0.01 0.4391 5.448 0.153 0.1993 0.5316 0.2138 0.115 0.5175 1.843 0.3811 0.479 0.2044 0.2669 2.433 0.3585 0.6762 1.537 1.743 0.7427 0.7593 0.6591 2.387 1.513 3.206 1.63 0.3433 0.09767 1.115 0.3913 0.2988 1.019 0.7811 0.4544 0.09392 0.26 1.398 1.005 1.187 0.3361 0.4068 0.3809 0.8003 0.5168 0.1538 0.3695 0.2763 5.064 0.4246 0.3724 0.2348 2.892 0.4379 0.5008 0.2431 0.5531 2.117";
	public static final String FREQ_STR = "0.03 0.07 0.04 0.06 0.01 0.09 0.02 0.08 0.05 0.05 0.033 0.067 0.1 0.03 0.03 0.04 0.05 0.05 0.05 0.05";
	
	
	
	@Test
    public void testSVSComparison() throws Exception {
		
		System.out.println("\n--------------------------------------------");
		System.out.println("Performing testSVSComparison");
		
		SVSGeneralSubstitutionModel svsModel = new SVSGeneralSubstitutionModel();
		SingleTransitionResub cherryResub = new SingleTransitionResub();
		getMatrices(svsModel, cherryResub, true, true, true);
		
		
		// Get probabilities
		double startHeight = 0.5;
		double endHeight = 0.1;
		double rate = 1.0;
		double[] outMatrixResub = new double[20*20];
		double[] outMatrixSVS = new double[20*20];
		cherryResub.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixResub);
		svsModel.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixSVS);
		
		
		// Confirm that all rows sum to 1
		int k = 0;
		for (int i = 0; i < 20; i ++) {
			double rowSum = 0;
			for (int j = 0; j < 20; j ++) {
				rowSum += outMatrixResub[k];
				k++;
			}
			
			Assert.assertEquals(rowSum, 1.0, PRECISION);
			System.out.println("row " + i + " sums to " + rowSum);
			
		}
		
		
		
		// Confirm that resub gives same answer as svs when the branch is under the transition
		System.out.println("Confirming that SVS and resub give the same answer on branches younger than the transition...");
		Assert.assertArrayEquals(outMatrixSVS, outMatrixResub, PRECISION);
		
		
		
		// Confirm that resub and svs give different answers when the branch crosses the transition
		startHeight = 1.1;
		endHeight = 0.9;
		cherryResub.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixResub);
		svsModel.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixSVS);
		
		
		// Confirm that all rows sum to 1
		k = 0;
		for (int i = 0; i < 20; i ++) {
			double rowSum = 0;
			for (int j = 0; j < 20; j ++) {
				rowSum += outMatrixResub[k];
				k++;
			}
			
			Assert.assertEquals(rowSum, 1.0, PRECISION);
			System.out.println("row " + i + " sums to " + rowSum);
			
		}
		
		
		System.out.println("Confirming that SVS and resub give the different answers on branches that cross the transition...");
		int ndiff = 0;
		for (int i = 0; i < 20; i ++) {
			if (Math.abs(outMatrixSVS[i] - outMatrixResub[i]) > PRECISION) {
				ndiff++;
			}
		}
		Assert.assertTrue(ndiff > 0);
		
		System.out.println("Test passed!");
		System.out.println("--------------------------------------------\n");
		
	}
	
	
	
	@Test
    public void test19x19() throws Exception {
		
		
		final int lastState = 19;
		
		// Does the 20x20 matrix with a missing row/column give the same exponentiation as a 19x19?
		System.out.println("\n--------------------------------------------");
		System.out.println("Performing test19x19");
		
		
		// Build the 20x20 matrix
		SVSGeneralSubstitutionModel svsModel = new SVSGeneralSubstitutionModel();
		SingleTransitionResub cherryResub = new SingleTransitionResub();
		getMatrices(svsModel, cherryResub, true, true, false);
		
		
		// Build an equivalent 19x19
		SVSGeneralSubstitutionModel svs19 = new SVSGeneralSubstitutionModel();
		get19x19Matrix(svs19);
		
		
		// Do the two models have the same frequencies?
		double[] freqs1 = cherryResub.getFrequencies();
		double[] freqs2 = svs19.getFrequencies();
		for (int i = 0; i < 19; i ++) {
			double f1 = freqs1[i];
			double f2 = freqs2[i];
			Assert.assertEquals(f1, f2, PRECISION);
			System.out.println("pi" + i + " = " + f1 + " = " + f2);
		}
		
		
		// Get probabilities of a branch above the transition
		double startHeight = 2.0;
		double endHeight = 1.75;
		double rate = 1.0;
		double[] outMatrixResub = new double[20*20];
		double[] outMatrix19x19 = new double[19*19];
		cherryResub.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixResub);
		svs19.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrix19x19);
		
		
		// Confirm that the 19x19 rows sum to 1
		int k = 0;
		for (int i = 0; i < 19; i ++) {
			double rowSum = 0;
			for (int j = 0; j < 19; j ++) {
				rowSum += outMatrix19x19[k];
				k++;
			}
			Assert.assertEquals(rowSum, 1.0, PRECISION);
			System.out.println("row " + i + " sums to " + rowSum);
			
		}
		
		
		// Confirm that two matrices are the same even through they are different dimensions
		k = 0;
		int l = 0;
		for (int i = 0; i < 20; i ++) {
			for (int j = 0; j < 20; j ++) {
				
				if (i == lastState) {
					
					// Do nothing
					
				}else if (j == lastState) {
					double p = outMatrixResub[k];
					Assert.assertEquals(p, 0.0, PRECISION);
					//System.out.println(i + " dead state has p entry " + p);
					System.out.println("dead state p(" + i + "," + j + ") = " + 0);
				}else {
					double p1 = outMatrixResub[k];
					double p2 = outMatrix19x19[l];
					
					Assert.assertEquals(p1, p2, PRECISION);
					System.out.println("p(" + i + "," + j + ") = " + p1 + " = " + p2);
					
					l++;
				}
				
				k++;
			}
			
			
		}
				
		
		// Confirm that the two converge to the same equilibrium distribution
		startHeight = 10000.0;
		endHeight = 2.0;
		cherryResub.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixResub);
		svs19.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrix19x19);
		
		k = 0;
		for (int i = 0; i < 20; i ++) {
			
			
			for (int j = 0; j < 20; j ++) {
				
				if (i != lastState) {
					double expectedP = freqs1[j];
					double p = outMatrixResub[k];
					Assert.assertEquals(expectedP, p, PRECISION);
					System.out.println(i + " " + j + " equilibrium frequency is " + expectedP + " = " + p);
				}
				k++;
			}
			
		}
		
		
		// Confirm that the equilibrium frequencies are not reached at after time 0
		startHeight = 2.0;
		endHeight = startHeight - 0.01;
		cherryResub.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrixResub);
		svs19.getTransitionProbabilities(null, startHeight, endHeight, rate, outMatrix19x19);
		
		k = 0;
		for (int i = 0; i < 20; i ++) {
			
			
			for (int j = 0; j < 20; j ++) {
				
				if (j != lastState && i != lastState) {
					double expectedP = freqs1[j];
					double p = outMatrixResub[k];
					Assert.assertTrue(Math.abs(expectedP-p) > PRECISION);
					System.out.println(i + " " + j + " frequency is NOT the equilibrium " + p + " != " + expectedP);
				}
				k++;
			}
			
		}
				
		
		
		// Confirm that the 19x19 rows sum to 1
		k = 0;
		for (int i = 0; i < 19; i ++) {
			double rowSum = 0;
			for (int j = 0; j < 19; j ++) {
				rowSum += outMatrix19x19[k];
				k++;
			}
			Assert.assertEquals(rowSum, 1.0, PRECISION);
			System.out.println("row " + i + " sums to " + rowSum);
			
		}
		
		
		// Confirm that two matrices are the same even through they are different dimensions
		k = 0;
		l = 0;
		for (int i = 0; i < 20; i ++) {
			for (int j = 0; j < 20; j ++) {
				
				if (i == lastState) {
					
					// Do nothing
					
				}else if (j == lastState) {
					double p = outMatrixResub[k];
					Assert.assertEquals(p, 0.0, PRECISION);
					//System.out.println(i + " dead state has p entry " + p);
					System.out.println("dead state p(" + i + "," + j + ") = " + 0);
				}else {
					double p1 = outMatrixResub[k];
					double p2 = outMatrix19x19[l];
					
					Assert.assertEquals(p1, p2, PRECISION);
					System.out.println("p(" + i + "," + j + ") = " + p1 + " = " + p2);
					
					l++;
				}
				
				k++;
			}
			
			
		}
		
		

		
		System.out.println("Test passed!");
		System.out.println("--------------------------------------------\n");
		
	}
	
	
	// Build a 19x19 matrix with the final state (Y) dropped
	private void get19x19Matrix(SVSGeneralSubstitutionModel svsModel) {
		
		
		// Build rates from 19x19
		String rates_str_19 = "";
		String[] bits = RATES_STR.split(" ");
		int k = 0;
		for (int i = 0; i < 20; i ++) {
			for (int j = i; j < 20; j ++) {
				
				//System.out.println(i + " " + j + " " + k);
				if (i == j) continue;
				
				if (i != 19 && j != 19) {
					rates_str_19 += bits[k] + " ";
					//System.out.println(bits[k]);
				}
				
				k++;
			}
		}
		
		RealParameter rates = new RealParameter(rates_str_19);
		//rates.initByName("dimension", 171);
		
		
		// Frequencies with just 19 elements
		String[] bits2 = FREQ_STR.split(" ");
		double[] freqs19 = new double[19];
		String freqs_str_19 = "";
		double fsum = 0;
		for (int i = 0; i < 19; i ++) {
			freqs19[i] = Double.parseDouble(bits2[i]);
			fsum += freqs19[i];
		}
		for (int i = 0; i < 19; i ++) {
			freqs_str_19 += "" + (freqs19[i] / fsum) + " ";
		}
		
		RealParameter freqs = new RealParameter(freqs_str_19);
		//freqs.initAndValidate();
		Frequencies frequencies = new Frequencies();
		frequencies.initByName("frequencies", freqs);
		
		
		// Rate indicator
		BooleanParameter indicator = new BooleanParameter("true");
		indicator.initByName("dimension", 171);
		
		
		// SVS
		svsModel.initByName(
								"symmetric", true,
								"rates", rates,
								"frequencies", frequencies,
								"rateIndicator", indicator
							);
				
		
	}
		
	
	

	private void getMatrices(SVSGeneralSubstitutionModel svsModel, SingleTransitionResub cherryResub, boolean useResubVal, boolean expandVal, boolean refineVal) {
		

		// Rates
		RealParameter rates = new RealParameter(RATES_STR);
		rates.initByName("dimension", 190);
		
		
		// Frequencies
		RealParameter freqs = new RealParameter(FREQ_STR);
		freqs.initAndValidate();
		Frequencies frequencies = new Frequencies();
		frequencies.initByName("frequencies", freqs);
		
		
		// Rate indicator
		BooleanParameter indicator = new BooleanParameter("true");
		indicator.initByName("dimension", 190);
		
		
		// SVS
		svsModel.initByName(
								"symmetric", true,
								"rates", rates,
								"frequencies", frequencies,
								"rateIndicator", indicator
							);
		
		
		
		// Resub parameters
		RealParameter h1 = new RealParameter("1.0");
		RealParameter pi = new RealParameter("0.5");
		BooleanParameter useResub = new BooleanParameter("" + useResubVal);
		BooleanParameter refine = new BooleanParameter("" + refineVal);
		BooleanParameter expand = new BooleanParameter("" + expandVal);
		
		
		// Make the resub model
		
		cherryResub.initByName(
								"state1", "W",
								"state2", "Y",
								"h1", h1,
								"pi", pi,
								"substModel", svsModel,
								"frequencies", frequencies,
								"useResub", useResub,
								"refine", refine,
								"expand", expand
							);
		
		
		
	}
	
	

}
