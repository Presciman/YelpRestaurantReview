import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Random;
//import org.apache.commons.math.special.Gamma;


public class FLDA extends TopicModel {
	public HashMap<String,Integer> wordMap;
	public HashMap<Integer,String> wordMapInv;

	public HashMap<Integer,int[]> iToVectorCache;
	public HashMap<int[],Integer> vectorToICache;

	public String priorPrefix;

	public int[][] docs;
	public int[] docsC;
	public int[][] docsZ;
	public int[][][] docsZZ;

	public int[][] nDZ;
	public int[] nD;
	public int[][] nZW;
	public int[] nZ;

	public int K;
	public int D;
	public int W;
	public int[] Z;
	public int Z0; // size of Cartesian product
	public int[] Zsub;

	public double alphaB;
	public double[][] alphaZ;	
	public double[][][] alphaDZ;
	public double[][] priorDZ;
	public double[] alphaNorm;

	public double omegaB;
	public double[] omegaW;
	public double[][][] omegaZW;
	public double[][] priorZW;
	public double[] omegaNorm;
	public double[] etaW;
	public double[][][] etaZW;

	public double[] beta;
	
	public double stepSizeADZ;
	public double stepSizeAZ;
	public double stepSizeAB;
	public double stepSizeW;
	public double stepSizeWB;
	public double stepSizeB;

	public double delta0;
	public double delta1;

	public double sigmaA;
	public double sigmaAB;
	public double sigmaW;
	public double sigmaWB;

	public int likelihoodFreq;
	public int blockFreq;

	public FLDA(int[] z, double sigmaA0, double sigmaAB0, double sigmaW0, double sigmaWB0, double stepSizeADZ0, double stepSizeAZ0, double stepSizeAB0, double stepSizeW0, double stepSizeWB0, double stepSizeB0, double delta00, double delta10, double alphaB0, double omegaB0, int likelihoodFreq0, int blockFreq0, String prefix) {
		K = z.length;
		Z = new int[K];
		Zsub = new int[K];

		Z0 = 1;
		for (int i = K-1; i >= 0; i--) {
			Z[i] = z[i];
			Zsub[i] = Z0;
			Z0 *= Z[i];
		}

		iToVectorCache = new HashMap<Integer, int[]>();
		for (int x = 0; x < Z0; x++) {
			int[] z0 = iToVector0(x);
			iToVectorCache.put(x, z0);
		}

		sigmaA = sigmaA0;
		sigmaAB = sigmaAB0;
		sigmaW = sigmaW0;
		sigmaWB = sigmaWB0;
		stepSizeADZ = stepSizeADZ0;
		stepSizeAZ = stepSizeAZ0;
		stepSizeAB = stepSizeAB0;
		stepSizeW = stepSizeW0;
		stepSizeWB = stepSizeWB0;
		stepSizeB = stepSizeB0;
		delta0 = delta00;
		delta1 = delta10;
		alphaB = alphaB0;
		omegaB = omegaB0;

		likelihoodFreq = likelihoodFreq0;
		blockFreq = blockFreq0;
		priorPrefix = prefix;
	}
	
	public void initialize() {
		System.out.println("Initializing...");
		Random r = new Random();

		System.out.println("sigmaA = "+sigmaA);
		System.out.println("sigmaW = "+sigmaW);
		System.out.println("stepSizeADZ = "+stepSizeADZ);
		System.out.println("stepSizeAZ = "+stepSizeAZ);
		System.out.println("stepSizeAB = "+stepSizeAB);
		System.out.println("stepSizeW = "+stepSizeW);
		System.out.println("stepSizeWB = "+stepSizeWB);
		System.out.println("stepSizeB = "+stepSizeB);
		System.out.println("delta0 = "+delta0);
		System.out.println("delta1 = "+delta1);
		System.out.println("alphaB = "+alphaB);
		System.out.println("omegaB = "+omegaB);

		alphaZ = new double[K][];
		alphaDZ = new double[K][][];
		for (int i = 0; i < K; i++) {
			alphaZ[i] = new double[Z[i]];
			alphaDZ[i] = new double[Z[i]][D];
		}
		priorDZ = new double[D][Z0];
		alphaNorm = new double[D];

		omegaW = new double[W];
		omegaZW = new double[K][][];
		priorZW = new double[Z0][W];
		omegaNorm = new double[Z0];
		for (int i = 0; i < K; i++) omegaZW[i] = new double[Z[i]][W];

		beta = new double[Z0];

		docsZ = new int[D][];
		docsZZ = new int[D][][];

		nDZ = new int[D][Z0];
		nD = new int[D];
		nZW = new int[Z0][W];
		nZ = new int[Z0];
	
		for (int w = 0; w < W; w++) {
			omegaW[w] = etaW[w];
		}
		for (int i = 0; i < K; i++) {
			for (int z = 0; z < Z[i]; z++) {
				for (int w = 0; w < W; w++) {
					omegaZW[i][z][w] = etaZW[i][z][w];
				}
			}
		}
		for (int i = 0; i < Z0; i++) {
			omegaNorm[i] = 0;
			for (int w = 0; w < W; w++) {
				priorZW[i][w] = priorW(w, i);
				omegaNorm[i] += priorZW[i][w];
			}
		}
		
		for (int d = 0; d < D; d++) {
			alphaNorm[d] = 0;

			for (int i = 0; i < Z0; i++) {
				priorDZ[d][i] = priorA(d, i); 
				alphaNorm[d] += priorDZ[d][i];
			}
		}

		for (int d = 0; d < D; d++) { 
			docsZ[d] = new int[docs[d].length];
			docsZZ[d] = new int[docs[d].length][Z0];
			
			for (int n = 0; n < docs[d].length; n++) {
				int w = docs[d][n];

				//int z = r.nextInt(Z0); // sample uniformly
				//docsZ[d][n] = z;

				// sample initial value from word priors
				int z = -1;
				double[] p = new double[Z0];
				double pTotal = 0;
				double pMax = 0;
				for (int i = 0; i < Z0; i++) {	
					p[i] = priorZW[i][w];
					pTotal += p[i];
					if (p[i] > pMax) {
						pMax = p[i]; // record max in case you want to initialize to max
					}
				}
				double u = r.nextDouble() * pTotal;
				double v = 0.0;
				for (int i = 0; i < Z0; i++) {
					v += p[i];
					if (v > u) {
						z = i;
						break;
					}
				}
				docsZ[d][n] = z;
				
				// update counts

				nZW[z][w] += 1;	
				nZ[z] += 1;				
				nDZ[d][z] += 1;
				nD[d] += 1;
			}
		}
	}

	public int[] iToVector0(int x) {
		int[] z = new int[K];

		for (int i = 0; i < K; i++) {
			z[i] = x / Zsub[i];
			x = x % Zsub[i];
		}

		return z;
	}

	public int[] iToVector(int x) {
		return iToVectorCache.get(x);
	}

	public int vectorToI(int[] z) {
		int x = 0;

		for (int i = 0; i < K; i++) {
			x += z[i] * Zsub[i];
		}

		return x;
	}

	// returns the alpha_dz prior given all the parameters
	public double priorA(int d, int x) {
		double weight = alphaB;
		int[] z = iToVector(x);

		for (int i = 0; i < K; i++) {
			weight += alphaZ[i][z[i]] + alphaDZ[i][z[i]][d];
		}

		double b = logistic(beta[x]);
		return b * Math.exp(weight);
	}

	// returns the omega_zw prior given all the parameters
	public double priorW(int w, int x) {
		double weight = omegaB + omegaW[w];
		int[] z = iToVector(x);

		for (int i = 0; i < K; i++) {
			weight += omegaZW[i][z[i]][w];
		}

		return Math.exp(weight);
	}

	public double logistic(double x) {
		return 1.0 / (1.0 + Math.exp(-1.0*x));
	}

	// derivative of logistic function
	public double dlogistic(double x) {
		return logistic(x) * (1.0 - logistic(x));
	}

	public void updateWeights(int iter) {
		updateWeightsW(iter);
		updateWeightsA(iter);
	}

	// gradient ascent update on all alpha params
	public void updateWeightsA(int iter) {
		double sigma = sigmaA;

		double[] gradientBeta = new double[Z0];
		double gradientB = 0;
		double[][] gradientZ = new double[K][];
		for (int i = 0; i < K; i++) {
			gradientZ[i] = new double[Z[i]];
		}

		for (int d = 0; d < D; d++) {
			double[][] gradientDZ = new double[K][];
			for (int i = 0; i < K; i++) {
				gradientDZ[i] = new double[Z[i]];
			}

			for (int x = 0; x < Z0; x++)
			{
				int[] z = iToVector(x);

				double dg1 = digamma0(alphaNorm[d]);
				double dg2 = digamma0(alphaNorm[d] + docs[d].length);
				double dgW1 = digamma0(priorDZ[d][x] + nDZ[d][x]);
				double dgW2 = digamma0(priorDZ[d][x]);
				
				double gradientLL = priorDZ[d][x]*(dg1-dg2+dgW1-dgW2);

				for (int i = 0; i < K; i++) {
					gradientZ[i][z[i]] += gradientLL;
					gradientDZ[i][z[i]] += gradientLL;
				}
				gradientB += gradientLL;

				// exp(...)*P(beta)*(1.0-P(beta))
				// priorDZY = exp(...)*P(beta)
				// where P is logistic
				gradientBeta[x] += gradientLL * (1.0-logistic(beta[x]));
			}

			for (int i = 0; i < K; i++) {
				for (int z = 0; z < Z[i]; z++) {
					gradientDZ[i][z] += -alphaDZ[i][z][d] / (sigma*sigma); 
	
					alphaDZ[i][z][d] = alphaDZ[i][z][d] + stepSizeADZ*gradientDZ[i][z];
					//if (d % 1000 == 0) System.out.println(d+" alphaDZ["+i+"]["+z+"]["+d+"] = "+alphaDZ[i][z][d]);
				}
			}
		}	

		for (int i = 0; i < K; i++) {
			for (int z = 0; z < Z[i]; z++) {
				gradientZ[i][z] += -alphaZ[i][z] / (sigma*sigma); 
				alphaZ[i][z] = alphaZ[i][z] + (stepSizeAZ)*gradientZ[i][z];
				System.out.println("alphaZ_"+i+","+z+" "+alphaZ[i][z]);
			}
		}

		gradientB += -alphaB / (sigmaAB*sigmaAB);
		alphaB = alphaB + (stepSizeAB)*gradientB;
		System.out.println("alphaB "+alphaB);

		if (iter < 100) return;  // don't update beta until after small burn-in
		for (int x = 0; x < Z0; x++) {

			gradientBeta[x] += (delta0-1.0) * dlogistic(beta[x]) / logistic(beta[x]);
			gradientBeta[x] += (delta1-1.0) * (-1.0*dlogistic(beta[x])) / (1.0-logistic(beta[x]));

			beta[x] = beta[x] + (stepSizeB)*gradientBeta[x];
			System.out.println("beta_"+x+" "+beta[x]+" (b="+logistic(beta[x])+")");
		}
	}

	// gradient ascent update on all omega params
	public void updateWeightsW(int iter) {
		if (iter < 100) return; // don't update omega until after small burn-in

		double sigma = sigmaW;

		double gradientB = 0;

		for (int w = 0; w < W; w++) {
			double gradientW = 0;
			double[][] gradientZW = new double[K][];
			for (int i = 0; i < K; i++) gradientZW[i] = new double[Z[i]];

			for (int x = 0; x < Z0; x++)
			{
				int[] z = iToVector(x);

				double dg1 = digamma0(omegaNorm[x]);
				double dg2 = digamma0(omegaNorm[x] + nZ[x]);
				double dgW1 = digamma0(priorZW[x][w] + nZW[x][w]);
				double dgW2 = digamma0(priorZW[x][w]);
				
				double gradientLL = priorZW[x][w]*(dg1-dg2+dgW1-dgW2);

				for (int i = 0; i < K; i++) {
					gradientZW[i][z[i]] += gradientLL;
				}
				gradientW += gradientLL;
				gradientB += gradientLL;
			}

			for (int i = 0; i < K; i++) {
				for (int z = 0; z < Z[i]; z++) {
					gradientZW[i][z] += -(omegaZW[i][z][w] - etaZW[i][z][w]) / (sigma*sigma); 
					omegaZW[i][z][w] = omegaZW[i][z][w] + stepSizeW*gradientZW[i][z];
				}
			}

			gradientW += -(omegaW[w] - etaW[w]) / (sigma*sigma); 
			omegaW[w] = omegaW[w] + (stepSizeW)*gradientW;
		}
	
		gradientB += -omegaB / (sigmaWB*sigmaWB); 
		omegaB = omegaB + (stepSizeWB)*gradientB;
	}

	// the E and M steps, for one iteration
	public void doSampling(int iter) {
		// sample z values for all the tokens
		if (iter % blockFreq == 0) {
			for (int d = 0; d < D; d++) {
				for (int n = 0; n < docs[d].length; n++) {
					sample(d, n);
				}
			}
		} else {
			for (int d = 0; d < D; d++) {
				for (int n = 0; n < docs[d].length; n++) {
					sampleInd(d, n);
				}
			}
		}

		// 1 iteration of gradient ascent (change 1 to x to do x iterations)
		for (int i = 0; i < 1; i++) updateWeights(iter);

		// compute the priors with the new params and update the cached prior variables 
		alphaNorm = new double[D];
		omegaNorm = new double[Z0];
		for (int z = 0; z < Z0; z++) {
			for (int d = 0; d < D; d++) {
				priorDZ[d][z] = priorA(d, z);
				alphaNorm[d] += priorDZ[d][z];
				//if (d % 1000 == 0) System.out.println("alphaNorm"+d+" "+alphaNorm[d]);
			}

			for (int w = 0; w < W; w++) {
				priorZW[z][w] = priorW(w, z);
				omegaNorm[z] += priorZW[z][w];
			}
			System.out.println("omegaNorm"+z+" "+omegaNorm[z]);
		}

		if (iter % likelihoodFreq == 0) {
			System.out.println("Log-likelihood: "+computeLL());
		}

		// collect samples (docsZZ) 
		if (burnedIn) {
			for (int d = 0; d < D; d++) {
				for (int n = 0; n < docs[d].length; n++) { 
					int x = docsZ[d][n];
	
					docsZZ[d][n][x] += 1;
				}
			}
		}


	}

	public void sample(int d, int n) {
		int w = docs[d][n];
		int tuple = docsZ[d][n];
		
		// decrement counts

		nZW[tuple][w] -= 1;	
		nZ[tuple] -= 1;
		nDZ[d][tuple] -= 1;
	
		// sample new tuple value 

		double[] p = new double[Z0];
		double pTotal = 0;
	
		for (int z = 0; z < Z0; z++) {
			p[z] = (nDZ[d][z] + priorDZ[d][z]) *
				  (nZW[z][w] + priorZW[z][w]) / (nZ[z] + omegaNorm[z]);

			pTotal += p[z];
		}

		Random r = new Random();
		double u = r.nextDouble() * pTotal;
		
		double v = 0.0;
		for (int z = 0; z < Z0; z++) {
			v += p[z];
			
			if (v > u) {
				tuple = z;
				break;
			}
		}
		
		// increment counts

		nZW[tuple][w] += 1;	
		nZ[tuple] += 1;
		nDZ[d][tuple] += 1;
		
		// set new assignments

		docsZ[d][n] = tuple;
	}

	public void sampleInd(int d, int n) {
		int w = docs[d][n];
		int tuple = docsZ[d][n];
		
		// decrement counts

		nZW[tuple][w] -= 1;	
		nZ[tuple] -= 1;
		nDZ[d][tuple] -= 1;
	
		int[] z = iToVector(tuple);
		int[] zNew = new int[K];

		// sample new value for each factor
		for (int k = 0; k < K; k++) {
			double[] p = new double[Z[k]];
			double pTotal = 0;
	
			int[] zk = z.clone();

			for (int j = 0; j < Z[k]; j++) {
				zk[k] = j;
				int x = vectorToI(zk);	

				p[j] = (nDZ[d][x] + priorDZ[d][x]) *
					  (nZW[x][w] + priorZW[x][w]) / (nZ[x] + omegaNorm[x]);

				pTotal += p[j];
			}

			Random r = new Random();
			double u = r.nextDouble() * pTotal;
		
			double v = 0.0;
			for (int j = 0; j < Z[k]; j++) {
				v += p[j];
			
				if (v > u) {
					zNew[k] = j;
					break;
				}
			}
		}

		tuple = vectorToI(zNew);

		// increment counts

		nZW[tuple][w] += 1;	
		nZ[tuple] += 1;
		nDZ[d][tuple] += 1;
		
		// set new assignments

		docsZ[d][n] = tuple;
	}

	// computes the log-likelihood of the corpus
	// this marginalizes over the hidden variables but not the parameters
	// (i.e. we condition on the current MAP estimates of \theta and \phi)
	public double computeLL() {
		double LL = 0;

		for (int d = 0; d < D; d++) {
			for (int n = 0; n < docs[d].length; n++) { 
				int w = docs[d][n];

				double tokenLL = 0;

				// marginalize over z

				for (int z = 0; z < Z0; z++) {
					tokenLL += (nDZ[d][z] + priorDZ[d][z]) / (nD[d] + alphaNorm[d])*
				  		(nZW[z][w] + priorZW[z][w]) / (nZ[z] + omegaNorm[z]);
				}

				LL += Math.log(tokenLL);
			}
		}

		return LL;
	}

	public void readDocs(String filename) throws Exception {
		System.out.println("Reading input...");
		
		wordMap = new HashMap<String,Integer>();
		wordMapInv = new HashMap<Integer,String>();
		
		FileReader fr = new FileReader(filename);
		BufferedReader br = new BufferedReader(fr); 

		String s;
	
		D = 0;
		while((s = br.readLine()) != null) {
			D++;
		}

		docs = new int[D][];
		docsC = new int[D];

		fr = new FileReader(filename);
		br = new BufferedReader(fr); 

		int d = 0;
		while ((s = br.readLine()) != null) {
			String[] tokens = s.split("\\s+");
			
			int N = tokens.length;
			
			docs[d] = new int[N-1];
			docsC[d] = Integer.parseInt(tokens[0]); // the first token is a document ID

			for (int n = 1; n < N; n++) {
				String word = tokens[n];
				
				int key = wordMap.size();
				if (!wordMap.containsKey(word)) {
					wordMap.put(word, new Integer(key));
					wordMapInv.put(new Integer(key), word);
				}
				else {
					key = ((Integer) wordMap.get(word)).intValue();
				}
				
				docs[d][n-1] = key;
			}
			
			d++;
		}
		
		br.close();
		fr.close();
		
		W = wordMap.size();

		System.out.println(D+" documents");
		System.out.println(W+" word types");

		readPriors(priorPrefix);
	}

	public void readPriors(String prefix) throws Exception {
		etaW = new double[W];
		etaZW = new double[K][][];
		for (int k = 0; k < K; k++) {
			etaZW[k] = new double[Z[k]][W];
		}

		if (prefix == "") return;

		FileReader fr = new FileReader(prefix+".txt");
		BufferedReader br = new BufferedReader(fr); 

		String s;
		while ((s = br.readLine()) != null) {
			String[] tokens = s.split("\\s+");

			String word = tokens[0];
			if (!wordMap.containsKey(word)) continue;
			int w = wordMap.get(word);

			double weight = Double.parseDouble(tokens[1]);

			etaW[w] = weight;
		}

		br.close();
		fr.close();

		for (int k = 0; k < K; k++) {
			for (int z = 0; z < Z[k]; z++) {
				fr = new FileReader(prefix+k+"_"+z+".txt");
				br = new BufferedReader(fr); 

				while ((s = br.readLine()) != null) {
					String[] tokens = s.split("\\s+");

					String word = tokens[0];
					if (!wordMap.containsKey(word)) continue;
					int w = wordMap.get(word);

					double weight = Double.parseDouble(tokens[1]);
	
					etaZW[k][z][w] = weight;
				}

				br.close();
				fr.close();
			}
		}
	}


	public void writeOutput(String filename) throws Exception {
		FileWriter fw = new FileWriter(filename+".assign");
		BufferedWriter bw = new BufferedWriter(fw); 		

		for (int d = 0; d < D; d++) {
			bw.write(docsC[d]+" ");
			
			for (int n = 0; n < docs[d].length; n++) {
				String word = wordMapInv.get(docs[d][n]);

				//bw.write(word+":"+docsZ[d][n]+" "); // only current sample
				bw.write(word);  // for multiple samples
				for (int zz = 0; zz < Z0; zz++) {
					bw.write(":"+docsZZ[d][n][zz]);
				}
				bw.write(" ");
			}
			bw.newLine();
		}
		
		bw.close();
		fw.close();

		for (int i = 0; i < K; i++) {
			fw = new FileWriter(filename+".omegaZW"+i);
			bw = new BufferedWriter(fw); 		
			for (int w = 0; w < W; w++) {
				String word = wordMapInv.get(w);
				bw.write(word);
				for (int z = 0; z < Z[i]; z++) {
					bw.write(" "+omegaZW[i][z][w]);
				}
				bw.newLine();
			}
			bw.close();
			fw.close();
		}

		fw = new FileWriter(filename+".omegaW");
		bw = new BufferedWriter(fw); 		
		for (int w = 0; w < W; w++) {
			String word = wordMapInv.get(w);
			bw.write(word);
			bw.write(" "+omegaW[w]);
			bw.newLine();
		}
		bw.close();
		fw.close();

		fw = new FileWriter(filename+".omegaB");
		bw = new BufferedWriter(fw); 		
		bw.write(""+omegaB);
		bw.close();
		fw.close();

		for (int i = 0; i < K; i++) {
			fw = new FileWriter(filename+".alphaZ"+i);
			bw = new BufferedWriter(fw); 		
			for (int z = 0; z < Z[i]; z++) {
				bw.write(alphaZ[i][z]+"");
				bw.newLine();
			}
			bw.close();
			fw.close();
		}

		fw = new FileWriter(filename+".alphaB");
		bw = new BufferedWriter(fw); 		
		bw.write(""+alphaB);
		bw.close();
		fw.close();

		fw = new FileWriter(filename+".beta");
		bw = new BufferedWriter(fw); 		
		for (int x = 0; x < Z0; x++) {
			bw.write(""+logistic(beta[x]));
			bw.newLine();
		}
		bw.close();
		fw.close();
	}

	// Approximation to the digamma function, from Radford Neal.
	// can also use Gamma.digamma() from commons
	public static double digamma0(double x)
	{
		//return Gamma.digamma(x);

		double r = 0.0;

		while (x <= 5.0) {
			r -= 1.0 / x;
			x += 1.0;
		}

		double f = 1.0 / (x * x);
		double t = f * (-1 / 12.0 + f * (1 / 120.0 + f * (-1 / 252.0 + f * (1 / 240.0 + f * (-1 / 132.0 + f * (691 / 32760.0 + f * (-1 / 12.0 + f * 3617.0 / 8160.0)))))));
		return r + Math.log(x) - 0.5 / x + t;
	}
}
