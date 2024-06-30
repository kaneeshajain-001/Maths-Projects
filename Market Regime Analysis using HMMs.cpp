#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
// #include <chrono>

using namespace std;
// using namespace std::chrono; 

double f(long double y, long double mean, long double var) {
    long double num = exp(-0.5*pow((y - mean) / sqrt(var), 2));
    long double den = sqrt(2*M_PI*var);
    return num / den;
}


int main() {
	long double epsilon = pow(10,-8);
	
    vector<vector<long double>> A(2, vector<long double>(2));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            cin >> A[i][j];
        }
    }
    
    vector<long double> pi(2);
    cin >> pi[0] >> pi[1];
    
    int T;
    cin >> T;
    
    vector<long double> Y(T);
    for (int i = 0; i < T; ++i) {
        cin >> Y[i];
    }
	
	//auto start = high_resolution_clock::now();

    vector<long double> mean = {0, 0};
    vector<long double> var = {1, 1};

    vector<vector<long double>> alpha(2, vector<long double>(T));
    vector<vector<long double>> beta(2, vector<long double>(T));
    vector<vector<long double>> gamma(2, vector<long double>(T));
    vector<vector<vector<long double>>> e(2, vector<vector<long double>>(2, vector<long double>(T - 1)));
    vector<vector<long double>> prob(2, vector<long double>(T));
	
	vector<vector<long double>> old_A = A;
	vector<long double> old_pi = pi;
	vector<long double> old_mean = mean;
	vector<long double> old_var = var;
	
	for (int i = 0; i < T; ++i) {
        prob[0][i] = f(Y[i], mean[0], var[0]);
        prob[1][i] = f(Y[i], mean[1], var[1]);
    }
	
    for (int i = 0; i < 2; ++i) {
        alpha[i][0] = pi[i]*prob[i][0];
    }

    for (int t = 1; t < T; ++t) {
        for (int i = 0; i < 2; ++i) {
            long double sum_all = 0;
            for (int j = 0; j < 2; ++j) {
                sum_all += alpha[j][t - 1] * A[j][i];
            }
            alpha[i][t] = prob[i][t]*sum_all;
        }
    }

    for (int i = 0; i < 2; ++i) {
        beta[i][T - 1] = 1;
    }

    for (int t = T - 2; t >= 0; --t) {
        for (int i = 0; i < 2; ++i) {
            beta[i][t] = 0;
            for (int j = 0; j < 2; ++j) {
                beta[i][t] += beta[j][t+1]*A[i][j]*prob[j][t+1];
            }
        }
    }

    for (int t = 0; t < T; ++t) {
        long double sum_all = 0;
        for (int j = 0; j < 2; ++j) {
            sum_all += alpha[j][t]*beta[j][t];
        }
        for (int i = 0; i < 2; ++i) {
            gamma[i][t] = alpha[i][t]*beta[i][t] / sum_all;
        }
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int t = 0; t < T - 1; ++t) {
                long double num = alpha[i][t]*A[i][j]*beta[j][t+1]*prob[j][t+1];
                long double den = 0;
                for (int k = 0; k < 2; ++k) {
                    for (int w = 0; w < 2; ++w) {
                        den += alpha[k][t]*A[k][w]*beta[w][t+1]*prob[w][t+1];
                    }
                }
                e[i][j][t] = num / den;
            }
        }
    }

    for (int i = 0; i < 2; ++i) {
        pi[i] = gamma[i][0];
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            long double num = 0;
            long double den = 0;
            for (int t = 0; t < T - 1; ++t) {
                num += e[i][j][t];
                den += gamma[i][t];
            }
            A[i][j] = num / den;
        }
    }

    for (int i = 0; i < 2; ++i) {
        long double num = 0;
        long double den = 0;
        for (int t = 0; t < T; ++t) {
            num += Y[t] * gamma[i][t];
            den += gamma[i][t];
        }
        mean[i] = num / den;
    }

    for (int i = 0; i < 2; ++i) {
        long double num = 0;
        long double den = 0;
        for (int t = 0; t < T; ++t) {
            num += pow((Y[t] - mean[i]), 2)*gamma[i][t];
            den += gamma[i][t];
        }
        var[i] = num / den;
    }
	
    long double max_diff_A = 0;
	long double max_diff_pi = 0;
	long double max_diff_mean = 0;
	long double max_diff_var = 0;

    
    for (int i=0; i<2; ++i) {
    	for (int j=0; j<2; ++j) {
    		if (abs(old_A[i][j] - A[i][j]) > max_diff_A) {
    			max_diff_A = abs(old_A[i][j] - A[i][j]);
			}
		}
	}
	
	for (int i=0; i<2; ++i) {
		if (abs(old_pi[i] - pi[i]) > max_diff_pi) {
			max_diff_pi = abs(old_pi[i] - pi[i]);
		}
		if (abs(old_mean[i] - mean[i]) > max_diff_mean) {
			max_diff_mean = abs(old_mean[i] - mean[i]);
		}
		if (abs(old_var[i] - var[i]) > max_diff_var) {
			max_diff_var = abs(old_var[i] - var[i]);
		}
	}
	
	/*cout<<"first iteration"<<endl;
    cout<<"mean diff"<<" "<<max_diff_mean<<endl;
    cout<<"var diff"<<" "<<max_diff_var<<endl;
    cout<<"a diff"<<" "<<max_diff_A<<endl;
    cout<<"pi diff"<<" "<<max_diff_pi<<endl;*/
    
    
	int loop = 0;
    while(max_diff_mean >= epsilon or max_diff_var >= epsilon or max_diff_A >= epsilon or max_diff_pi >= epsilon) {
    	loop++;
		old_A = A;
		old_pi = pi;
		old_mean = mean;
		old_var = var;
		
		for (int i = 0; i < T; ++i) {
        	prob[0][i] = f(Y[i], mean[0], var[0]);
        	prob[1][i] = f(Y[i], mean[1], var[1]);
    	}
		
	    for (int i = 0; i < 2; ++i) {
	        alpha[i][0] = pi[i]*prob[i][0];
	    }
	
	    for (int t = 1; t < T; ++t) {
	        for (int i = 0; i < 2; ++i) {
	            long double sum_all = 0;
	            for (int j = 0; j < 2; ++j) {
	                sum_all += alpha[j][t - 1] * A[j][i];
	            }
	            alpha[i][t] = prob[i][t]*sum_all;
	        }
	    }
	
	    for (int i = 0; i < 2; ++i) {
	        beta[i][T - 1] = 1;
	    }
	
	    for (int t = T - 2; t >= 0; --t) {
	        for (int i = 0; i < 2; ++i) {
	            beta[i][t] = 0;
	            for (int j = 0; j < 2; ++j) {
	                beta[i][t] += beta[j][t+1]*A[i][j]*prob[j][t+1];
	            }
	        }
	    }
	
	    for (int t = 0; t < T; ++t) {
	        long double sum_all = 0;
	        for (int j = 0; j < 2; ++j) {
	            sum_all += alpha[j][t]*beta[j][t];
	        }
	        for (int i = 0; i < 2; ++i) {
	            gamma[i][t] = alpha[i][t]*beta[i][t] / sum_all;
	        }
	    }
	
	    for (int i = 0; i < 2; ++i) {
	        for (int j = 0; j < 2; ++j) {
	            for (int t = 0; t < T - 1; ++t) {
	                long double num = alpha[i][t]*A[i][j]*beta[j][t+1]*prob[j][t+1];
	                long double den = 0;
	                for (int k = 0; k < 2; ++k) {
	                    for (int w = 0; w < 2; ++w) {
	                        den += alpha[k][t]*A[k][w]*beta[w][t+1]*prob[w][t+1];
	                    }
	                }
	                e[i][j][t] = num / den;
	            }
	        }
	    }
	
	    for (int i = 0; i < 2; ++i) {
	        pi[i] = gamma[i][0];
	    }
	
	    for (int i = 0; i < 2; ++i) {
	        for (int j = 0; j < 2; ++j) {
	            long double num = 0;
	            long double den = 0;
	            for (int t = 0; t < T - 1; ++t) {
	                num += e[i][j][t];
	                den += gamma[i][t];
	            }
	            A[i][j] = num / den;
	        }
	    }
	
	    for (int i = 0; i < 2; ++i) {
	        long double num = 0;
	        long double den = 0;
	        for (int t = 0; t < T; ++t) {
	            num += Y[t] * gamma[i][t];
	            den += gamma[i][t];
	        }
	        mean[i] = num / den;
	    }
	
	    for (int i = 0; i < 2; ++i) {
	        long double num = 0;
	        long double den = 0;
	        for (int t = 0; t < T; ++t) {
	            num += pow((Y[t] - mean[i]), 2)*gamma[i][t];
	            den += gamma[i][t];
	        }
	        var[i] = num / den;
	    }
		
	    max_diff_A = 0;
		max_diff_pi = 0;
		max_diff_mean = 0;
		max_diff_var = 0;
	    
	    for (int i=0; i<2; ++i) {
	    	for (int j=0; j<2; ++j) {
	    		if (abs(old_A[i][j] - A[i][j]) > max_diff_A) {
	    			max_diff_A = abs(old_A[i][j] - A[i][j]);
				}
			}
		}
		
		for (int i=0; i<2; ++i) {
			if (abs(old_pi[i] - pi[i]) > max_diff_pi) {
				max_diff_pi = abs(old_pi[i] - pi[i]);
			}
			if (abs(old_mean[i] - mean[i]) > max_diff_mean) {
				max_diff_mean = abs(old_mean[i] - mean[i]);
			}
			if (abs(old_var[i] - var[i]) > max_diff_var) {
				max_diff_var = abs(old_var[i] - var[i]);
			}
		}
    }
    
    //cout<<"Loop: "<<loop<<endl;
    
    
	for (int t = 0; t < T; ++t) {
        if (mean[0] > mean[1]) {
            if (gamma[0][t] > gamma[1][t]) {
                cout << "Bull" << endl;
            } else {
                cout << "Bear" << endl;
            }
        } else {
            if (gamma[0][t] > gamma[1][t]) {
                cout << "Bear" << endl;
            } else {
                cout << "Bull" << endl;
            }
        }
    }
    
	/*auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<milliseconds>(stop - start); 

    cout << "Execution time: " << duration.count() << " milliseconds" << endl; */
    
    return 0;

}

