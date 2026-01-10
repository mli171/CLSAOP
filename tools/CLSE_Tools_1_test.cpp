#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double pmvnorm_cpp(arma::vec& lb, arma::vec& ub, arma::vec& mu, arma::vec& lowertrivec, double abseps = 1e-3){
  
  int n = lb.n_elem;
  int nu = 0;
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate
  
  double* lb_ = lb.memptr();                        // lower bound;
  double* ub_ = ub.memptr();                        // upper bound;
  double* correlationMatrix = lowertrivec.memptr(); // array of correlation coefficients;
  int* infin = new int[n];                          // Integer, array of integration limits flag;
  double* mu_ = mu.memptr();                        // array of non-centrality parameters;
  
  // if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  // if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  // if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  // if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)];
  
  for (int i = 0; i < n; ++i) {
    if(lb(i) == R_NegInf){
      infin[i] = 0;
    }else if (ub(i) == R_PosInf){
      infin[i] = 1;
    }else{
      infin[i] = 2;
    }
  }
  
  // return values
  double error;
  double value;
  int inform;
  
  // Fortran function details check:
  // https://github.com/cran/mvtnorm/blob/master/src/mvt.f
  mvtnorm_C_mvtdst(&n, &nu, lb_, ub_,
                   infin, correlationMatrix, mu_,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   delete[] (infin);
                   
                   return value;
}

// // [[Rcpp::export]]
// arma::mat multikappa_rcpp(int i, int j,
//                           double msti, double mstj,
//                           arma::vec EXwi, arma::vec EXwj,
//                           arma::vec seg, double rho) {
//   int K = seg.n_elem - 1;
//   arma::mat res(K - 1, K - 1);
//   
//   int lag = std::abs(i - j);
//   double corrval = std::pow(rho, lag);
//   // NumericMatrix corr(2, 2);
//   // corr(0, 0) = 1.0;
//   // corr(1, 1) = 1.0;
//   // corr(0, 1) = corrval;
//   // corr(1, 0) = corrval;
//   // NumericVector mean = NumericVector::create(msti, mstj);
//   // Function pmvnorm("pmvnorm", Environment::namespace_env("mvtnorm"));
//   
//   double p;
//   arma::vec mean = {msti, mstj};
//   arma::vec lb(2), ub(2), correle(1);
//   correle(0) = corrval;
//   
//   for (int k = 0; k < K - 1; ++k) {
//     for (int kp = 0; kp < K - 1; ++kp) {
//       lb = {seg(k), seg(kp)};
//       ub = {seg(k+1), seg(kp+1)};
//       res(k, kp) = pmvnorm_cpp(lb, ub, mean, correle, 1e-3);
//       // NumericVector lower = NumericVector::create(seg(k), seg(kp));
//       // NumericVector upper = NumericVector::create(seg(k + 1), seg(kp + 1));
//       // NumericVector p = pmvnorm(_["lower"] = lower,
//       //                           _["upper"] = upper,
//       //                           _["mean"] = mean,
//       //                           _["corr"] = corr);
//       // res(k, kp) = p[0];
//     }
//   }
//   
//   res = res - EXwi * EXwj.t();
//   // Rcout << "Mo 1 = " << res(k, kp) << std::endl;
//   return res;
// }

// [[Rcpp::export]]
arma::mat multikappa_rcpp(int i, int j,
                          double msti, double mstj,
                          arma::vec EXwi, arma::vec EXwj,
                          arma::vec seg, double rho) {
  int K = seg.n_elem - 1;
  arma::mat res(K - 1, K - 1);
  
  int lag = std::abs(i - j);
  double corrval = std::pow(rho, lag);
  NumericMatrix corr(2, 2);
  corr(0, 0) = 1.0;
  corr(1, 1) = 1.0;
  corr(0, 1) = corrval;
  corr(1, 0) = corrval;
  NumericVector mean = NumericVector::create(msti, mstj);
  Function pmvnorm("pmvnorm", Environment::namespace_env("mvtnorm"));
  
  for (int k = 0; k < K - 1; ++k) {
    for (int kp = 0; kp < K - 1; ++kp) {
      NumericVector lower = NumericVector::create(seg(k), seg(kp));
      NumericVector upper = NumericVector::create(seg(k + 1), seg(kp + 1));
      NumericVector p = pmvnorm(_["lower"] = lower,
                                _["upper"] = upper,
                                _["mean"] = mean,
                                _["corr"] = corr);
      res(k, kp) = p[0];
    }
  }
  
  res = res - EXwi * EXwj.t();
  return res;
}

// [[Rcpp::export]]
List multiInnov_cpp(NumericVector mst,
                    arma::mat EXW,
                    arma::vec seg,
                    double rho,
                    int numcoef,
                    int mylag=0,
                    double mytol = 1e-6) {
  
  int Ts = mst.size();
  int K = seg.n_elem - 1;
  
  std::vector<arma::mat> V(Ts);
  arma::cube HQ(K - 1, K - 1, Ts * numcoef, fill::zeros);
  
  //--------------- t=0 ---------------//
  arma::vec EXW0 = EXW.row(0).t();
  V[0] = multikappa_rcpp(1, 1, mst[0], mst[0], EXW0, EXW0, seg, rho);
  
  //--------------- t=1 ---------------//
  arma::vec EXW1 = EXW.row(1).t();
  HQ.slice(0) = multikappa_rcpp(2, 1, mst[1], mst[0], EXW1, EXW0, seg, rho) * arma::inv(V[0]); //k=0
  arma::mat tmp = HQ.slice(0) * V[0] * HQ.slice(0).t();
  V[1] = multikappa_rcpp(2, 2, mst[1], mst[1], EXW1, EXW1, seg, rho) - tmp;
  
  //--------------- t>=2 ---------------//
  // int mylag = 0;
  arma::mat tmpcheck = arma::mat(K - 1, K - 1).fill(100.0);
  for (int t = 2; t < Ts; ++t) {
    arma::vec EXWt = EXW.row(t).t();
    if (arma::any(arma::vectorise(tmpcheck) > mytol) && mylag == 0) {
      // k = 0
      HQ.slice((t-1)*numcoef + t-1) = multikappa_rcpp(t + 1, 1, mst[t], mst[0], EXWt, EXW0, seg, rho) * arma::inv(V[0]);
      // Rcpp::Rcout << "HQ.slice((t-1)*numcoef + t-1) = " << HQ.slice((t-1)*numcoef + t-1) << std::endl;
      
      // k > 0
      for (int k = 1; k < t; ++k) {
        arma::mat acc(K - 1, K - 1, arma::fill::zeros);
        for (int j = 0; j < k; ++j) {
          acc += HQ.slice((t-1)*numcoef + t -j-1) * V[j] * HQ.slice((k-1)*numcoef + k-j-1).t();
        }
        arma::vec EXWk = EXW.row(k).t();
        HQ.slice((t-1)*numcoef + t-k-1) = (multikappa_rcpp(t + 1, k + 1, mst[t], mst[k], EXWt, EXWk, seg, rho) - acc) * arma::inv(V[k]);
        // Rcpp::Rcout << "HQ.slice((t-1)*numcoef + t-k-1) = " << HQ.slice((t-1)*numcoef + t-k-1) << std::endl;
      }
      // Prediction MSE
      arma::mat acc(K - 1, K - 1, arma::fill::zeros);
      for (int j = 0; j < t; ++j) {
        acc += HQ.slice((t-1)*numcoef + t-j-1) * V[j] * HQ.slice((t-1)*numcoef + t-j-1).t();
      }
      V[t] = multikappa_rcpp(t + 1, t + 1, mst[t], mst[t], EXWt, EXWt, seg, rho) - acc;
      
      tmpcheck = arma::abs(HQ.slice((t-1-1)*numcoef + t-1-1));
      if (HQ.slice((t-1)*numcoef + t-1).max() < mytol){
        // set lags to save computation time
        mylag = t;
      }
    } else {
      if(t-mylag < 0){
        // k = 0
        HQ.slice((t-1)*numcoef + t-1) = multikappa_rcpp(t + 1, 1, mst[t], mst[0], EXWt, EXW0, seg, rho) * arma::inv(V[0]);
        // Rcpp::Rcout << "HQ.slice((t-1)*numcoef + t-1) = " << HQ.slice((t-1)*numcoef + t-1) << std::endl;
        
        // k > 0
        for (int k = 1; k < t; ++k) {
          arma::mat acc(K - 1, K - 1, arma::fill::zeros);
          for (int j = 0; j < k; ++j) {
            acc += HQ.slice((t-1)*numcoef + t -j-1) * V[j] * HQ.slice((k-1)*numcoef + k-j-1).t();
          }
          arma::vec EXWk = EXW.row(k).t();
          HQ.slice((t-1)*numcoef + t-k-1) = (multikappa_rcpp(t + 1, k + 1, mst[t], mst[k], EXWt, EXWk, seg, rho) - acc) * arma::inv(V[k]);
          // Rcpp::Rcout << "HQ.slice((t-1)*numcoef + t-k-1) = " << HQ.slice((t-1)*numcoef + t-k-1) << std::endl;
        }
        // Prediction MSE
        arma::mat acc(K - 1, K - 1, arma::fill::zeros);
        for (int j = 0; j < t; ++j) {
          acc += HQ.slice((t-1)*numcoef + t-j-1) * V[j] * HQ.slice((t-1)*numcoef + t-j-1).t();
        }
        V[t] = multikappa_rcpp(t + 1, t + 1, mst[t], mst[t], EXWt, EXWt, seg, rho) - acc;
        
        tmpcheck = arma::abs(HQ.slice((t-1-1)*numcoef + t-1-1));
      }else{
        // Only to cover the effective coefficients (lag), the rest are less than 1e-05
        for (int k=t-mylag; k<t; ++k) {
          // for (int k=t-mylag; k<6; ++k) {
          //  Rcpp::Rcout << "====================== k=" << k << std::endl;
          arma::mat acc(K - 1, K - 1, arma::fill::zeros);
          for (int j = std::max(0, k - mylag); j < k; ++j) {
            if (t-j < mylag) {
              // Rcpp::Rcout << "t=" << t << std::endl;
              // Rcpp::Rcout << "t-j=" << t-j << std::endl;
              // Rcpp::Rcout << "HQ.slice((t-1)*numcoef + t-j-1)=" << HQ.slice((t-1)*numcoef + t-j-1)*1000000000 << std::endl;
              // Rcpp::Rcout << "HQ.slice((k-1)*numcoef + k-j-1)=" << HQ.slice((k-1)*numcoef + k-j-1) << std::endl;
              acc += HQ.slice((t-1)*numcoef + t-j-1) * V[j] * HQ.slice((k-1)*numcoef + k-j-1).t();
              // Rcpp::Rcout << "acc=" << acc << std::endl;
            }
          }
          arma::vec EXWk = EXW.row(k).t();
          HQ.slice((t-1)*numcoef + t-k-1) = (multikappa_rcpp(t + 1, k + 1, mst[t], mst[k], EXWt, EXWk, seg, rho) - acc) * arma::inv(V[k]);
        }
        
        arma::mat acc(K - 1, K - 1, arma::fill::zeros);
        for (int j = std::max(0, t - mylag); j < t; ++j) {
          if (t - j < mylag) {
            acc += HQ.slice((t-1)*numcoef + t-j-1) * V[j] * HQ.slice((t-1)*numcoef + t-j-1).t();
          }
        }
        V[t] = multikappa_rcpp(t + 1, t + 1, mst[t], mst[t], EXWt, EXWt, seg, rho) - acc;
        
        if ((t - 1 < Ts) && (mylag > 0) && (mylag <= numcoef)) {
          tmpcheck = arma::abs(HQ.slice((t - 1) * numcoef + mylag - 1));
        }
      }
    }
  }
  
  return List::create(Named("HQ") = HQ,
                      Named("V") = V,
                      Named("mylag") = mylag);
}



// [[Rcpp::export]]
Rcpp::List MultiOneStepPred_cpp(const arma::mat& EXW,
                                const arma::mat& Xw,
                                const arma::cube& HQ,
                                int numcoef,
                                int mylag) {
  
  int Ts = Xw.n_rows;
  int K1 = Xw.n_cols;
  
  arma::mat Err = Xw - EXW;
  arma::mat PredErr(Ts, K1, arma::fill::zeros);
  
  // for (int t = 0; t < Ts - 1; ++t) {
  //   int predlag = (mylag > 0) ? mylag : t;
  //   
  //   for (int j = 1; j <= std::min(t + 1, predlag); ++j) {
  //     arma::vec diff = Err.row(t + 1 - j).t() - PredErr.row(t + 1 - j).t();
  //     arma::mat HQtj = HQ.slice(t * numcoef + j - 1);
  //     PredErr.row(t + 1) += (HQtj * diff).t();
  //   }
  // }
  
  for (int t = 1; t < Ts - 1; ++t) {
    int predlag = (mylag > 0) ? mylag : t;
    int block = (t - 1) * numcoef;
    
    for (int j = 1; j <= std::min(t + 1, predlag); ++j) {
      arma::vec diff = Err.row(t + 1 - j).t() - PredErr.row(t + 1 - j).t();
      arma::mat HQtj = HQ.slice(block + j - 1);
      PredErr.row(t + 1) += (HQtj * diff).t();
    }
  }
  
  arma::mat ErrErr = Err - PredErr;
  arma::mat PredX = PredErr + EXW;
  
  return Rcpp::List::create(Rcpp::Named("ErrErr") = ErrErr,
                            Rcpp::Named("PredX") = PredX);
}