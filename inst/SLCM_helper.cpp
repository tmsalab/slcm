#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

// input: 1. K number of classes; 2. c, class index;  
// output: attribute profile, e.g. 3 = (0,1,1)
arma::mat bijectionvector(unsigned int K, arma::vec c) {
  unsigned int N = c.size();
  arma::mat vv=arma::zeros<arma::mat>(N,K);
  unsigned int tmp;
  for (unsigned int n = 0;n<N;n++){
    tmp = c(n);
    for(unsigned int k=0;k<K;k++){
      if(tmp/pow(2,K-k-1)>=1){vv(n,k)=1;}
      tmp = tmp - vv(n,k)*pow(2,K-k-1);
      
    }
  }
  return vv;
}

// input: 1. index, list of indexes; 2. P, classes : 0 - (P-1);  
// output: number of indexes in each class - 1, e.g. c(1,2,3), 4 ==> (-1, 0, 0, 0).
arma::vec tableC(arma::ivec index, unsigned int P) {
  arma::vec counts = -1*arma::ones<arma::vec>(P);
  
  unsigned int n = index.size();
  
  for (int i = 0; i < n; i++) {
    counts(index(i))++;
  }
  
  return counts;
}



double nek_check(arma::vec count, unsigned int K){
  count.shed_row(0);
  arma::vec row_sum = 2*arma::ones<arma::vec>(K) - count.rows(0,K-1);
  row_sum.elem(arma::find(row_sum < 0) ).fill(0.0);
  if(max(row_sum) < 1){
    return 1;
  }
  for (unsigned int a1=0; a1< std::min(count(K),row_sum(0))+1; a1++){
    unsigned int a2 = count(K)-a1;
    for (unsigned int b1=0; b1< std::min(count((K+1)),row_sum(0)-a1)+1; b1++){
      unsigned int b2 = count((K+1)) - b1;
      for (unsigned int c1=0; c1 < std::min(count((K+2)),row_sum(1)-a2)+1; c1++){
        unsigned int c2 = count((K+2)) - c1;
        unsigned int d1 = row_sum(0) - a1 - b1;
        if(d1>count((pow(2,K)-2))){
          continue;
        }
        int tmp = row_sum(1)-a2-c1;
        int tmp1 = 0;
        int d2 = std::max(tmp,tmp1);
        if(d2>(count(pow(2,K)-2)-d1)){
          continue;
        }
        unsigned int d3 = count((pow(2,K)-2))-d1-d2;
        if((b2+c2+d3)>(row_sum(2)-1)){
          return 1;
        }
      }
    }
  }
  return 0;
}


double check_D(arma::mat Delta, arma::mat k3){
  arma::mat D = Delta;
  unsigned int P = D.n_cols+1;
  unsigned int flag = 0;
  unsigned int K = log(P)/log(2);
  arma::mat D1 = D.cols(0,K-1);
  arma::vec col_sum = sort(arma::sum(D1,0).t(),"descend");
  
  if(min(col_sum)<3){
    return 0;
  }
  
  arma::vec tmp1 = 2*arma::linspace<arma::vec>(K,1,K)+1;
  arma::vec diff = col_sum - tmp1;
  if(min(diff)>-1){
    return 1.0;
  }
  
  k3 = join_rows(arma::zeros<arma::vec>(K),k3);
  D1 = join_cols(D1,k3.t());
  
  arma::vec tmp2 = arma::logspace(log(2.0)/log(10.0)*(K-1), 0, K);
  arma::ivec tmp3 = arma::conv_to<arma::ivec>::from(round(tmp2));
  arma::ivec index = arma::conv_to<arma::ivec>::from(round(D1*tmp2));
  
  arma::ivec ref = arma::conv_to<arma::ivec>::from(round(k3.t()*tmp2));
  arma::vec count1 = tableC(index,P);
  arma::uvec uref=arma::conv_to<arma::uvec>::from(ref);
  arma::vec count = count1.rows(uref);
  
  if(count((pow(2,K)-1))>0){
    count((pow(2,K)-1)) = count((pow(2,K)-1))-1;
    return(nek_check(count,K));
  }
  
  for (unsigned int l = 0; l<K; l++){
    unsigned int k = tmp3(l);
    
    if((count1(k)*count1((P-k-1)))>0){
      arma::vec count2 = count1;
      count2(k) = count1(k)-1;
      count2((P-k-1)) = count2((P-k-1))-1;
      flag = nek_check(count2.rows(uref),K);
      if(flag==1){
        return(flag);
        break;
      }
    }
  }
  
  for (unsigned int i =0; i<(K-1); i++){
    for (unsigned int j = i+1; j<K; j++){
      if((count((K+i+1))*count((K+j+1))>0)){
        arma::vec count2 = count;
        count2((K+i+1)) = count((K+i+1)) -1;
        count2((K+j+1)) = count((K+j+1)) -1;
        flag = nek_check(count2,K);
        if(flag==1){
          return(flag);
          break;
        }
        
      }
    }
  }
  for (unsigned int k=0; k<K; k++){
    count((k+1))--;
  }
  flag = nek_check(count,K);
  
  return(flag);
}


// check if the Delta matrix satisfy the generically identifiability condition
double gen_identify_check(const arma::mat Delta, arma::mat k3){
  
  arma::mat D = Delta;
  double min_c = check_D(D,k3);
  arma::vec r_sum = arma::sum(D,1);
  double min_r = (arma::min(r_sum)>0);
  
  return (min_c+min_r>1);
}


// check if the Delta matrix satisfy the strict identifiability condition
double identify_check(const arma::mat Delta){
  
  arma::mat D = Delta;
  unsigned int K = log(D.n_cols+1)/log(2.0);
  
  // make a matrix to check the identity matrix
  arma::mat ones_zero_on_diag = -1*arma::ones<arma::mat>(K,K);
  ones_zero_on_diag.diag().fill(0);
  
  // find the rows with row sum 0 in the interation terms
  arma::vec row_sum = (arma::sum(D.cols(K,pow(2,K)-2),1));
  arma::mat Dsub = D.rows(find(row_sum==0));
  
  // check if 2 identity matrices exist
  arma::mat I_check = Dsub.cols(0,K-1)*ones_zero_on_diag;
  arma::mat I_count = arma::zeros<arma::mat>(Dsub.n_rows,K);
  I_count.elem(arma::find(I_check > -1) ).fill(1.0);
  arma::vec n_ek = (arma::sum(I_count,0)).t();
  
  // check column sum and row sum
  arma::vec c_sum = (arma::sum(D.cols(0,K-1),0)).t();
  arma::vec r_sum = arma::sum(D,1);
  
  
  double min_c = (arma::min(c_sum)>2);
  double min_r = (arma::min(r_sum)>0);
  double min_ek = (arma::min(n_ek)>1);
  
  return (min_c+min_r+min_ek>2);
}



// one random sample normal(mu, sigma)
// [[Rcpp::export]]
arma::vec  myrnorm(arma::vec mu, arma::vec sigma) {
  unsigned int n = mu.size(); 
  arma::vec a = as<arma::vec>(Rcpp::rnorm(n))%sigma + mu;
  return (a);
}


//one random sample multinomial(ps)
// output: from 0 to class-1
unsigned int rmultinomial(const arma::vec& ps){
  unsigned int C = ps.n_elem;
  arma::vec pps = ps/sum(ps);
  double u = R::runif(0,1);
  arma::vec cps = cumsum(pps);
  arma::vec Ips = arma::zeros<arma::vec>(C);
  Ips.elem(arma::find(cps < u) ).fill(1.0);
  
  return sum(Ips);
}


//a list of random sample bernoulli(prob)
// output: from 0 to 1
arma::vec rbernoulli (arma::vec prob){
  unsigned int n = prob.n_elem;
  Rcpp::NumericVector u(n);
  u = Rcpp::runif(n); 
  arma::vec z = arma::zeros<arma::vec>(n);
  
  for (unsigned int i = 0;i < n;i++){
    if(u(i) < prob(i))
    {z(i) = 1;}
    else 
    {z(i) = 0;}
  }
  
  return z;
}

// one random sample bernoulli(prob)
// output: from 0 to 1
int rbernoulli_one (double prob){
  double u;
  u = R::runif(0,1); 
  int z;
  if(u < prob){
    z = 1;}
  else{z = 0;}
  
  return z;
}





// one random sample of dirichlet(deltas)
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec& deltas){
  unsigned int C = deltas.n_elem;
  arma::vec Xgamma(C);
  
  //generating gamma(deltac,1)
  for(unsigned int c = 0;c < C;c++){
    Xgamma(c) = R::rgamma(deltas(c),1.0);
  }
  return Xgamma/sum(Xgamma);
}


// random restricted normal according to y: if y = 0, Z < c, if y = 1, Z > c;
//input: lists of mu, sigma, y, c (cutoff)
// [[Rcpp::export]]
arma::vec res_norm(arma::vec mu, arma::vec sigma, arma::vec y,arma::vec c){
  unsigned int n = mu.size();
  Rcpp::NumericVector Y = wrap(y);
  Rcpp::NumericVector u = Rcpp::runif(n);
  Rcpp::NumericVector sd_mu = wrap((c - mu) / sigma);
  Rcpp::NumericVector r = Rcpp::pnorm(sd_mu); // quantile of c
  Rcpp::NumericVector z = Rcpp::qnorm((r + u*(1-r)) * Y + u * r * (1-Y));
  arma::vec Z =  as<arma::vec>(z);
  Z = Z%sigma + mu;
  Z(find(Z > 999 )) = c(find(Z > 999));
  Z(find(Z < -999))= c(find(Z <- 999));
  return Z;}


// one sample of random restricted normal according to y: if y = 0, Z < c, if y = 1, Z > c;
//input: mu, sigma, y, c (cutoff)
double res_norm1(double mu,double sigma, double y,double c){
  Rcpp::NumericVector Y = wrap(y);
  Rcpp::NumericVector u = Rcpp::runif(1);
  Rcpp::NumericVector sd_mu = wrap((c-mu)/sigma);
  Rcpp::NumericVector r = Rcpp::pnorm(sd_mu);
  Rcpp::NumericVector z = Rcpp::qnorm((r+u*(1-r))*Y+u*r*(1-Y));
  double Z = z(0);
  Z = Z * sigma + mu;
  if(Z>999){Z = c;}
  if(Z<-999){Z = c;}
  return Z;}


// get design matrix from attribute profiles: 
// input: terms, the mapping from attribute profile to design matrix
// [[Rcpp::export]]
arma::mat  des_mat(arma::mat a, arma::mat terms) {
  unsigned int K = a.n_cols;
  unsigned int N = a.n_rows;
  unsigned int P = pow(2,K);
  
  arma::mat A = arma::ones<arma::mat>(N,P);
  
  
  for (unsigned int p = 1;p < P;p++){
    arma::vec index = terms.col(p-1);
    arma::uvec non_zero_ind = find(index);
    
    unsigned int n = non_zero_ind.size();
    
    for (unsigned int j=0;j<n;j++){
      A.col(p) = A.col(p)%a.col(non_zero_ind(j));
    }
  }
  
  return A;
}

// n random sample of multivarite normal (mu, sigma)
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// change Delta to Q (no intercept)
// input : delta , terms, output: Q
arma::mat toQ(arma::mat D, arma::mat terms){
  arma::mat absD = abs(D);
  arma::mat Q = absD*terms.t();
  Q.elem(find(Q>0)).fill(1.0);
  return Q;}


// change beta to beta (with intercept)
// input : beta, terms, output: prob (ordered as terms)
arma::mat toProb(arma::mat beta, arma::mat terms){
  unsigned int P = beta.n_cols;
  unsigned int J = beta.n_rows;
  unsigned int K = log(P)/log(2);
  arma::mat alpha = bijectionvector(K,arma::linspace<arma::vec>(0,(P-1),P));
  arma::mat O = des_mat(alpha,terms);
  Rcpp::NumericMatrix mu = wrap(O*beta.t());
  arma::mat theta(P,J);
  for (unsigned int j=0;j < J;j++){
    Rcpp::NumericVector tmp = Rcpp::pnorm(mu(_,j));
    theta.col(j) = as<arma::vec>(tmp);
  }
  
  return theta;}


// update Piic: the probability of each person i in each class c
// input: mu = x\beta; Pi: mixing weights; Y: observed value.
// [[Rcpp::export]]
arma::mat update_piic(const Rcpp::NumericMatrix & mu, const arma::rowvec & Pi,const arma::mat & Y) {
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  unsigned int C = Pi.size();
  Rcpp::NumericMatrix prob(C,J);
  arma::mat Piic(N,J);
  
  for (unsigned int c = 0;c<C;c++){
    prob(c,_)=Rcpp::pnorm(mu(c,_));
  }
  
  arma::mat Prob = as<arma::mat>(prob);
  Piic = Y*log(Prob.t()+0.0000000001)+(1-Y)*log(1-Prob.t()+0.0000000001);
  for (unsigned int n = 0;n<N;n++){
    Piic.row(n)  = Piic.row(n)+log(Pi);
    Piic.row(n) = exp(Piic.row(n))/sum(exp(Piic.row(n)));
  }
  return Piic; }

// update a  according to Piic: class from 0 to P - 1.
// [[Rcpp::export]]
arma::vec update_a(const arma::mat & piic){
  unsigned int P = piic.n_cols;
  unsigned int N = piic.n_rows;
  arma::vec a_class(N);
  arma::vec prob(P);
  
  for (unsigned int i=0;i<N;i++){
    prob = (piic.row(i)).t();
    a_class(i) = rmultinomial(prob);
  }
  return a_class;
}


// update Z  based on mu = X\beta and Y
// [[Rcpp::export]]
arma::mat  update_Z (const arma::mat & mu, const arma::mat & Y) {
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  arma::mat Z(N,J);
  
  for (unsigned int j = 0; j< J; j++){
    Z.col(j)= res_norm(mu.col(j),arma::ones<arma::vec>(N),Y.col(j),arma::zeros<arma::vec>(N));
  }
  return Z;
}




// update mu, sigma  based on Z, X, beta, sigma_B
// [[Rcpp::export]]
Rcpp::List  update_mu_sigma (const arma::mat & Z, const arma::mat & A, const arma::mat & Beta, double sigma_B) {
  unsigned int J = Beta.n_rows;
  unsigned int P = A.n_cols;
  
  arma::mat W(J,P);
  arma::vec Sigma2_p(P);
  arma::mat Mu(J,P);
  
  for (unsigned int p = 0; p < P; p++){
    arma::mat tmp1 = Beta;
    arma::mat tmp2 = A;
    tmp1.shed_col(p);
    tmp2.shed_col(p);
    double sigma2_p = 1.0/(dot(A.col(p),A.col(p))+1.0/(sigma_B*sigma_B));
    Sigma2_p(p) = sigma2_p;
    
    arma::vec mu_p(J);
    
    for (unsigned int j = 0; j < J; j++){
      arma::mat Z_til = Z.col(j) - tmp2*tmp1.row(j).t();
      mu_p(j) = dot(A.col(p),Z_til)*sigma2_p;
    }
    Mu.col(p) = mu_p; 
    
  }
  return Rcpp::List::create(
    Rcpp::Named("Mu")=Mu,
    Rcpp::Named("Sigma2_p")=Sigma2_p
  );}



double truncate(arma::mat terms, const arma::rowvec & beta, unsigned int p){
  unsigned int P = terms.n_cols+1;
  unsigned int K = log(P)/log(2.0);
  
  arma::mat Q = toQ(beta.cols(1,P-1),terms);
  arma::mat O = des_mat(Q,terms);
  arma::mat OO = des_mat(bijectionvector(K,arma::linspace<arma::vec>(0,(P-1),P)),terms);
  arma::rowvec beta_1 = beta;
  beta_1(0) = 0;
  beta_1(p) = 0;
  arma::vec theta = OO*beta_1.t();
  arma::vec candi = -theta(find(OO.col(p)>0));
  double tmp = dot(O.row(0),beta_1.row(0));
  double max_prob = max(theta(find(OO.col(p)<1))) - tmp;
  double B = 0.5*((max(candi)+max_prob)+std::abs(max(candi)-max_prob));
  
  return B;}

// update W 
// [[Rcpp::export]]
double update_W (double mu, double sigma, double B, double w,double sigma_B) {
  
  Rcpp::NumericVector tmp1 = wrap((mu-B)/sigma);
  tmp1 = Rcpp::pnorm(tmp1);
  arma::vec tmp2 = as<arma::vec>(tmp1);
  
  double tmp3 = mu*mu/(2.0*sigma*sigma);
  if(tmp3>50){tmp3=50;}
  
  Rcpp::NumericVector tmp4 = wrap(-B/sigma_B);
  tmp4 = Rcpp::pnorm(tmp4);
  arma::vec tmp5 = as<arma::vec>(tmp4);
  
  double tmp8 = w*sigma/sigma_B*exp(tmp3)*tmp2(0)*(1.0/tmp5(0));
  double W = tmp8/(tmp8+1-w);
  
  
  return W;
}


// [[Rcpp::export]]
Rcpp::List update_D (arma::mat D, const arma::mat & terms, arma::mat Beta, 
                     const arma::mat & Mu, const arma::vec & Sigma_p,
                     double w, double sigma_B){
  // unpdate Delta matrix with strict identifiability conditions
  unsigned int J = D.n_rows;
  unsigned int P = D.n_cols+1;
  arma::mat D_new = D;
  
  Beta.col(0) = myrnorm(Mu.col(0),Sigma_p(0)*arma::ones<arma::vec>(J));
  
  for (unsigned int j=0; j<J;j++){
    arma::rowvec beta = Beta.row(j);
    for (unsigned int p=1; p<P;p++){
      arma::mat D0 = D;
      D0(j,p-1) = 1 - D(j,p-1);
      
      
      double flag = identify_check(D0);
      double b = truncate(terms,beta,p);
      if(D(j,p-1)>0){flag=flag*(b<=0);}
      
      if (flag>0){
        double W = update_W(Mu(j,p), Sigma_p(p), b, w,sigma_B);
        D(j,p-1) = rbernoulli_one(W);
      }
      
      
      
      if(D(j,p-1)>0){
        double beta_p = res_norm1(Mu(j,p), Sigma_p(p),1,b);
        beta(p) = beta_p;
        Beta(j,p) = beta_p;}
      
      if(D(j,p-1)<1){
        beta(p) = 0;
        Beta(j,p) = 0;}
    }
    
  }
  return Rcpp::List::create(Rcpp::Named("Delta")=D,
                            Rcpp::Named("Beta")=Beta);}



// update Delta and beta based on generic identifiability conditions
// [[Rcpp::export]]
Rcpp::List gen_update_D (arma::mat D, const arma::mat & terms, arma::mat Beta, 
                         const arma::mat & Mu, const arma::vec & Sigma_p,
                         double w, double sigma_B, const arma::mat & k3){
  // unpdate Delta matrix
  unsigned int J = D.n_rows;
  unsigned int P = D.n_cols+1;
  
  //update intercepts
  Beta.col(0) = myrnorm(Mu.col(0),Sigma_p(0)*arma::ones<arma::vec>(J));
  
  //update beta row by row 
  for (unsigned int j = 0; j < J;j++){
    // current beta vector (for j-th questions)
    arma::rowvec beta = Beta.row(j);
    for (unsigned int p = 1; p < P;p++){
      
      
      arma::mat D0 = D;
      D0(j,p-1) = 1 - D(j,p-1);
      
      
      
      // check if violates identifiable condtions
      //double flag = 1;
      double flag =  gen_identify_check(D0,k3);
      // obtain the lower bound
      double b = truncate(terms,beta,p);
      // if lower bound > 0, cannot make beta to be 0, i.e., Delta should = 1.
      flag = flag * (b<=0);
      
      if (flag>0){
        // update W
        //double W = update_W(Mu(j,p), Sigma_p(p), -99999, w,sigma_B);
        double W = update_W(Mu(j,p), Sigma_p(p), b, w,sigma_B);
        // update D
        D(j,p-1) = rbernoulli_one(W);
      }
      
      // update beta
      if(D(j,p-1)>0){
        //double beta_p = res_norm1(Mu(j,p), Sigma_p(p),1,-9999);
        double beta_p = res_norm1(Mu(j,p), Sigma_p(p),1,b);
        beta(p) = beta_p;
        Beta(j,p) = beta_p;}
      
      if(D(j,p-1)<1){
        beta(p) = 0;
        Beta(j,p) = 0;}
    }
    
  }
  return Rcpp::List::create(Rcpp::Named("Delta")=D,
                            Rcpp::Named("Beta")=Beta);}



// generate random delta from identifiable space
// [[Rcpp::export]]
arma::mat random_D(unsigned int J,unsigned int K, double w){
  unsigned int P = pow(2,K);
  //Generate two identity matrices
  arma::mat I_K(K,K,arma::fill::eye);
  arma::mat Two_I_K = arma::join_cols(I_K,I_K);
  arma::mat zeros_2K = arma::zeros<arma::mat>((2*K),(P-K-1));
  arma::mat Q0 = arma::join_rows(Two_I_K ,zeros_2K);
  
  //generate Q1
  unsigned int R_Q1 = J-2*K;
  arma::mat U1 = arma::randu<arma::mat>(R_Q1,K);
  arma::mat Q1 = arma::zeros<arma::mat>(R_Q1,K);
  
  //fix elements so columns are nonzero
  arma::vec row_ks = arma::randi<arma::vec>(K,arma::distr_param(0,R_Q1-1));
  for(unsigned int k=0;k<K;k++){
    Q1(row_ks(k),k) = 1;
  }
  
  Q1.elem(arma::find(U1 < w)).fill(1.0);
  
  //Generating the remaining elements of Q in Q2 
  arma::mat U2 = arma::randu<arma::mat>(R_Q1,(P-K-1));
  arma::mat Q2 = arma::zeros<arma::mat>(R_Q1,(P-K-1));
  Q2.elem(arma::find(U2 < w)).fill(1.0);
  
  arma::mat Q = arma::join_rows(Q1,Q2);
  arma::vec r_sum = arma::sum(Q,1);
  arma::uvec r_ind = find(r_sum==0);
  unsigned int r_num = r_ind.size();
  
  //fix elements so rows are nonzero
  if(r_num>0){
    arma::vec col_ks = arma::randi<arma::vec>(r_num,arma::distr_param(0,(P-2)));
    for(unsigned int i=0;i<r_num;i++){
      Q(r_ind(i),col_ks(i)) = 1;
    }
  }
  
  
  Q = arma::join_cols(Q0,Q);
  arma::vec one_J = arma::ones<arma::vec>(J);
  
  //Q
  arma::uvec p = arma::linspace<arma::uvec>(0,(J-1),J);
  for(unsigned int j=0;j<J;j++){
    p(j)=j;
  }
  p = arma::shuffle(p);
  return Q.rows(p);
}






// generate simulation data
// [[Rcpp::export]]
Rcpp::List  Data_gen (unsigned int N, arma::mat Q, double rho,arma::mat terms) {
  unsigned int K = Q.n_cols;
  unsigned int J = Q.n_rows;
  unsigned int P = pow(2,K);
  
  // Theta Matrix 
  arma::mat prob_incrs(J,K);
  arma::vec row_sum = (arma::sum(Q,1));
  arma::rowvec item_j(P);
  
  for (unsigned int j=0;j<J;j++){
    double incres = 0.6/row_sum(j);
    item_j = Q.row(j);
    item_j.elem(find(item_j)).fill(incres);
    prob_incrs.row(j) = item_j;
  }
  
  arma::mat alpha = join_cols(arma::zeros<arma::rowvec>(K) , terms.t()) ;
  arma::mat theta = alpha*prob_incrs.t()+0.2;
  
  
  
  
  // A matrix 
  arma::mat a = arma::zeros<arma::mat>(N,K);
  arma::mat A_sigma = rho*arma::ones<arma::mat>(K,K)+(1-rho)*arma::eye<arma::mat>(K,K);
  // Environment package_env ("package:MASS"); 
  // Function mvrnorm = package_env["mvrnorm"];
  // Rcpp::NumericMatrix a_norm_r = mvrnorm(N,arma::zeros<arma::vec>(K),A_sigma);
  // arma::mat a_norm = as<arma::mat>(a_norm_r);
  // a_norm = R::mvrnorm(N, rep(0), A_sigma);
  arma::mat a_norm = mvrnormArma(N,arma::zeros<arma::vec>(K),A_sigma);
  a.elem(find(a_norm>0)).fill(1.0);
  
  
  //Y matrix
  arma::mat Y(N,J);
  arma::vec theta_i;
  for (unsigned int i=0;i<N;i++){
    theta_i = prob_incrs*a.row(i).t()+0.2;
    Y.row(i) = rbernoulli(theta_i).t();
  }
  return Rcpp::List::create(Rcpp::Named("theta")=theta,
                            Rcpp::Named("a")=a,
                            Rcpp::Named("Y")=Y);}