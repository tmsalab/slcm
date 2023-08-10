#include <RcppArmadillo.h>
// #include <edmcore.h>

template <class T> bool is_needle_in_haystack(T x, unsigned int needle)
{
    arma::uvec m = arma::find(x == needle);

    if (m.n_elem > 0) {
        return true;
    }

    return false;
}

// Attribute Routines ----

//' Generate a vector to map binary vector to integers
//'
//' Converts class into a bijection to integers
//'
//' @param M      Number of Options
//' @param K      Number of Attributes
//'
//' @return
//'
//' Return a matrix containing the class table
//'
//' @noRd
// [[Rcpp::export]]
arma::vec bijectionvector(unsigned int K)
{
    arma::vec vv(K);
    for (unsigned int k = 0; k < K; k++) {
        vv(k) = pow(2, K - k - 1);
    }
    return vv;
}

//' Generate a vector to map polytomous vector to integers
//'
//' Converts class into a bijection to integers
//'
//' @param M      Number of Responses
//' @param K      Number of Attributes
//'
//' @return
//'
//' Return a matrix containing the class table
//'
//' @noRd
template <class T> T gen_bijectionvector(unsigned int K, unsigned int M)
{
    T vv(K);
    for (unsigned int k = 0; k < K; k++) {
        vv(k) = pow(M, K - k - 1.0);
    }
    return vv;
}

//' Create the K inverse bijection of attribute vectors
//'
//' Converts the class into a bijection.
//'
//' @param nClass Number of Attribute Classes
//' @param M      Number of Options
//' @param K      Number of Attributes
//'
//' @return
//'
//' Return a matrix containing the class table
//'
//' @noRd
template <class T>
T inv_gen_bijectionvector(unsigned int K, unsigned int M, double CL)
{
    T alpha(K);
    for (unsigned int k = 0; k < K; k++) {
        double Mpow = pow(M, K - k - 1);
        double ak = 0.;
        while (((ak + 1.) * Mpow <= CL) & (ak < M)) {
            ak += 1.;
        }
        alpha(k) = ak;
        CL = CL - Mpow * alpha(k);
    }
    return alpha;
}

//' Create a K by M^K table of attribute vectors
//'
//' @param nClass Number of Attribute Classes
//' @param M      Number of Responses
//' @param K      Number of Attributes
//'
//' @return
//'
//' Return a matrix containing the class table
//'
//' @noRd
// [[Rcpp::export]]
arma::mat CL_gen_invbijection_table(unsigned int K, unsigned int M,
                                    unsigned int nClass)
{
    arma::mat CLtable(K, nClass);
    for (unsigned int cc = 0; cc < nClass; cc++) {
        CLtable.col(cc) = inv_gen_bijectionvector<arma::vec>(K, M, cc);
    }
    return CLtable;
}

//' Create a permutation table based on indices
//'
//' Generates a vector containing permutation indices.
//'
//' @param nClass Number of Attribute Classes
//' @param K      Number of Attributes
//' @param M      Number of Responses
//' @param order  Order is 1, main-effects, 2 main-effects + interactions.
//'               Highest level of interactions you want.
//' @param vv     Vector of a bijection
//' @param perm   Permutations
//'
//' @return
//' 
//' A `vec` containing indices permutated.
//' 
//' @noRd
// [[Rcpp::export]]
arma::vec permuteAtableIndices(unsigned int nClass, unsigned int K,
                               unsigned int M, unsigned int order,
                               const arma::vec &vv, const arma::vec &perm)
{
    arma::vec vvperm(K);
    arma::vec fullorigperm = arma::linspace(0, nClass - 1, nClass);
    for (unsigned int k = 0; k < K; k++) {
        vvperm(k) = vv(perm(k));
    }
    arma::vec model(nClass);
    arma::vec fullpermindices(nClass);
    for (unsigned int cr = 0; cr < nClass; cr++) {
        arma::vec alpha_r = inv_gen_bijectionvector<arma::vec>(K, M, cr);
        double nof0s = 0.;
        for (unsigned int k = 0; k < K; k++) {
            nof0s += 1. * (alpha_r(k) == 0);
        }
        model(cr) = 1. * (nof0s > double(K - order) - 1);
        arma::vec alpha_perm(K);
        fullpermindices(cr) = arma::accu(alpha_r % vvperm);
    }
    arma::uvec finalcols = find(model == 1);
    arma::vec origperm = fullorigperm(finalcols);
    arma::vec reducedpermindices = fullpermindices(finalcols);
    arma::vec permindices(origperm.n_elem);
    for (unsigned int p = 0; p < origperm.n_elem; p++) {
        double origval = origperm(p);
        for (unsigned int pp = 0; pp < origperm.n_elem; pp++) {
            if (origval == reducedpermindices(pp)) {
                permindices(p) = pp;
            }
        }
    }
    return permindices;
}

// RNG Distribution Routines ----

//' Generate Multinomial Random Variable
//'
//' Sample a multinomial random variable for given probabilities.
//'
//' @param ps A `vector` for the probability of each category.
//' @param M  Number of Categories.
//'
//' @return
//' A `vector` from a multinomial with probability `ps`.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @noRd
// [[Rcpp::export]]
double rmultinomial(const arma::vec &ps, unsigned int M)
{
    double u = R::runif(0, 1);
    arma::vec cps = cumsum(ps);
    arma::vec Ips = arma::zeros<arma::vec>(M);
    Ips.elem(arma::find(cps < u)).fill(1.0);

    return sum(Ips);
}

//' Generate Truncated Normal Random Variable
//'
//' Sample a truncated random variable for given probabilities.
//'
//' @param mean  Mean of the distribution
//' @param sd    Standard deviation of the distribution
//' @param w     Upper bound limits??
//' @param b_lb  Lower bound limits??
//' @param M     Number of Responses??
//'
//' @return
//' A single random `double` value from a truncated normal.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @noRd
// [[Rcpp::export]]
double rTruncNorm_b(double mean, double sd, double w)
{
    // Upper bound?
    double p1 = R::pnorm(mean, 0, sd, 1, 0);
    double p0 = 1 - p1;
    double uZ = R::runif(0, 1);

    // Either upper or lower bound based on whether w is 0 or 1
    double pz = w * (p0 + uZ * p1) + (1 - w) * (uZ * p0);
    double Z = R::qnorm(pz, mean, sd, 1, 0);
    return (Z);
}

// [[Rcpp::export]]
double rTruncNorm_lb(double mean, double sd, double b_lb)
{
    double p0 = R::pnorm(b_lb, mean, sd, 1, 0);
    double p1 = 1 - p0;
    double uZ = R::runif(0, 1);
    double Z = R::qnorm(p0 + uZ * p1, mean, sd, 1, 0);
    return (Z);
}

// [[Rcpp::export]]
double rTruncNorm(double mean, double sd, double w, const arma::vec &ps,
                  unsigned int M)
{
    double uZ = R::runif(0, 1);
    double p0 = ps(w);
    double p1 = ps(w + 1);
    double pz = p0 + uZ * (p1 - p0);
    double Z = R::qnorm(pz, mean, sd, 1, 0);
    return Z;
}

//' Generate Dirichlet Random Variable
//'
//' Sample a Dirichlet random variable.
//'
//' @param deltas A `vector` of Dirichlet parameters.
//'
//' @return
//' A `vector` from a Dirichlet.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @noRd
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec &deltas)
{
    unsigned int C = deltas.n_elem;
    arma::vec Xgamma(C);

    // generating gamma(deltac,1)
    for (unsigned int c = 0; c < C; c++) {
        Xgamma(c) = R::rgamma(deltas(c), 1.0);
    }
    return Xgamma / sum(Xgamma);
}

// generate random J by K Q matrix
// [[Rcpp::export]]
arma::mat random_Q(unsigned int J, unsigned int K)
{
    unsigned int nClass = pow(2, K);
    arma::vec vv = gen_bijectionvector<arma::vec>(K, 2);
    arma::vec Q_biject(J);
    Q_biject(arma::span(0, K - 1)) = vv;
    Q_biject(arma::span(K, 2 * K - 1)) = vv;
    arma::vec Jm2K =
        arma::randi<arma::vec>(J - 2 * K, arma::distr_param(1, nClass - 1));
    Q_biject(arma::span(2 * K, J - 1)) = Jm2K;
    Q_biject = arma::shuffle(Q_biject);
    arma::mat Q(J, K);
    for (unsigned int j = 0; j < J; j++) {
        arma::vec qj = inv_gen_bijectionvector<arma::vec>(K, 2, Q_biject(j));
        Q.row(j) = qj.t();
    }
    return Q;
}

//' Generate tables to store the results during iterations
//'
//' Generate tables to store the results during iterations
//'
//' @param nClass   Number of Attribute Classes
//' @param M        Number of Responses
//' @param K        Number of Attributes
//' @param order    Order is 1, main-effects, 2 main-effects + interactions.
//'               Highest level of interactions you want.
//'
//' @return
//'
//' Return a list containing the tables for different parameters
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List GenerateAtable(unsigned int nClass, unsigned int K, unsigned int M,
                          unsigned int order)
{

    arma::mat FullAtable(nClass, nClass);
    arma::mat FullLBtable(nClass, nClass);
    arma::mat FullDtoQtable(K, nClass);
    arma::mat Fulladjtable(nClass, nClass);
    arma::vec model(nClass);

    // Setup storage for main effect indices
    arma::uvec model_main_effects(nClass);

    // Construct a vv bijection
    arma::uvec vv_bijection = gen_bijectionvector<arma::uvec>(K, M);

    for (unsigned int cr = 0; cr < nClass; cr++) {
        arma::vec alpha_r = inv_gen_bijectionvector<arma::vec>(K, M, cr);

        double nof0s = 0.;
        for (unsigned int k = 0; k < K; k++) {
            nof0s += 1. * (alpha_r(k) == 0);
            FullDtoQtable(k, cr) = 1. * (alpha_r(k) > 0);
        }

        // Checking if the class is in the bijection vector
        // If integer is corresponding with main effect.
        model_main_effects(cr) = 1*is_needle_in_haystack(vv_bijection, cr);

        // Flags coefficients with no greater order
        model(cr) = 1. * (nof0s > double(K - order) - 1);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            arma::vec alpha_c = inv_gen_bijectionvector<arma::vec>(K, M, cc);
            double mindiff = arma::min(alpha_r - alpha_c);
            FullAtable(cr, cc) = 1. * (mindiff > -1);
            double maxdiff = arma::accu(abs(alpha_c - alpha_r));
            FullLBtable(cr, cc) = 1. * (maxdiff == 1) * (mindiff < 0) +
                                  2. * (mindiff > -1) * (maxdiff != 0);
            Fulladjtable(cr, cc) = 1. * (maxdiff == 1);
        }
    }


    // Columns that should be retained under the appropriate model order
    arma::uvec finalcols = find(model == 1);

    // Columns corresponding to the main effects
    // and is present within the model
    arma::uvec maineffectcols = find(model_main_effects == 1 && model == 1);

    arma::mat Atable = FullAtable.cols(finalcols);
    arma::mat LBtable = FullLBtable.cols(finalcols);
    arma::mat DtoQtable = FullDtoQtable.cols(finalcols);
    arma::mat adjtable = Fulladjtable.submat(finalcols, finalcols);

    return Rcpp::List::create(
        Rcpp::Named("Atable", Atable), 
        Rcpp::Named("LBtable", LBtable),
        Rcpp::Named("finalcols", finalcols),
        Rcpp::Named("DtoQtable", DtoQtable), 
        Rcpp::Named("adjtable", adjtable),
        Rcpp::Named("maineffectcols", maineffectcols)
        );
}

// Mapping of Q to Delta
// [[Rcpp::export]]
arma::mat ETAmat(unsigned int K, unsigned int J, unsigned int M,
                 const arma::mat &Q)
{
    double nClass = pow(2, K);
    arma::mat ETA(J, nClass);
    for (unsigned int cc = 0; cc < nClass; cc++) {
        arma::vec alpha_c = inv_gen_bijectionvector<arma::vec>(K, M, cc);
        for (unsigned int j = 0; j < J; j++) {
            arma::rowvec qj = Q.row(j);
            arma::uvec comparisons = qj.t() <= alpha_c;
            double compare = arma::prod(comparisons);
            ETA(j, cc) = compare;
        }
    }
    return ETA;
}

// Translate Beta to Theta
// [[Rcpp::export]]
arma::mat BetatoTheta(unsigned int J, unsigned int nClass,
                      const arma::mat &beta, const arma::mat &Atable)
{
    arma::mat BAp = beta * Atable.t();
    arma::mat theta(J, nClass);
    for (unsigned int j = 0; j < J; j++) {
        for (unsigned int cc = 0; cc < nClass; cc++) {
            theta(j, cc) = R::pnorm(BAp(j, cc), .0, 1., 1, 0);
        }
    }
    return theta;
}

// [[Rcpp::export]]
double computeLB(unsigned int nClass, const arma::mat &LBtable,
                 const arma::mat &Atable, unsigned int p,
                 const arma::rowvec &betaj, double Bp,
                 const arma::uvec &Bindices)
{
    double Bplb = -100.;
    if (p < nClass - 1) {
        arma::vec lbmax(2);
        arma::vec LBp = LBtable.col(p);
        arma::uvec lbinds = find(LBp == 1); // find lower classes
        arma::rowvec ap = Atable.row(Bindices(p));
        arma::rowvec betajpeq0 = betaj;
        betajpeq0(p) = 0.;
        double gamp = arma::accu(ap % betajpeq0);
        arma::vec gams = (Atable.rows(lbinds)) * betaj.t();
        lbmax(0) = arma::max(gams) - gamp;
        arma::uvec lbinds2 = find(LBp == 2); // find greater or equal classes
        arma::vec gamdiff = gamp - (Atable.rows(lbinds2)) * betajpeq0.t();
        lbmax(1) = arma::max(gamdiff);
        Bplb = arma::max(lbmax);
    }
    if (p == nClass - 1) { // need this because there is no 2 in LBtable for
        // last coefficient
        arma::vec LBp = LBtable.col(p);
        arma::uvec lbinds = find(LBp == 1); // find lower classes
        arma::rowvec ap = Atable.row(Bindices(p));
        arma::rowvec betajpeq0 = betaj;
        betajpeq0(p) = 0.;
        double gamp = arma::accu(ap % betajpeq0);
        arma::vec gams = (Atable.rows(lbinds)) * betaj.t();
        Bplb = arma::max(gams) - gamp;
    }
    return Bplb;
}


// [[Rcpp::export]]
double identify_check(const arma::mat Q)
{
    unsigned int K = Q.n_cols;
    unsigned int J = Q.n_rows;

    arma::mat ones_zero_on_diag = -1 * arma::ones<arma::mat>(K, K);
    arma::vec zeros_K = arma::zeros<arma::vec>(K);
    ones_zero_on_diag.diag() = zeros_K;
    arma::vec c_sum = (arma::sum(Q, 0)).t();
    arma::vec r_sum = arma::sum(Q, 1);
    arma::mat I_check = Q * ones_zero_on_diag;
    arma::mat I_count = arma::zeros<arma::mat>(J, K);
    I_count.elem(arma::find(I_check > -1)).fill(1.0);
    arma::vec n_ek = (arma::sum(I_count, 0)).t();

    double min_c = (arma::min(c_sum) > 2);
    double min_r = (arma::min(r_sum) > 0);
    double min_ek = (arma::min(n_ek) > 1);

    return (min_c + min_r + min_ek > 2);
}

// [[Rcpp::export]]
arma::mat eta_dina_matrix(const arma::mat &Q)
{

    // Set up attributes for populating eta.
    unsigned int K = Q.n_cols, J = Q.n_rows;

    // Compute the total number of attribute classes C = 2^K
    double nClass = pow(2, K);

    // Setup storage for the eta
    arma::mat eta_jc(J, nClass);

    for (unsigned int cc = 0; cc < nClass; ++cc) {
        arma::vec alpha_c = inv_gen_bijectionvector<arma::vec>(K, 2, cc);

        for (unsigned int j = 0; j < J; ++j) {
            arma::rowvec qj = Q.row(j);
            // Switch to as_scalar
            double compare = arma::as_scalar(qj * alpha_c - qj * qj.t());
            eta_jc(j, cc) = (compare >= 0);
        }
    }

    return eta_jc;
}

// [[Rcpp::export]]
arma::mat Q_prime_matrix(unsigned int K, const arma::mat &Atable,
                         const arma::vec &vv)
{

    // The number of rows for Q prime must always be p
    // P here is found under the number of columns not the rows of the
    // A table
    // Nuance: Model is not saturated, we may need to update the code elsewhere
    arma::mat Q_prime(Atable.n_cols, K);

    // Setup a basis vector
    arma::vec e_k = arma::zeros<arma::vec>(K);

    for (unsigned int k = 0; k < K; ++k) {

        // Fill the k-th position with 1
        e_k(k) = 1;

        // Perform a dot product
        unsigned int col_ind = arma::dot(e_k, vv);

        // Retrieve from the Atable the appropriate column based off of the
        // matrix bijection vector vv.
        Q_prime.col(k) = Atable.col(col_ind);

        // Reset the k-th position for next iteration
        e_k(k) = 0;
    }

    return Q_prime;
}

// [[Rcpp::export]]
double pnorm_ln_lower_tail(double B_p_lowerbound, double sigma_var_jp)
{
    return R::pnorm(B_p_lowerbound / sqrt(sigma_var_jp), 0.0, 1.0, 1, 1);
}

// [[Rcpp::export]]
double pnorm_ln_upper_tail(double B_p_lowerbound, double sigma_var_jp)
{
    return R::pnorm(B_p_lowerbound / sqrt(sigma_var_jp), 0.0, 1.0, 0, 1);
}

// [[Rcpp::export]]
void update_slipping_guessing(double &slipping, double &guessing,
                              const arma::mat &ab_tilde)
{
    // Global Slipping and Guessing update across parameters

    // Initialize on a uniform
    double slipping_unif = R::runif(0.0, 1.0);
    double guessing_unif = R::runif(0.0, 1.0);

    // Retrieve old slipping value for item j
    double slipping_old = slipping;

    // Draw guessing conditioned upon slipping - 1
    double ab_g1 = ab_tilde(0, 1); // eta 0, delta 1
    double ab_g0 = ab_tilde(0, 0); // eta 0, delta 0

    double pg = R::pbeta(1.0 - slipping_old, ab_g1 + 1., ab_g0 + 1., 1, 0);

    double guessing_new =
        R::qbeta(guessing_unif * pg, ab_g1 + 1., ab_g0 + 1., 1, 0);

    // Draw slipping conditioned upon guessing
    double ab_s1 = ab_tilde(1, 1); // eta_prime 1, delta 1
    double ab_s0 = ab_tilde(1, 0); // eta_prime 1, delta 0

    double ps = R::pbeta(1.0 - guessing_new, ab_s0 + 1., ab_s1 + 1., 1, 0);

    double slipping_new =
        R::qbeta(slipping_unif * ps, ab_s0 + 1., ab_s1 + 1., 1, 0);

    // Update slipping and guessing for all items
    slipping = slipping_new;
    guessing = guessing_new;
}

// [[Rcpp::export]]
double parm_update_nomiss(unsigned int N, unsigned int J, unsigned int K,
                          unsigned int nClass, unsigned int M,
                          const arma::mat &Y, arma::mat &BETA, arma::mat &TAU,
                          arma::vec &CLASS, arma::vec &pis, arma::mat &DELTA,
                          double omega, const arma::vec &vv,
                          const arma::uvec &main_effect_cols,
                          const arma::mat &CLtable, const arma::mat &Atable,
                          const arma::mat &LBtable, const arma::uvec &Bindices,
                          const arma::mat &qtable, unsigned int P,
                          const arma::vec &l1, double m0, 
                          arma::vec &d0, double a0, double b0)
{
    double u;
    arma::mat ApA = arma::zeros<arma::mat>(P, P);
    arma::mat ApZ = arma::zeros<arma::mat>(P, J);

    // Rows are Eta and Columns Delta
    arma::mat ab_tilde = arma::zeros<arma::mat>(2, 2);

    arma::cube PY_a(J, nClass, M + 1);
    PY_a.slice(0) = arma::zeros<arma::mat>(J, nClass);
    PY_a.slice(M) = arma::ones<arma::mat>(J, nClass);
    arma::mat ABETA(J, nClass);
    arma::vec ABETA_sqnorm = arma::zeros<arma::vec>(nClass);

    // cumulative probs Y given alpha
    for (unsigned int cc = 0; cc < nClass; cc++) {
        arma::rowvec a_alpha = Atable.row(cc);
        for (unsigned int j = 0; j < J; j++) {
            double aBj = arma::accu(a_alpha % BETA.row(j));
            ABETA(j, cc) = aBj;
            ABETA_sqnorm(cc) += aBj * aBj;
            for (unsigned int m = 0; m < M - 1; m++) {
                PY_a(j, cc, m + 1) = R::pnorm(TAU(j, m), aBj, 1., 1, 0);
            }
        }
    }

    for (unsigned int i = 0; i < N; i++) {
        arma::rowvec Yi = Y.row(i);
        double class_i = CLASS(i);
        d0(class_i) += -1.;

        // update alphas|B
        arma::vec numerator(nClass);
        double denominator = 0.;
        for (unsigned int cc = 0; cc < nClass; cc++) {
            double picc = 1.;
            for (unsigned int j = 0; j < J; j++) {
                double Yij = Yi(j);
                picc *= (PY_a(j, cc, Yij + 1.) - PY_a(j, cc, Yij));
            }
            numerator(cc) = picc * d0(cc);
            denominator += picc * d0(cc);
        }
        arma::vec pai = numerator / denominator;
        class_i = rmultinomial(pai, nClass);
        CLASS(i) = class_i;
        arma::rowvec a_alpha = Atable.row(class_i);
        d0(class_i) += 1.;

        arma::vec Zi(J);
        // arma::vec ai = CLtable.col(class_i);
        arma::vec aiB = ABETA.col(class_i);

        for (unsigned int j = 0; j < J; j++) {
            double Yij = Yi(j);
            double aiBj = aiB(j);
            arma::vec pYij = PY_a.tube(j, class_i);
            double Zijk = rTruncNorm(aiBj, 1., Yij, pYij, M);
            Zi(j) = Zijk;
        }
        ApA += a_alpha.t() * a_alpha;
        ApZ += a_alpha.t() * Zi.t();
    }

    // update pis
    pis = rDirichlet(d0);

    // update Q, BETA, this version uses MH for q. v2_1 uses Gibbs
    for (unsigned int j = 0; j < J; j++) {
        arma::vec ApZj = ApZ.col(j);
        arma::rowvec betaj = BETA.row(j);

        // ETA 2^K to P matrix
        // Retrieve: 1 x P
        arma::rowvec delta_j = DELTA.row(j);

        // P is just 0 to 2^(K-1) [saturated model]
        arma::uvec idx = arma::linspace<arma::uvec>(0, P - 1, P);

        for (unsigned int p = 0; p < P; p++) {

            // Update the Delta ----

            // Initialize the delta_jp to 1 because the intercept is always
            // held at 1 since we want it in the model.
            double delta_jp = 1;

            // Setup variables to compute the lower bound
            // This is only computed when we are not on the intercept update.
            double B_p = 0;
            double B_p_lowerbound = 0;

            double vforb = l1(p); // permanently active since delta_jp = 1
            double ApAp_p = ApA(p, p);
            double sigma_var_p = 1. / (ApAp_p + vforb);
            double ApA_p = arma::accu(ApA.row(p) % betaj) - ApAp_p * betaj(p);
            double mu_p = sigma_var_p * (ApZj(p) - ApA_p);
            double sigma_sd_p = sqrt(sigma_var_p);

            // Silently the intercept is going to be 1 automatically, so
            // no update will be performed
            if (p > 0) {

                // Compute the lower bound
                B_p = betaj(p);
                B_p_lowerbound =
                    computeLB(nClass, LBtable, Atable, p, betaj, B_p, Bindices);

                if (B_p_lowerbound <= 0) {

                    // Obtain delta_jp values
                    double delta_jp_previous = delta_j(p);
                    double delta_jp_proposed = 1.0 - delta_jp_previous;

                    // Check for generic identifiability
                    bool is_delta_identified = false;

                    // See if the item is in a main effect position
                    // If so, we need to perform an identifiability check on
                    // Delta with only main effects present.
                    if (is_needle_in_haystack(main_effect_cols, p) == true) {

                        // Retrieve entire delta matrix
                        arma::mat DELTA_star = DELTA;
                        DELTA_star(j, p) = delta_jp_proposed;
                        arma::mat DELTA_star_sub = DELTA_star.cols(main_effect_cols);

                        is_delta_identified = identify_check(DELTA_star_sub);

                        // is_delta_generic_identified =
                        //     edmcore::is_q_generic_identified(DELTA_star_sub);

                        // if (is_delta_generic_identified == true) {
                        //     DELTA_star(j, p) = 0;
                        //     DELTA_star_sub = DELTA_star.cols(main_effect_cols);
                        //
                        //     is_delta_generic_identified =
                        //         edmcore::is_q_generic_identified(
                        //             DELTA_star_sub);
                        //
                        //     if (is_delta_generic_identified == false) {
                        //         delta_jp = delta_j(p);
                        //     }
                        } else {
                            delta_jp = delta_j(p);
                        }


                    if (is_delta_identified == true) {

                        //  Obtain the variance
                        double sigma_var_beta = 1.0 / l1(p);

                        // Computes Psi(-L/sigma_b)^(-1) => ln(x) =>
                        // -ln(Psi(-L/sigma_b))
                        double ln_prior_normalizing =
                            -1.0 * pnorm_ln_lower_tail(-1.0 * B_p_lowerbound,
                                                       sigma_var_beta);

                        // Computes Psi((u_jp - L)/sigma_p) => ln(x) => ln((u_jp
                        // - L)/sigma_p)
                        double ln_full_conditional_normalizing =
                            pnorm_ln_lower_tail(
                                (mu_p - B_p_lowerbound) / sigma_sd_p, 1.0);

                        // Computes sigma_p / sigma_b => ln(x) ...
                        double ln_ratio_sd_sigmas =
                            log(sigma_sd_p / sqrt(sigma_var_beta));

                        // exp(1/2 * mu^2 / sigma^2_p)
                        double ln_exp_term = 0.5 * (mu_p * mu_p) / sigma_var_p;

                        double ln_integrated_augmented_likelihood =
                            ln_prior_normalizing + ln_ratio_sd_sigmas +
                            ln_exp_term + ln_full_conditional_normalizing;

                        // Compute the natural log ratios for both components in
                        // delta MH update
                        double ln_ratio_deltas = log(omega) - log(1 - omega);

                        double ln_acceptance_delta_update =
                            (delta_jp_proposed - delta_jp_previous) *
                            (ln_integrated_augmented_likelihood +
                             ln_ratio_deltas);

                        // Perform the Metropolis Hastings Update
                        // Sampling Delta_jp from a truncated normal using MH
                        u = R::runif(0, 1);

                        double flag_delta_update =
                            1. * (ln_acceptance_delta_update > log(u));

                        // Based on MH update change the delta or keep it as-is
                        delta_jp = delta_jp_proposed * flag_delta_update +
                                   (1. - flag_delta_update) * delta_jp_previous;
                        delta_j(p) = delta_jp;
                        DELTA(j, p) = delta_jp;
                    }
                } else {
                    delta_jp = delta_j(p);
                }

            } //

            if (p == 0) {
                mu_p += sigma_var_p * vforb * m0;
                betaj(p) = R::rnorm(mu_p, sigma_sd_p);
            } else {

                // Sample from truncated normal only if delta_jp is 1.
                if (delta_jp == 1) {

                    double Bp = betaj(p);
                    double Bplb = computeLB(nClass, LBtable, Atable, p, betaj,
                                            Bp, Bindices);
                    // betaj(p) = rTruncNorm_b(mub,sqrt(vb),1.);
                    if ((Bplb - mu_p) / sigma_sd_p > 4.264891) {
                        betaj(p) = Bplb;
                    } else {
                        betaj(p) = rTruncNorm_lb(mu_p, sigma_sd_p, Bplb);
                    }

                } else { // indicates delta_jp == 0
                    betaj(p) = 0.0;
                }
            }
        }
        BETA.row(j) = betaj;
    }

    // Add up elements in delta and subtract intercept (e.g. J)
    double sum_delta_without_intercept = arma::accu(DELTA) - DELTA.n_rows;

    double omeganew =
        R::rbeta(sum_delta_without_intercept + 1.0,
                 double((J - 1) * nClass) - sum_delta_without_intercept + 1.0);
    return omeganew;
}

// [[Rcpp::export]]
arma::mat kappa_initialize(unsigned int M, unsigned int J)
{

    arma::mat KAP0(J, M - 1);
    (KAP0.col(0)).fill(.0);
    if (M > 2) {
        for (unsigned int j = 0; j < J; j++) {
            for (unsigned int m = 1; m < M - 1; m++) {
                KAP0(j, m) = KAP0(j, m - 1) + R::runif(1, 2);
            }
        }
    }
    return KAP0;
}

// [[Rcpp::export]]
double RLCMm2ll(unsigned int J, unsigned int N, unsigned int nClass,
                const arma::mat &Y, const arma::mat &theta,
                const arma::vec &pis)
{
    arma::cube THETA(J, nClass, 2);
    THETA.slice(0) = arma::ones<arma::mat>(J, nClass) - theta;
    THETA.slice(1) = theta;
    double m2ll = .0;
    for (unsigned int i = 0; i < N; i++) {
        arma::rowvec Yi = Y.row(i);
        double pYi = 0.;
        for (unsigned int cc = 0; cc < nClass; cc++) {
            double pYi_a = 1.;
            for (unsigned int j = 0; j < J; j++) {
                pYi_a *= THETA(j, cc, Yi(j));
            }
            pYi += pYi_a * pis(cc);
        }
        m2ll += log(pYi);
    }
    return -2. * m2ll;
}

// [[Rcpp::export]]
arma::mat convert_to_q(const arma::mat &structure_matrix,
                       const arma::mat &betas)
{

    // Sum up the values of the relevant coefficients
    // J x K
    arma::mat q_tilde = betas * structure_matrix.t();

    unsigned int J = q_tilde.n_rows, K = q_tilde.n_cols;

    arma::mat q_hat(J, K);

    // Normalize
    for (unsigned k = 0; k < K; ++k) {
        q_hat.col(k) = q_tilde.col(k) / arma::sum(q_tilde, 1);
    }

    // Decision threshold check and conversion to integer matrix
    arma::umat q_decided = (q_hat >= 0.5);

    return arma::conv_to<arma::mat>::from(q_decided);
}

// [[Rcpp::export]]
arma::mat q_to_delta(const arma::mat& Q,
                     const arma::mat& Q_prime,
                     unsigned int M) {

    unsigned int K = Q.n_cols;

    // Table of deltas from Iteration 1 back.
    arma::mat ETA_prime = eta_dina_matrix(Q_prime);

    // Create a bijection vector
    arma::vec vv = gen_bijectionvector<arma::vec>(K, M);

    // J x 2^K (or lower order with J x 2^P)
    arma::mat DELTA(Q.n_rows, Q_prime.n_rows);

    for(unsigned int j = 0; j < Q.n_rows; ++j) {
        // Translate to an integer, pick out column of ETA prime and store it
        // as a row in DELTA
        DELTA.row(j) = ETA_prime.col( dot(Q.row(j), vv) ).t();
    }

    return DELTA;
}


// [[Rcpp::export]]
Rcpp::List slcm_cpp(const arma::mat &Y,
                    unsigned int K, unsigned int M,
                    unsigned int order,
                    const arma::vec &psi_invj, 
                    double m0, double bq,
                    unsigned int burnin = 1000,
                    unsigned int chain_length = 10000)
{
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int nClass = pow(M, K);

    // Chain length
    unsigned int chain_p_burn = chain_length + burnin;

    unsigned int tmburn;
    arma::vec vv = gen_bijectionvector<arma::vec>(K, M);

    // Construct the Atable for the full model
    arma::mat Atable_full = Rcpp::as<arma::mat>(GenerateAtable(nClass, K, M, K)[0]);

    // generating Atable and qtable
    Rcpp::List designmats = GenerateAtable(nClass, K, M, order);
    arma::mat Atable = Rcpp::as<arma::mat>(designmats[0]);
    arma::mat LBtable = Rcpp::as<arma::mat>(designmats[1]);
    arma::uvec Bindices = Rcpp::as<arma::uvec>(designmats[2]);
    arma::mat qtable = Rcpp::as<arma::mat>(designmats[3]);

    // Handles the determination of column information for the model under
    // a specified order's main effects.
    arma::uvec main_effect_cols = Rcpp::as<arma::uvec>(designmats[5]);

    // Ensure order is set correctly
    unsigned int P = Bindices.n_elem;

    // Setup for the delta component
    arma::mat Q_prime = Q_prime_matrix(K, Atable_full, vv);
    Q_prime = Q_prime.rows(Bindices);

    arma::mat CLtable = CL_gen_invbijection_table(K, M, nClass);

    arma::mat TAU = kappa_initialize(M, J);
    arma::mat Q = random_Q(J, K);

    // Relate the randomly sampled Q to DELTA through a realization of the
    // the Q prime/ideal matrix under order M.
    arma::mat DELTA = q_to_delta(Q, Q_prime, M);

    // Sample BETA from DELTA
    arma::mat BETA= arma::randu<arma::mat>(J,P);
    for(unsigned int j=0;j<J;j++){
        arma::rowvec betaj=BETA.row(j);
        arma::rowvec deltaj=DELTA.row(j);
        for(unsigned int p=0;p<P;p++){
            double delta_jp = deltaj(p);
            double vforb = psi_invj(p); // l1s
            double vb=1./vforb;
            double mub=0.;
            mub +=1.*(p==0)*vb*vforb*m0;
            double sdb=sqrt(vb);
            if (p == 0) {
                betaj(p) = R::rnorm(mub, sdb);
            } else {

                // Sample from truncated normal only if delta_jp is 1.
                if (delta_jp == 1) {

                    double Bp = betaj(p);
                    double Bplb = computeLB(nClass, LBtable, Atable, p, betaj,
                                            Bp, Bindices);
                    // betaj(p) = rTruncNorm_b(mub,sqrt(vb),1.);
                    if ((Bplb - mub) / sdb > 4.264891) {
                        betaj(p) = Bplb;
                    } else {
                        betaj(p) = rTruncNorm_lb(mub, sdb, Bplb);
                    }

                } else { // indicates delta_jp == 0
                    betaj(p) = 0.0;
                }
            }
        }
        BETA.row(j)=betaj;
    }
    
    // Initialize omega
    double omega(.5);
    
    // Saving chain output
    arma::cube BETAS(J, P, chain_length);
    arma::cube THETAS(J, nClass, chain_length);
    arma::mat CLs(N, chain_length);
    arma::mat PIs(nClass, chain_length);
    arma::vec omegas(chain_length);

    // Tabulation matrices
    arma::mat DELTA_tab = arma::zeros<arma::mat>(J, P);
    arma::mat Q_tab = arma::zeros<arma::mat>(J, K);

    // Initialize class membership
    arma::vec CLASS =
        arma::randi<arma::vec>(N, arma::distr_param(0, nClass - 1));
    double a0(1.), b0(1.);
    
    arma::vec d0 = arma::ones<arma::vec>(nClass);
    arma::vec pis = rDirichlet(d0);
    
    for (unsigned int i = 0; i < N; ++i) {
        double class_i = CLASS(i);
        d0(class_i) += 1.;
    }
    
    // Store -2LL value for each chain iteration.
    arma::vec m2LLs(chain_length);

    // Start Markov chain
    for (unsigned int t = 0; t < chain_p_burn; t++) {
      
        // Obtain an update for all parameters. 
        // Save the omega result back into the original 
        omega = parm_update_nomiss(N, J, K, nClass, M, Y, BETA, TAU, CLASS, pis,
                                   DELTA, omega, vv, main_effect_cols, CLtable,
                                   Atable, LBtable, Bindices, qtable, P,
                                   psi_invj, m0, d0, a0, b0);

        if (t > burnin - 1) {
            tmburn = t - burnin;
            arma::mat theta = BetatoTheta(J, nClass, BETA, Atable);
            // Compute -2ll for DIC criteria
            double m2ll = RLCMm2ll(J, N, nClass, Y, theta, pis);
            BETAS.slice(tmburn) = BETA;
            CLs.col(tmburn) = CLASS;
            THETAS.slice(tmburn) = theta;
            PIs.col(tmburn) = pis;
            DELTA_tab += DELTA;
            Q_tab += BETA * qtable.t();
            omegas(tmburn) = omega;
            m2LLs(tmburn) = m2ll;
        }
    }

    // Storage for obtaining an estimate for Q based on Delta
    arma::mat q_hat(J, K);

    // Use the heuristic to calculate the Q estimate. See SLCM paper.
    for (unsigned k = 0; k < K; ++k) {
        q_hat.col(k) = Q_tab.col(k) / arma::max(Q_tab, 1);
    }

    // Decision threshold check and conversion to integer matrix
    // True -> 1, False -> 0
    arma::umat q_decided = (q_hat >= 0.5);
    
    // Dichotomized Q estimate matrix as 1.0/0.0
    arma::mat Q_transformed = arma::conv_to<arma::mat>::from(q_decided);

    // Export estimates
    Rcpp::List estimates =
        Rcpp::List::create(
            Rcpp::Named("delta", DELTA_tab / (chain_length)),
            Rcpp::Named("q", Q_transformed),
            Rcpp::Named("tau", TAU)
        );

    // Export chain diagnostics
    Rcpp::List chain =
        Rcpp::List::create(
            Rcpp::Named("beta", BETAS),
            Rcpp::Named("theta", THETAS),
            Rcpp::Named("class", CLs),
            Rcpp::Named("pi", PIs),
            Rcpp::Named("omega", omegas),
            Rcpp::Named("m2ll", m2LLs)
        );

    // Return model object
    return Rcpp::List::create(
        Rcpp::Named("estimates", estimates),
        Rcpp::Named("chain", chain)
    );
}
